# -*- coding: utf-8 -*-

import logging
import os

import nibabel as nib
import numpy as np
import pandas as pd

from .datasets import (check_donors, fetch_fsaverage5, fetch_fsnative,
                       WELL_KNOWN_IDS)
from .samples_ import ONTOLOGY
from .utils import labeltable_to_df, load_gifti, first_entry
from . import matching, transforms

LGR = logging.getLogger('abagen')
BACKGROUND = [
    'unknown', 'corpuscallosum',
    'Background+FreeSurfer_Defined_Medial_Wall', '???'
]


def _reorient_image(image, orientation='RAS'):
    """
    Re-orients `image` to `orientation`

    Parameters
    ----------
    image : niimg_like
        Image to be re-oriented
    orientation : str or tuple-of-str
        Orientation, drawing from options ('L', 'R')('I', 'S')('P', 'S').
        Default: 'RAS'

    Returns
    -------
    reoriented : niimg_like
        Re-oriented image
    """

    orig_ornt = nib.io_orientation(image.affine)
    targ_ornt = nib.orientations.axcodes2ornt(orientation)
    transform = nib.orientations.ornt_transform(orig_ornt, targ_ornt)
    image = image.as_reoriented(transform)

    return image


def leftify_atlas(atlas):
    """
    Zeroes out all ROIs in the right hemisphere of volumetric `atlas`

    Assumes that positive X values indicate the right hemisphere (e.g., RAS+
    orientation) and that the X-origin is in the middle of the brain

    Parameters
    ----------
    atlas : str or niimg-like
        Filepath to or in-memory loaded image

    Returns
    -------
    atlas : niimg-like
        Loaded image with right hemisphere zeroed out
    """

    atlas = check_img(atlas)

    # get original orientation then orient to RAS+
    orient = nib.orientations.ornt2axcodes(nib.io_orientation(atlas.affine))
    atlas = _reorient_image(atlas, 'RAS')

    # get ijk corresponding to zero-point
    i = np.squeeze(transforms.xyz_to_ijk([0, 0, 0], atlas.affine))[0]

    # zero out all positive voxels; img is RAS+ so positive = right hemisphere
    data = np.array(atlas.dataobj, copy=True)
    data[int(i):] = 0

    out = atlas.__class__(data, atlas.affine, header=atlas.header)

    return _reorient_image(out, orient)


def annot_to_gifti(atlas):
    """
    Converts FreeSurfer-style annotation file `atlas` to in-memory GIFTI image

    Parameters
    ----------
    annot : os.PathLike
        Surface annotation file (.annot)

    Returns
    -------
    gifti : nib.gifti.GiftiImage
        Converted gifti image
    """

    labels, ctab, names = nib.freesurfer.read_annot(atlas)

    darr = nib.gifti.GiftiDataArray(labels, intent='NIFTI_INTENT_LABEL',
                                    datatype='NIFTI_TYPE_INT32')
    labeltable = nib.gifti.GiftiLabelTable()
    for key, label in enumerate(names):
        (r, g, b), a = (ctab[key, :3] / 255), (1.0 if key != 0 else 0.0)
        glabel = nib.gifti.GiftiLabel(key, r, g, b, a)
        glabel.label = label.decode()
        labeltable.labels.append(glabel)

    return nib.GiftiImage(darrays=[darr], labeltable=labeltable)


def _relabel(labels, minval=0, bgval=None):
    """
    Relabels `labels` so that they're consecutive

    Parameters
    ----------
    labels : (N,) array_like
        Labels to be re-labelled
    minval : int, optional
        What the new minimum value of the labels should be. Default: 0
    bgval : int, optional
        What the background value should be; the new labels will start at
        `minval` but the first value of these labels (i.e., labels == `minval`)
        will be set to `bgval`. Default: None

    Returns
    ------
    labels : (N,) np.ndarray
        New labels
    """

    labels = np.unique(labels, return_inverse=True)[-1] + minval
    if bgval is not None:
        labels[labels == minval] = bgval
    return labels


def relabel_gifti(atlas, background=BACKGROUND, offset=None):
    """
    Updates GIFTI images so label IDs are consecutive across hemispheres

    Parameters
    ----------
    atlas : (2,) tuple-of-str
        Surface label files in GIFTI format (lh.label.gii, rh.label.gii)
    background : list-of-str, optional
        If provided, a list of IDs in `atlas` that should be set to 0 (the
        presumptive background value). Other IDs will be shifted so they are
        consecutive (i.e., 0--N). Default: `abagen.images.BACKGROUND`
    offset : int, optional
        What the lowest value in `atlas[1]` should be not including background
        value. If not specified it will be purely consecutive from `atlas[0]`.
        Default: None

    Returns
    -------
    relabelled : (2,) tuple-of-nib.gifti.GiftiImage
        Re-labelled `atlas` files
    """

    out = tuple()
    minval = 0
    for hemi in atlas:
        # get necessary info from file
        img = load_gifti(hemi)
        data = img.agg_data()
        labels = img.labeltable.labels
        lt = {v: k for k, v in img.labeltable.get_labels_as_dict().items()}

        # get rid of labels we want to drop
        if background is not None:
            for val in background:
                idx = lt.get(val, 0)
                if idx == 0:
                    continue
                data[data == idx] = 0
                labels = [f for f in labels if f.key != idx]

        # reset labels so they're consecutive and update label keys
        data = _relabel(data, minval=minval, bgval=0)
        ids = np.unique(data)
        for n, i in enumerate(ids):
            labels[n].key = i
        minval = len(ids) - 1 if offset is None else int(offset) - 1

        # make new gifti image with updated information
        darr = nib.gifti.GiftiDataArray(data, intent='NIFTI_INTENT_LABEL',
                                        datatype='NIFTI_TYPE_INT32')
        labeltable = nib.gifti.GiftiLabelTable()
        labeltable.labels = labels
        img = nib.GiftiImage(darrays=[darr], labeltable=labeltable)
        out += (img,)

    return out


def check_img(img):
    """
    Very basic checker that loads `img`` and ensures it's 3D/int

    Parameters
    --------
    img : str or niimg-like
        Filepath to or in-memory loaded image

    Returns
    -------
    img : niimg-like
        Loaded 3D/int image
    """

    if isinstance(img, (str, os.PathLike)) and os.path.exists(img):
        img = nib.load(img)
    elif not isinstance(img, nib.spatialimages.SpatialImage):
        raise TypeError('Provided image must be an existing filepath or a '
                        'pre-loaded niimg-like object')

    # ensure 3D or squeezable to 3D
    img = nib.funcs.squeeze_image(img)
    if len(img.shape) != 3:
        raise ValueError('Provided image must be 3D')

    # check if atlas data is int or castable to int
    # if image is arrayproxy convert it to an array for speed-up
    data = np.asarray(img.dataobj)
    cast = nib.is_proxy(img.dataobj)
    if img.header.get_data_dtype().kind not in ['i', 'u']:
        idata = data.astype('int32')
        cast = np.allclose(idata, data)
        data = idata
        if not cast:
            raise ValueError('Provided image should have integer values or '
                             'be safely castable to int without data loss')
    if cast:
        img = img.__class__(data, img.affine, header=img.header)
        img.header.set_data_dtype(np.int32)

    return img


def check_surface(atlas):
    """
    Very basic checker that loads provided surface `atlas`

    Parameters
    ----------
    atlas : (2,) tuple
        Surface files in GIFTI format (lh, rh)

    Returns
    -------
    atlas : (N,) np.ndarray
        Loaded parcellation information from provided `atlas`
    atlas_info : pd.DataFrame or None
        If provided `atlas` files have valid GIFTI label tables then this will
        be a dataframe with information about the loaded parcellation; if not,
        this value will be None
    """

    # if we're not dealing with a len-2 tuple of GIFTI, check if it's a single
    # string (and barf) or just return it (hoping that it's array_like)
    if len(atlas) != 2:
        if isinstance(atlas, (str, os.PathLike)):
            raise TypeError('Must provide a tuple of surface atlases')
        return atlas, None
    for img in atlas:
        if (not isinstance(img, nib.GiftiImage)
                and not (img.endswith('.gii') or img.endswith('.gii.gz'))):
            raise TypeError('Provided surface atlases must be in GIFTI format')

    adata, labs = [], []
    for hemi in atlas:
        hemi = load_gifti(hemi)
        data = np.squeeze(hemi.agg_data())
        if data.ndim > 1:
            raise ValueError('Provided GIFTIs must have only a single, vector '
                             'data array')
        # if data aren't integer then they might not be a parcellation.
        # check that they're safely castable and, if so, use the casted int
        if data.dtype.char not in ('i'):
            idata = data.astype('int32')
            if np.allclose(idata, data):
                data = idata
            else:
                raise ValueError('Provided GIFTIs do not seem to be valid '
                                 'label.gii files')
        adata.append(data)
        ldict = hemi.labeltable.get_labels_as_dict()
        labs.append({k: ldict.get(k) for k in np.unique(data)})

    # we need each hemisphere to have unique values so they don't get averaged
    # check to see if the two hemispheres have more than 1 overlapping value
    # (assume exactly one for a background value of 0 for e.g., medial wall)
    offset = len(np.intersect1d(*adata))
    if offset > 1:
        offset = len(np.unique(adata))
        adata[1] += offset
        labs[1] = {k + offset: v for k, v in labs[1].items()}

    adata = np.hstack(adata)
    atlas_info = labeltable_to_df(labs)

    return adata, atlas_info


def check_atlas(atlas, atlas_info=None, geometry=None, space=None, donor=None,
                data_dir=None):
    """
    Checks that `atlas` is a valid atlas

    Parameters
    ----------
    atlas : niimg-like object or (2,) tuple-of-GIFTI
        Parcellation image or surface, where voxels / vertices belonging to a
        given parcel are identified with a unique integer ID
    atlas_info : {os.PathLike, pandas.DataFrame, None}, optional
        Filepath or dataframe containing information about `atlas`. Must have
        at least columns ['id', 'hemisphere', 'structure'] containing
        information mapping `atlas` IDs to hemisphere (i.e., "L", "R", "B") and
        broad structural class (i.e.., "cortex", "subcortex/brainstem",
        "cerebellum", "white matter", or "other"). Default: None
    geometry : (2,) tuple-of-GIFTI, optional
        Surfaces files defining geometry of `atlas`, if `atlas` is a tuple of
        GIFTI images. Default: None
    space : {'fsaverage', 'fsnative', 'fslr'}, optional
        If `geometry` is supplied, what space files are in. Default: None
    donor : str, optional
        If specified, indicates which donor the specified `atlas` belongs to.
        Only relevant when `atlas` is surface-based, to ensure the correct
        geometry files are fetched. Default: None (i.e., group-level atlas)
    data_dir : str, optional
        Directory where donor-specific FreeSurfer data should be downloaded and
        unpacked. Only used if provided `donor` is not None. Default: $HOME/
        abagen-data

    Returns
    -------
    atlas : :obj:`abagen.AtlasTree`
        AtlasTree object with information about `atlas` and functionality for
        labelling coordinates
    """

    if isinstance(atlas, matching.AtlasTree):
        if atlas_info is not None:
            atlas.atlas_info = atlas_info
        return atlas

    try:
        atlas = check_img(atlas)
        coords = triangles = None
    except TypeError:
        atlas, info = check_surface(atlas)
        # backwards compatibility for `donor` keyword
        if geometry is None and donor is None:
            geometry = fetch_fsaverage5()
            space = 'fsaverage5'
        elif geometry is None and donor is not None:
            geometry = fetch_fsnative(donor, data_dir=data_dir)
            space = 'fsnative'
        elif geometry is not None and space is None:
            raise ValueError('If providing geometry files space parameter '
                             'must be specified')
        coords, triangles = check_geometry(geometry, space, donor=donor,
                                           data_dir=data_dir)
        if atlas_info is None and info is not None:
            atlas_info = info

    atlas = matching.AtlasTree(atlas, coords=coords, triangles=triangles)

    if atlas_info is not None:
        atlas.atlas_info = atlas_info

    return atlas


def check_geometry(surface, space, donor=None, data_dir=None):
    """
    Loads geometry `surface` files and transforms coordinates in `space`

    Parameters
    ----------
    surface : (2,) tuple-of-GIFTI
        Surface geometry files in GIFTI format (lh, rh)
    space : {'fsaverage', 'fsnative', 'fslr'}
        What space `surface` files are in; used to apply appropriate transform
        to MNI152 space. If 'fsnative' then `donor` must be supplied as well
    donor : str, optional
        If specified, indicates which donor the specified `surface` belongs to
    data_dir : str, optional
        Directory where donor-specific FreeSurfer data exists (or should be
        downloaded and unpacked). Only used if provided `donor` is not None.
        Default: $HOME/abagen-data

    Returns
    -------
    coords : (N, 3) np.ndarray
        Coordinates from `surface` files
    triangles : (T, 3) np.ndarray
        Triangles from `surface` files
    """

    if len(surface) != 2:
        raise TypeError('Must provide a tuple of geometry files')

    # fsaverage5, fsaverage6, etc
    if 'fsaverage' in space and space != 'fsaverage':
        space = 'fsaverage'
    space_opts = ('fsaverage', 'fsnative', 'fslr')
    if space not in space_opts:
        raise ValueError(f'Provided space must be one of {space_opts}.')
    if space == 'fsnative' and donor is None:
        raise ValueError('Specified space is "fsnative" but no donor ID '
                         'supplied')

    try:
        coords, triangles = map(list, zip(*[
            load_gifti(img).agg_data() for img in surface
        ]))
    except TypeError:
        coords, triangles = map(list, zip(*[i for i in surface]))

    triangles[-1] += coords[0].shape[0]
    coords, triangles = np.row_stack(coords), np.row_stack(triangles)

    if space == 'fsaverage':
        coords = transforms.fsaverage_to_mni152(coords)
    elif space == 'fsnative':
        coords = transforms.fsnative_to_xyz(coords, donor, data_dir=data_dir)

    return coords, triangles


def check_atlas_info(atlas_info, labels):
    """
    Checks whether provided `atlas_info` is correct format for processing

    Parameters
    ----------
    atlas_info : str or pandas.DataFrame
        Filepath or dataframe containing information about atlas. Must have
        at least columns 'id', 'hemisphere', and 'structure' containing
        information mapping atlas IDs to hemisphere (i.e., "L", "R", "B") and
        broad structural class (i.e.., "cortex", "subcortex/brainstem",
        "cerebellum", "white matter", or "other").
    labels : array_like
        List of parcel IDs that should be present in `atlas_info`

    Returns
    -------
    atlas_info : pandas.DataFrame
        Loaded dataframe with information on atlas
    """

    valid_structures = list(ONTOLOGY.value_set('structure'))
    hemi_swap = {
        'lh': 'L', 'LH': 'L', 'l': 'L', 'left': 'L',
        'rh': 'R', 'RH': 'R', 'r': 'R', 'right': 'R',
        'bilateral': 'B'
    }
    struct_swap = {
        'subcortex': 'subcortex/brainstem', 'brainstem': 'subcortex/brainstem'
    }
    expected_cols = ['hemisphere', 'structure']

    # load info, if not already
    if not isinstance(atlas_info, pd.DataFrame):
        try:
            atlas_info = pd.read_csv(atlas_info)
        except (ValueError, TypeError):
            pass

    try:
        atlas_info = atlas_info.copy()
        if 'id' in atlas_info.columns:
            atlas_info = atlas_info.set_index('id')
    except AttributeError:
        raise TypeError('Provided `atlas_info` must be a filepath or pandas.'
                        'DataFrame. Please confirm inputs and try again.')

    try:
        assert all(c in atlas_info.columns for c in expected_cols)
        assert 'id' == atlas_info.index.name
        assert len(np.setdiff1d(labels, atlas_info.index)) == 0
    except AssertionError:
        raise ValueError('Provided `atlas_info` does not have adequate '
                         'information on all `labels`. Please confirm '
                         'that atlas_info has columns [\'id\', '
                         '\'hemisphere\', \'structure\'], and that the region '
                         'IDs listed in `atlas_info` account for all those '
                         'found in atlas.')

    try:
        atlas_info['hemisphere'] = atlas_info['hemisphere'].replace(hemi_swap)
        hemi_diff = np.setdiff1d(atlas_info['hemisphere'], ['L', 'R', 'B'])
        assert len(hemi_diff) == 0
    except AssertionError:
        raise ValueError('Provided `atlas_info` has invalid values in the'
                         '\'hemisphere\' column. Only the following values '
                         'are allowed: {}. Invalid value(s): {}'
                         .format(['L', 'R', 'B'], hemi_diff))

    try:
        atlas_info['structure'] = atlas_info['structure'].replace(struct_swap)
        struct_diff = np.setdiff1d(atlas_info['structure'], valid_structures)
        assert len(struct_diff) == 0
    except AssertionError:
        raise ValueError('Provided `atlas_info` has invalid values in the'
                         '\'structure\' column. Only the following values are '
                         'allowed: {}. Invalid value(s): {}'
                         .format(valid_structures, struct_diff))

    return atlas_info


def coerce_atlas_to_dict(atlas, donors, atlas_info=None, data_dir=None):
    """
    Coerces `atlas` to dict with keys `donors`

    If already a dictionary, confirms that `atlas` has entries for all values
    in `donors`

    Parameters
    ----------
    atlas : niimg-like object
        A parcellation image in MNI space, where each parcel is identified by a
        unique integer ID
    donors : array_like
        Donors that should have entries in returned `atlas` dictionary
    atlas_info : os.PathLike or pandas.DataFrame, optional
        Filepath to or pre-loaded dataframe containing information about
        `atlas`. Must have at least columns 'id', 'hemisphere', and 'structure'
        containing information mapping atlas IDs to hemisphere (i.e, "L", "R")
        and broad structural class (i.e., "cortex", "subcortex/brainstem",
        "cerebellum"). If provided, this will constrain matching of tissue
        samples to regions in `atlas`. Default: None
    data_dir : str, optional
        Directory where data should be downloaded and unpacked. Only used if
        provided `atlas` is a dictionary of surface files. Default: $HOME/
        abagen-data

    Returns
    -------
    atlas : dict
        Dict where keys are `donors` and values are `atlas`. If a dict was
        provided it is checked to ensure
    group_atlas : bool
        Whether one atlas was provided for all donors (True) instead of
        donor-specific atlases (False)
    """

    donors = check_donors(donors)
    group_atlas = True

    try:
        atlas = {
            WELL_KNOWN_IDS.subj[donor]: check_atlas(atl, atlas_info,
                                                    donor=donor,
                                                    data_dir=data_dir)
            for donor, atl in atlas.items()
        }
        # if it's a group atlas they should all be the same object
        base = first_entry(atlas)
        group_atlas = all(base is atl for atl in atlas.values())
        missing = set(donors) - set(atlas)
        if len(missing) > 0:
            raise ValueError('Provided `atlas` does not have entry for all '
                             f'requested donors. Missing donors: {donors}.')
    except AttributeError:
        atlas = check_atlas(atlas, atlas_info)
        atlas = {donor: atlas for donor in donors}

    if group_atlas:
        LGR.info('Group-level atlas provided; using MNI coords for '
                 'tissue samples')
    else:
        LGR.info('Donor-specific atlases provided; using native coords for '
                 'tissue samples')

    # update group atlas status based on what was decided / derived
    for atl in atlas.values():
        atl.group_atlas = group_atlas

    return atlas, group_atlas
