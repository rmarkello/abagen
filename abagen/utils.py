# -*- coding: utf-8 -*-

import os
import itertools

import nibabel as nib
import numpy as np
import pandas as pd
from scipy.ndimage import center_of_mass
from scipy.spatial.distance import cdist
from scipy.stats import zscore

AGG_FUNCS = dict(
    mean=np.mean,
    median=np.median
)


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
    if isinstance(img, str) and os.path.exists(img):
        img = nib.load(img)
    elif not isinstance(img, nib.spatialimages.SpatialImage):
        raise TypeError('Provided image must be an existing filepath or a '
                        'pre-loaded niimg-like object')

    # ensure 3D or squeezable to 3D (this is a faster check than type casting)
    if len(img.shape) == 4 and img.shape[3] == 1:
        data = np.asarray(img.dataobj)
        affine = img.affine
        img = img.__class__(data[:, :, :, 0], affine, header=img.header)
    elif len(img.shape) != 3:
        raise ValueError('Provided image must be 3D')

    # check if atlas is int or castable to int
    if np.asarray(img.dataobj).dtype.kind != 'i':
        cast = all([float(f).is_integer() for f in np.unique(img.dataobj)])
        if not cast:
            raise ValueError('Provided image should have integer values or '
                             'be safely castable to integer')
        img = img.__class__(np.asarray(img.dataobj).astype(np.int32),
                            img.affine, header=img.header)
        img.header.set_data_dtype(np.int32)

    return img


def check_atlas_info(atlas, atlas_info, labels=None, validate=False):
    """
    Checks whether provided `info` on `atlas` is sufficient for processing

    Parameters
    ----------
    atlas : niimg-like object
        Parcellation image, where voxels belonging to a given parcel should be
        identified with a unique integer ID
    atlas_info : str or pandas.DataFrame
        Filepath or dataframe containing information about `atlas`. Must have
        at least columns 'id', 'hemisphere', and 'structure' containing
        information mapping atlas IDs to hemisphere (i.e., "L" or "R") and
        broad structural class (i.e.., "cortex", "subcortex", "cerebellum",
        "brainstem", "white matter", or "other").
    labels : array_like, optional
        List of values containing labels to compare between `atlas` and
        `atlas_info`, if they don't all match. If not specified this function
        will attempt to confirm that all IDs present in `atlas` have entries in
        `atlas_info` and vice versa. Default: None
    validate : bool, optional
        Whether to only validate (True) the provided `atlas` and `atlas_info`
        instead of returning (False) the validated dataframe. Default: False

    Returns
    -------
    atlas_info : pandas.DataFrame
        Loaded dataframe with information on atlas
    """

    from .samples import ONTOLOGY  # avoid circular imports

    atlas = check_img(atlas)
    ids = get_unique_labels(atlas) if labels is None else labels
    valid_structures = list(ONTOLOGY.value_set('structure'))
    hemi_swap = {
        'lh': 'L', 'LH': 'L', 'l': 'L',
        'rh': 'R', 'RH': 'R', 'r': 'R'
    }
    expected_cols = ['hemisphere', 'structure']

    # load info, if not already
    try:
        atlas_info = pd.read_csv(atlas_info)
    except ValueError:
        pass

    try:
        atlas_info = atlas_info.copy()
        if 'id' in atlas_info.columns:
            atlas_info = atlas_info.set_index('id')
    except AttributeError:
        raise TypeError('Provided atlas_info must be a filepath or pandas.'
                        'DataFrame. Please confirm inputs and try again.')

    try:
        assert all(c in atlas_info.columns for c in expected_cols)
        assert 'id' == atlas_info.index.name
        assert len(np.setdiff1d(ids, atlas_info.index)) == 0
    except AssertionError:
        raise ValueError('Provided atlas_info does not have adequate '
                         'information on supplied atlas. Please confirm '
                         'that atlas_info has columns [\'id\', '
                         '\'hemisphere\', \'structure\'], and that the region '
                         'IDs listed in atlas_info account for all those '
                         'found in atlas.')

    try:
        atlas_info['hemisphere'] = atlas_info['hemisphere'].replace(hemi_swap)
        hemi_diff = np.setdiff1d(atlas_info['hemisphere'], ['L', 'R'])
        assert len(hemi_diff) == 0
    except AssertionError:
        raise ValueError('Provided atlas_info has invalid values in the'
                         '\'hemisphere\' column. Only the following values '
                         'are allowed: {}. Invalid value(s): {}'
                         .format(['L', 'R'], hemi_diff))

    try:
        struct_diff = np.setdiff1d(atlas_info['structure'], valid_structures)
        assert len(struct_diff) == 0
    except AssertionError:
        raise ValueError('Provided atlas_info has invalid values in the'
                         '\'structure\' column. Only the following values are '
                         'allowed: {}. Invalid value(s): {}'
                         .format(valid_structures, struct_diff))

    if not validate:
        return atlas_info


def check_metric(metric):
    """
    Confirms `metric` is a valid aggregation metric

    Parameters
    ----------
    metric : {'mean', 'median'} or callable, optional
        Mechanism by which to reduce arrays. If a callable, should be able to
        accept an `N`-dimensional input and the `axis` keyword argument and
        return an `N-1`-dimensional output. Default: 'mean'

    Returns
    -------
    metric : callable
    """

    if not callable(metric):
        try:
            metric = AGG_FUNCS[metric]
        except KeyError:
            raise ValueError('Provided metric {0} is not valid. If supplied'
                             'as string, metric must be in {1}.'
                             .format(metric, list(AGG_FUNCS)))
    else:
        try:
            test = np.random.rand(10, 10, 10)
            assert metric(test, axis=0).shape == (10, 10)
            assert isinstance(metric(test), float)
        except (AssertionError, TypeError):
            raise TypeError('Provided metric {0} does not perform as '
                            'expected. Please ensure it accepts the `axis` '
                            'keyword argument and can reduce an `N`-'
                            'dimensional input to an `N-1`-dimensional output.'
                            .format(metric))

    return metric


def efficient_corr(x, y):
    """
    Computes correlation of matching columns in `x` and `y`

    Parameters
    ----------
    x, y : (N, M) array_like
        Input data arrays

    Returns
    -------
    corr : (M,) numpy.ndarray
        Correlations of columns in `x` and `y`
    """
    corr = np.sum(zscore(x, ddof=1) * zscore(y, ddof=1), axis=0) / (len(x) - 1)

    return corr


def get_unique_labels(label_image):
    """
    Returns all possible ROI labels from ``label_image``

    Parameters
    ----------
    label_image : niimg-like object
        ROI image, where each ROI is identified with a unique integer ID

    Returns
    -------
    labels : np.ndarray
        Integer labels of all ROIS found within ``label_image``
    """

    label_image = check_img(label_image)
    return np.trim_zeros(np.unique(label_image.dataobj)).astype(int)


def get_centroids(image, labels=None, image_space=False):
    """
    Finds centroids of ``labels`` in ``label_image``

    Parameters
    ----------
    label_image : niimg-like object
        3D image containing integer label at each point
    labels : array_like, optional
        List of values containing labels of which to find centroids. Default:
        all possible labels
    image_space : bool, optional
        Whether to return xyz (image space) coordinates for centroids based
        on transformation in ``label_image.affine``. Default: False

    Returns
    -------
    centroids : (N, 3) np.ndarray
        Coordinates of centroids for ROIs in input data
    """

    image = check_img(image)
    data = np.asarray(image.dataobj)

    # if no labels of interest provided, get all possible labels
    if labels is None:
        labels = np.trim_zeros(np.unique(data))

    # get centroids for all possible labels
    centroids = np.row_stack(center_of_mass(data, labels=data, index=labels))

    # return xyz if desired; otherwise, ijk
    if image_space:
        centroids = ijk_to_xyz(centroids, image.affine)

    return centroids


def closest_centroid(coord, centroids, return_dist=False):
    """
    Returns index of `centroids` closest to `coord`

    Computes Euclidean distances between `coord` and each of `centroids`,
    returning index of closest centroid

    Parameters
    ----------
    coord : (1, 3) array_like
        Coordinates of sample
    centroids : (N, 3) array_like
        Centroids of parcels (in same space as `coord`)
    return_dist : bool, optional
        Whether to also return distance of closest centroid

    Returns
    -------
    closest : int
        Index of closest centroid in `centroids` to `coord`
    distance : float
        Distance of closest centroid in `centroids` to `coord`. Only returned
        if `return_dist=True`
    """

    distances = np.squeeze(cdist(np.atleast_2d(coord), centroids))
    closest = distances.argmin(axis=0)

    if return_dist:
        return closest, distances[closest]

    return closest


def _check_coord_inputs(coords):
    """
    Confirms `coords` are appropriate shape for coordinate transform

    Parameters
    ----------
    coords : array-like

    Returns
    -------
    coords : (3, N) np.ndarray
    """
    coords = np.atleast_2d(coords).T
    if 3 not in coords.shape:
        raise ValueError('Input coordinates must be of shape (3 x N). '
                         'Provided coordinate shape: {}'.format(coords.shape))
    if coords.shape[0] != 3:
        coords = coords.T
    # add constant term to coords to make 4 x N
    coords = np.row_stack([coords, np.ones_like(coords[0])])
    return coords


def ijk_to_xyz(coords, affine):
    """
    Converts `coords` in cartesian space to `affine` space

    Parameters
    ----------
    coords : (N, 3) array_like
        Cartesian (ijk) coordinate values
    affine : (4, 4) array_like
        Affine matrix containing displacement + boundary

    Returns
    -------
    xyz : (N, 3) np.ndarray
        Provided ``coords`` in ``affine`` space
    """
    coords = _check_coord_inputs(coords)
    aff_coords = np.dot(affine, coords)[:3].T
    return aff_coords


def xyz_to_ijk(coords, affine):
    """
    Converts `coords` in `affine` space to cartesian space

    Parameters
    ----------
    coords : (N, 3) array_like
        Image coordinate (xyz) values
    affine : (4, 4) array_like
        Affine matrix containing displacement + boundary

    Returns
    ------
    ijk : (N, 3) numpy.ndarray
        Provided `coords` in cartesian space
    """
    coords = _check_coord_inputs(coords)
    ijk_coords = np.linalg.solve(affine, coords)[:3].T.astype(int)
    return ijk_coords


def expand_roi(coords, dilation=1, return_array=True):
    """
    Expands coordinates ``coords`` to include neighboring coordinates

    Computes all possible coordinates of distance ``dilation`` from ``coords``.
    Returns a generator (``itertools.product``) by default, but can return an
    array if desired (if ``return_array=True``).

    Parameters
    ----------
    coords : (1, 3) array_like
        List of ijk values for coordinate in 3D space
    dilation : int, optional
        How many neighboring voxels to expand around `coords`. Default: 1
    return_array : bool, optional
        Whether to return array instead of generator. Default: True

    Returns
    -------
    coords : (27, 3) generator
        Coordinates of expanded ROI
    """

    def expand(x, d=1):
        return np.arange((x - d), (x + d + 1), dtype=int)

    # return all combinations of coordinates
    coords = np.squeeze(coords)
    gen = itertools.product(*[expand(n, d=dilation) for n in coords])

    # coerce to array if desired
    if return_array:
        return np.asarray(list(gen))

    return gen
