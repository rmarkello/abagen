# -*- coding: utf-8 -*-
"""
Functions for cleaning and processing the AHBA microarray dataset
"""

from pkg_resources import resource_filename

from nibabel.volumeutils import Recoder
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from . import io, utils

# AHBA structure IDs corresponding to different brain parts
ONTOLOGY = Recoder(
    (('4008', 'cerebral cortex', 'cortex'),
     ('4275', 'cerebral nuclei', 'subcortex'),
     ('4391', 'diencephalon', 'subcortex'),
     ('9001', 'mesencephalon', 'subcortex'),
     ('4696', 'cerebellum', 'cerebellum'),
     ('9131', 'pons', 'brainstem'),
     ('9512', 'myelencephalon', 'brainstem'),
     ('9218', 'white matter', 'white matter'),
     ('9352', 'sulci & spaces', 'other'),
     ('4219', 'hippocampal formation', 'subcortex')),
    fields=('id', 'name', 'structure')
)


def update_mni_coords(annotation):
    """
    Replaces MNI coords in `annotation` with corrected coords from `alleninf`

    Parameters
    ----------
    annotation : str
        Annotation file from Allen Brain Institute. Optimally obtained by
        calling `abagen.fetch_microarray()` and accessing the `annotation`
        attribute on the resulting object

    Returns
    -------
    corrected : pandas.DataFrame
        Annotation data with corrected MNI coordinates

    References
    ----------
    Updated MNI coordinates taken from https://github.com/chrisfilo/alleninf,
    which is licensed under the BSD-3 (reproduced here):

    Copyright (c) 2018, Krzysztof Gorgolewski
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.
    * Neither the name of the copyright holder nor the names of its
      contributors may be used to endorse or promote products derived from
      this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
    LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.
    """

    coords = resource_filename('abagen', 'data/corrected_mni_coordinates.csv')
    coords = pd.read_csv(coords).rename(dict(corrected_mni_x='mni_x',
                                             corrected_mni_y='mni_y',
                                             corrected_mni_z='mni_z'),
                                        axis=1)
    coords = coords.set_index('well_id')

    annotation = io.read_annotation(annotation, copy=True)

    # basic check that all well_ids in annotation are present in coords
    # a future pandas update may cause this to raise a KeyError but we want
    # this to raise a KeyError NOW
    diff = np.setdiff1d(annotation['well_id'], coords.index)
    if len(diff) > 0:
        raise KeyError('Provided annotation file has well IDs that do not '
                       'exist in updated MNI coordinate file from `alleninf`. '
                       'Please check input annotation file and try again. '
                       'Unknown well IDs: {}'.format(diff))

    mni_coords = coords.loc[annotation.well_id]
    annotation[['mni_x', 'mni_y', 'mni_z']] = np.asarray(mni_coords)

    return annotation


def _get_struct(structure_path):
    """
    Gets overarching "structure" of region defined by `structure_path`

    Structure here is defined as being one of {'cortex', 'subcortex',
    'cerebellum', 'brainstem', 'white matter', 'other'}.

    Parameters
    ----------
    structure_path : str
        As obtained from an ontology file from Allen Brain Institute. Optimally
        obtained by calling `abagen.fetch_microarray()`, loading the
        `ontology` attribute on the resulting object, and querying the
        `structure_id_path` column from the resulting dataframe

    Returns
    -------
    structure : str
        Structure, or None if unable to identify a corresponding structure
    """

    structures = set(structure_path.strip('/').split('/'))
    ids = sorted(set(ONTOLOGY.value_set('id')) & structures,
                 key=lambda x: structure_path.index(x))

    try:
        return ONTOLOGY.structure[ids[-1]]
    except IndexError:
        return


def drop_mismatch_samples(annotation, ontology):
    """
    Removes samples from `annotation` whose coordinates do not match `ontology`

    Checks MNI coordinates specified in `annotation` and matches them to L/R
    hemisphere designation in `ontology`; samples who do not match (e.g., L
    hemisphere designation in `ontology` but X coordinate < 0) are removed.

    Optionally updates MNI coordinates in `annotation` (see `corrected`).

    Parameters
    ----------
    annotation : str
        Annotation file from Allen Brain Institute. Optimally obtained by
        calling `abagen.fetch_microarray()` and accessing the `annotation`
        attribute on the resulting object
    ontology : str
        Ontology file from Allen Brain Institute. Optimally obtained by
        calling `abagen.fetch_microarray()` and accessing the `ontology`
        attribute on the resulting object
    corrected : bool, optional
        Whether to use the "corrected" MNI coordinates shipped with the
        `alleninf` package instead of the coordinates provided with the AHBA
        data. Default: True

    Returns
    -------
    annot : pandas.DataFrame
        Annotation data
    """

    # read in data files
    annot = io.read_annotation(annotation)
    ont = io.read_ontology(ontology).set_index('id')

    # get hemisphere and structure path
    hemisphere = np.asarray(ont.loc[annot['structure_id'], 'hemisphere'])
    structure = np.asarray(ont.loc[annot['structure_id'], 'structure_id_path']
                              .apply(_get_struct))

    # add hemisphere + brain "structure" designation to annotation data and
    # only keep samples with consistent hemisphere + MNI coordinate designation
    annot = annot.assign(hemisphere=hemisphere, structure=structure) \
                 .query('(hemisphere == "L" & mni_x < 0) '
                        '| (hemisphere == "R" & mni_x > 0) '
                        '| (hemisphere.isna() & mni_x == 0)')

    return annot


def _assign_sample(sample, atlas, sample_info=None, atlas_info=None,
                   tolerance=2):
    """
    Determines which parcel `sample` belongs to in `atlas`

    Parameters
    ----------
    sample : (1, 3) array_like
        Coordinates (ijk) of microarray sample in `atlas` space
    atlas : niimg-like object
        ROI image, where each ROI should be identified with a unique
        integer ID
    sample_info : pandas.DataFrame
        A single row of an `annotation` file, corresponding to the given sample
    atlas_info : pandas.DataFrame,
        Dataframe containing information about the specified `atlas`. Must have
        _at least_ columns 'id', 'hemisphere', and 'structure' containing
        information mapping atlas IDs to hemisphere and broad structural class
        (i.e., "cortex", "subcortex", "cerebellum"). Default: None
    tolerance : int, optional
        Distance (in mm) that a sample must be from a parcel for it to be
        matched to that parcel. This is only considered if the sample is not
        directly within a parcel. Default: 2

    Returns
    -------
    label : int
        Parcel label of `sample`
    """

    # pull relevant info from atlas
    label_data = utils.check_img(atlas).get_data()

    # expand provided coordinates to include those w/i `tolerance` of `coords`
    # set a hard euclidean distance limit to account for different voxel sizes
    coords = utils.expand_roi(sample, dilation=tolerance, return_array=True)
    coords = coords[cdist(sample, coords).squeeze() < tolerance]

    # grab non-zero labels for expanded coordinates
    possible_labels = label_data[coords[:, 0], coords[:, 1], coords[:, 2]]
    nz_labels = possible_labels[possible_labels.nonzero()]
    labels, counts = np.unique(nz_labels, return_counts=True)

    # if atlas_info and sample_info are provided, drop potential labels who
    # don't match hemisphere or structural class defined in `sample_info`
    if atlas_info is not None and sample_info is not None:
        for old_label in labels:
            new_label = _check_label(old_label, sample_info, atlas_info)
            if old_label != new_label:
                nz_labels[nz_labels == old_label] = new_label
        labels, counts = np.unique(nz_labels[nz_labels.nonzero()],
                                   return_counts=True)

    # if there is still nothing in the vicinity, return 0
    if labels.size == 0:
        return 0
    # if there is only one ROI in the vicinity, use that
    elif labels.size == 1:
        return labels[0]

    # if more than one ROI in the vicinity, return the most frequent
    indmax, = np.where(counts == counts.max())
    if indmax.size == 1:
        return labels[indmax[0]]

    # if two or more parcels tied for neighboring frequency, use ROI
    # with closest centroid to `coords`
    centroids = utils.get_centroids(atlas, labels)
    return labels[utils.closest_centroid(sample, centroids)]


def _check_label(label, sample_info, atlas_info):
    """
    Checks that `label` defined by `sample_info` is coherent with `atlas_info`

    Parameters
    ----------
    label : int
        Tenative label for sample described by `sample_info`
    sample_info : pandas.DataFrame
        A single row of an `annotation` file, corresponding to the given sample
    atlas_info : pandas.DataFrame,
        Dataframe containing information about the atlas of interest. Must have
        _at least_ columns 'id', 'hemisphere', and 'structure' containing
        information mapping atlas IDs to hemisphere and broad structural class
        (i.e., "cortex", "subcortex", "cerebellum"). Default: None

    Returns
    -------
    label : int
        New label for sample
    """

    cols = ['hemisphere', 'structure']

    if label != 0:
        sample_info = sample_info[cols]
        atlas_info = atlas_info.loc[label][cols]
        if not np.all(sample_info.values == atlas_info.values):
            label = 0

    return label


def label_samples(annotation, atlas, atlas_info=None, tolerance=2):
    """
    Matches all microarray samples in `annotation` to parcels in `atlas`

    Attempts to place each sample provided in `annotation` into a parcel in
    `atlas`, where the latter is a 3D niimg-like object that contains parcels
    each idnetified by a unique integer ID.

    The function tries to best match samples in `annotation` to parcels defined
    in `atlas` by:

        1. Determining if the sample falls directly within a parcel,
        2. Checking to see if there are nearby parcels by slowly expanding the
           search space to include nearby voxels, up to a specified distance
           (specified via the `tolerance` parameter),
        3. Assigning the sample to the closest parcel if there are multiple
           nearby parcels, where closest is determined by the parcel centroid.

    If at any step a sample can be assigned to a parcel the matching process is
    terminated. If there is still no parcel for a given sample after this
    process the sample is provided a label of 0.

    Parameters
    ----------
    annotation : (S, 13) pandas.DataFrame
        Pre-loaded annotation information for a given AHBA donor
    atlas : niimg-like object
        A parcellation image in MNI space, where each parcel is identified by a
        unique integer ID
    atlas_info : pandas.DataFrame, optional
        Filepath to or pre-loaded dataframe containing information about
        `atlas`. Must have *at least* columns 'id', 'hemisphere', and
        'structure' containing information mapping atlas IDs to hemisphere and
        broad structural class (i.e., "cortex", "subcortex", "cerebellum").
        Default: None
    tolerance : int, optional
        Distance (in mm) that a sample must be from a parcel for it to be
        matched to that parcel. This is only considered if the sample is not
        directly within a parcel. Default: 2

    Returns
    -------
    labels : (S, 1) pandas.DataFrame
        Dataframe with parcel labels for each of `S` samples
    """

    # get annotation and atlas data
    annotation = io.read_annotation(annotation)
    atlas = utils.check_img(atlas)
    label_data, affine = atlas.get_data(), atlas.affine

    # load atlas_info, if provided
    if atlas_info is not None:
        atlas_info = utils.check_atlas_info(atlas, atlas_info)

    # get ijk coordinates for microarray samples and find labels
    g_ijk = utils.xyz_to_ijk(annotation[['mni_x', 'mni_y', 'mni_z']], affine)
    labelled = label_data[g_ijk[:, 0], g_ijk[:, 1], g_ijk[:, 2]]

    # if sample coordinates aren't directly inside a parcel, increment radius
    # around sample up to `tolerance` to try and find nearby parcels.
    # if still no parcel, then ignore this sample
    for idx in np.where(labelled == 0)[0]:
        label, tol = labelled[idx], 1
        while label == 0 and tol <= tolerance:
            label = _assign_sample(g_ijk[[idx]], atlas,
                                   sample_info=annotation.iloc[idx],
                                   atlas_info=atlas_info,
                                   tolerance=tol)
            tol += 1
        labelled[idx] = label

    return pd.DataFrame(labelled, dtype=int,
                        columns=['label'], index=annotation.index)


def mirror_samples(microarray, pacall, annotation, ontology, inplace=False):
    """
    Mirrors all tissue samples across hemispheres (L->R and R->L)

    Methodology follows [SA1]_ in that samples are mirrored bi-directionally
    rather than uni-directionally (R->L) as in [SA2]_.

    Parameters
    ----------
    microarray : list of str or pandas.DataFrame
        List of filepaths to microarray expression files from Allen Brain
        Institute (i.e., as obtained by calling :func:`abagen.fetch_microarray`
        and accessing the `microarray` attribute on the resulting object).
    pacall : list of str or pandas.DataFrame
        List of filepaths to pacall files from Allen Brain Institute (i.e., as
        obtained by calling :func:`abagen.fetch_microarray` and accessing the
        `pacall` attribute on the resulting object).
    annotation : list of str or pandas.DataFrame
        List of filepaths to annotation files from Allen Brain Institute (i.e.,
        as obtained by calling :func:`abagen.fetch_microarray` and accessing
        the `annotation` attribute on the resulting object).
    ontology : list of str or pandas.DataFrame
        List of filepaths to ontology files from Allen Brain Institute (i.e.,
        as obtained by calling :func:`abagen.fetch_microarray` and accessing
        the `ontology` attribute on the resulting object).
    inplace : bool, optional
        Whether to conserve memory by editing dataframes in-place instead of
        returning edited copies. Default: False

    Returns
    -------
    microarray,pacall,annotation : list of pandas.DataFrame
        Loaded input data with all samples duplicated across hemispheres

    References
    ----------
    .. [SA1] Gryglewski, G., Seiger, R., James, G. M., Godbersen, G. M.,
       Komorowski, A., Unterholzner, J., ... & Kasper, S. (2018). Spatial
       analysis and high resolution mapping of the human whole-brain
       transcriptome for integrative analysis in neuroimaging. NeuroImage, 176,
       259-267.

    .. [SA2] Romero-Garcia, R., Whitaker, K. J., Váša, F., Seidlitz, J., Shinn,
       M., Fonagy, P., ... & Vértes, P. E. (2018). Structural covariance
       networks are coupled to expression of genes enriched in supragranular
       layers of the human cortex. NeuroImage, 171, 256-267.
    """

    # FIXME: seems like there should be a better way to do this, but both str
    # and dataframes _are_ iterable just...not the way we want them to be here
    inp = [microarray, pacall, annotation, ontology]
    for n, i in enumerate(inp):
        if isinstance(i, (str, pd.DataFrame)):
            inp[n] = [i]

    # flipped will be a list of len-3 tuples: (microarray, pacall, annotation)
    flipped = [_mirror_samples(*i, inplace=inplace) for i in zip(*inp)]

    # unpack so that even if user assigns output to 1 variable they get a tuple
    microarray, pacall, annotation = [list(f) for f in zip(*flipped)]

    return microarray, pacall, annotation


def _mirror_samples(microarray, pacall, annotation, ontology, inplace=False):
    """
    Mirrors tissue samples across hemispheres for single donor

    microarray,pacall,annotation,ontology : str or pandas.DataFrame
        Filepath to {microarray,pacall,annotation,ontology} file
    inplace : bool, optional
        Whether to conserve memory by editing dataframes in-place instead of
        returning edited copies. Default: False

    Returns
    -------
    microarray,pacall,annotation : pandas.DataFrame
        Loaded input data with all samples duplicated across hemispheres
    """

    microarray = io.read_microarray(microarray, copy=not inplace)
    pacall = io.read_pacall(pacall, copy=not inplace)
    annotation = io.read_annotation(annotation, copy=not inplace)
    ontology = io.read_ontology(ontology)

    # take all lh and rh samples and flip x-coordinate
    # also update ontology information (structure_id/acronym/name) as this is
    # used when dropping mismatched samples later in the workflow
    lh = _mirror_ontology(annotation[annotation['mni_x'] < 0], ontology)
    rh = _mirror_ontology(annotation[annotation['mni_x'] > 0], ontology)
    for df in [lh, rh]:
        df['mni_x'] *= -1

    # grow microarray and pacall based on duplicated samples
    annotation = pd.concat([annotation, lh, rh])
    microarray = microarray.loc[:, annotation.index]
    pacall = pacall.loc[:, annotation.index]

    # now reset index / columns of all dataframes to be consecutive
    sids = pd.Series(range(len(annotation)), name='sample_id')
    microarray.columns = sids
    pacall.columns = sids
    annotation.index = sids

    return microarray, pacall, annotation


def _mirror_ontology(annotation, ontology):
    """
    Assumes all hemisphere assignments of structures in `annotation` are wrong

    Parameters
    ----------
    annotation,ontology : str
        Filepath to {annotation,ontology} file

    Returns
    -------
    annotation : pandas.DataFrame
        Loaded annotation with updated structure information
    """

    HEMI_SWAP = dict(L='R', R='L')

    annotation = io.read_annotation(annotation, copy=True)
    ontology = io.read_ontology(ontology)

    # structure IDs are specific to structure + hemisphere, so we can use this
    # to grab the hemisphere designation for each sample and flip that L<->R
    hemi = ontology.set_index('id') \
                   .loc[annotation['structure_id'], 'hemisphere'] \
                   .replace(HEMI_SWAP)

    # structure acronyms are distinct to structure but not to hemisphere, so we
    # can use this to grab all variations of a given structure and, with the
    # flipped hemisphere designations we just generated, find the ID/name of
    # the relevant structure in the new hemisphere
    acr = ontology.set_index(['acronym', 'hemisphere']) \
                  .loc[zip(annotation['structure_acronym'], hemi)] \
                  .reset_index()

    # update the original annotation with the new values
    annotation['structure_acronym'] = acr['acronym'].values
    annotation['structure_id'] = acr['id'].values
    annotation['structure_name'] = acr['name'].values

    return annotation
