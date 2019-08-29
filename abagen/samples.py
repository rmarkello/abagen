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
     ('9352', 'sulci & spaces', 'other')),
    fields=('id', 'name', 'structure')
)


def _replace_mni_coords(annotation):
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

    annotation = io.read_annotation(annotation)
    mni_coords = coords.loc[annotation.well_id]
    annotation[['mni_x', 'mni_y', 'mni_z']] = np.asarray(mni_coords)

    return annotation


def _get_struct(structure_path):
    """
    Gets overarching "structure" of region defined by `structure_path`

    Structure here is defined as one of ['cortex', 'subcortex', 'cerebellum']

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
        One of ['cortex', 'subcortex', 'cerebellum'] or None
    """

    structure_path = set(structure_path.split('/'))
    ids = list(set(ONTOLOGY.value_set('id')).intersection(structure_path))

    try:
        return ONTOLOGY.structure[ids[0]]
    except IndexError:
        return


def drop_mismatch_samples(annotation, ontology, corrected=True):
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
    keep : pandas.DataFrame
        Annotation data
    """

    # read in data files
    annot = io.read_annotation(annotation)
    ont = io.read_ontology(ontology)

    # add hemisphere + brain "structure" designation to annotation data
    hemi = dict(zip(ont['id'], ont['hemisphere']))
    struct = dict(zip(ont['id'], ont['structure_id_path'].apply(_get_struct)))
    annot = annot.assign(hemisphere=annot['structure_id'].replace(hemi),
                         structure=annot['structure_id'].replace(struct))

    # correct MNI coordinates, if requested
    if corrected:
        annot = _replace_mni_coords(annot)

    # only keep samples with consistent hemisphere + MNI coordinate designation
    keep = annot.query('(hemisphere == "L" & mni_x < 0)'
                       '| (hemisphere == "R" & mni_x > 0)')

    return keep


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
        `atlas`. Must have _at least_ columns 'id', 'hemisphere', and
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
