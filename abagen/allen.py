# -*- coding: utf-8 -*-
"""
Functions for mapping AHBA microarray dataset to atlases and and parcellations
in MNI space
"""
from functools import reduce
from nilearn._utils import check_niimg_3d
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from sklearn.utils import Bunch
from abagen import io, process, utils


def _assign_sample(sample, atlas, sample_info=None, atlas_info=None,
                   tolerance=3):
    """
    Determines what ROI `sample` belongs to in `atlas`

    Parameters
    ----------
    sample : (1, 3) array_like
        Coordinates (ijk) of microarray probe in `atlas` space
    atlas : niimg-like object
        ROI image, where each ROI should be identified with a unique
        integer ID
    tolerance : int, optional
        Distance that sample must be from ROI for it to be considered within
        that ROI. This is only used if the sample is not directly inside a ROI.
        Default: 3mm

    Returns
    -------
    label : int
        ROI label of `sample`
    """

    # pull relevant info from atlas
    label_data = check_niimg_3d(atlas).get_data()

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
            nz_labels[nz_labels == old_label] = new_label
        labels, counts = np.unique(nz_labels, return_counts=True)

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
        Tenative label for sample in `sample_info`
    sample_info : pandas.DataFrame
    atlas_info : pandas.DataFrame

    Returns
    -------
    label : int
        New tentative label
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
    Matches all samples in `annotation` to closest ROIs in `atlas_info`

    Attempts to place each sample provided in `annotation` into a labelled ROI
    in `atlas_info`, where the latter is a 3D niimg-like object that contains
    ROI data, each with a unique integer ID.

    The function tries to best match the microarray sample to the ROIs by:
    1. Determining if the sample falls directly within a ROI.
    2. Checking to see if there are nearby ROIs by slowly expanding radius
       around sample (up to radius = `tolerance`); if there are multiple nearby
       ROIs, determining which ROI is closer (i.e., by centroid of the ROI).
    3. If `user_centroids` is set, find the ROI with the closest centroid

    If there is still no ROI for a given sample, it is assigned as "NA" and
    should likely be ignored in further analyses.

    Parameters
    ----------
    annotation : (S, 13) pandas.DataFrame
        Annotation information, where `S` is samples
    atlas : niimg-like object
        ROI image, where each ROI should be identified with a unique integer ID
    atlas_info : pandas.DataFrame, optional
        Filepath to or pre-loaded dataframe containing information about
        `atlas`. Must have _at least_ columns 'id', 'hemisphere', and
        'structure' containing information mapping atlas IDs to hemisphere and
        broad structural class (i.e., "cortex", "subcortex", "cerebellum").
        Default: None
    tolerance : int, optional
        Distance (in mm) that a sample must be from an ROI boundary for it to
        be considered within that ROI. This is only used if the sample is not
        directly within a ROI. Default: 2

    Returns
    -------
    labels : (S, 1) pandas.DataFrame
        Dataframe with ROI labels
    """

    # get annotation and atlas data
    annotation = io.read_annotation(annotation)
    atlas = check_niimg_3d(atlas)
    label_data, affine = atlas.get_data(), atlas.affine

    # load atlas_info, if provided
    if atlas_info is not None:
        atlas_info = utils.check_atlas_info(atlas, atlas_info)

    # get ijk coordinates for microarray samples and find labels
    g_ijk = utils.xyz_to_ijk(annotation[['mni_x', 'mni_y', 'mni_z']], affine)
    labelled_samples = label_data[g_ijk[:, 0], g_ijk[:, 1], g_ijk[:, 2]]

    # if sample coordinates aren't directly inside a parcel, increment radius
    # around sample up to `tolerance` to try and find nearby parcels.
    # if still no parcel, then ignore this sample
    for idx in np.where(labelled_samples == 0)[0]:
        label, tol = labelled_samples[idx], 1
        while label == 0 and tol <= tolerance:
            label = _assign_sample(g_ijk[[idx]], atlas,
                                   sample_info=annotation.iloc[idx],
                                   atlas_info=atlas_info,
                                   tolerance=tol)
            tol += 1
        labelled_samples[idx] = label

    return pd.DataFrame(labelled_samples, dtype=int,
                        columns=['label'], index=annotation.index)


def group_by_label(microarray, sample_labels, labels=None, metric='mean'):
    """
    Averages expression data in `microarray` over samples with same label

    Parameters
    ----------
    microarray : (S, G) pandas.DataFrame
        Microarray expression data, where `S` is samples and `G` is genes
    sample_labels : (S, 1) pandas.DataFrame
        ROI labels for `S` samples, as returned by e.g., `label_samples()`
    labels : (L,) array_like, optional
        All possible labels for parcellation scheme (to account for possibility
        that some parcels have NO expression data). Default: None
    metric : str or func
        Mechanism by which to collapse across samples within an ROI. If str,
        should be in ['mean', 'median'].

    Returns
    -------
    gene_by_label : (L, G) pandas.DataFrame
        Microarray expression data
    """

    # get combination function
    metric = utils.check_metric(metric)

    # get missing labels
    if labels is not None:
        missing = np.setdiff1d(labels, sample_labels)
        labels = pd.DataFrame(columns=microarray.columns,
                              index=pd.Series(missing, name='label'))

    gene_by_label = (microarray.merge(sample_labels,
                                      left_index=True, right_index=True)
                               .groupby('label')
                               .aggregate(metric)
                               .append(labels)
                               .drop([0])
                               .sort_index()
                               .rename_axis('label'))

    return gene_by_label


def get_expression_data(files, atlas, atlas_info=None, *,
                        exact=True, tolerance=2, metric='mean',
                        ibf_threshold=0.5, corrected_mni=True,
                        return_counts=False):
    """
    Assigns microarray expression data in `files` to ROIs defined in `atlas`

    Parameters
    ----------
    files : dict
        Optimally obtained by calling `abagen.fetch_microarray()`, this should
        be a dict with keys ['annotation', 'microarray', 'ontology', 'pacall',
        'probes'], where corresponding values are lists of filepaths to
        downloaded AHBA data
    atlas : niimg-like object
        A parcellation image in MNI space, where each parcel is identified by a
        unique integer ID
    atlas_info : str or :class:`pandas.DataFrame`, optional
        Filepath to or pre-loaded dataframe containing information about
        `atlas`. Must have at least columns 'id', 'hemisphere', and 'structure'
        containing information mapping atlas IDs to hemisphere (i.e, "L", "R")
        and broad structural class (i.e., "cortex", "subcortex", "cerebellum").
        Default: None
    exact : bool, optional
        Whether to use exact matching of donor tissue samples to parcels in
        `atlas`. If True, this function will match tissue samples to parcels
        within `threshold` mm of the sample; any samples that are beyond
        `threshold` mm of a parcel will be discarded. This may result in some
        parcels having no assigned sample / expression data. If False, the
        default matching procedure will be performed and followed by a check
        for parcels with no assigned samples; any such parcels will be matched
        to the nearest sample (nearest defined as the sample with the closest
        Euclidean distance to the parcel centroid). Default: True
    tolerance : int, optional
        Distance (in mm) that a sample must be from a parcel for it to be
        matched to that parcel. Default: 3
    metric : str or func, optional
        Mechanism by which to collapse across donors, if input `files` provides
        multiple donor datasets. If a str, should be in ['mean', 'median']; if
        a function, should be able to accept an `N`-dimensional input and the
        `axis` keyword argument and return an `N-1`-dimensional output.
        Default: 'mean'
    ibf_threshold : [0, 1] float, optional
        Threshold for intensity-based filtering specifying. This number should
        specify the ratio of samples, across all supplied donors, for which a
        probe must have signal above background noise in order to be retained.
        Default: 0.5
    corrected_mni : bool, optional
        Whether to use the "corrected" MNI coordinates shipped with the
        `alleninf` package instead of the coordinates provided with the AHBA
        data when matching tissue samples to anatomical regions. Default: True
    return_counts : bool, optional
        Whether to return how many samples were assigned to each parcel in
        `atlas` for each donor. Default: False

    Returns
    -------
    expression : (R x G) :class:`pandas.DataFrame`
        Microarray expression for `R` regions in `atlas` across `G` genes,
        averaged across donors.
    counts : :class:`pandas.DataFrame`
        Number of samples assigned to each parcel in `atlas` (by donor, if
        multiple donors provided); only provided if `return_counts=True`
    """

    # coerce to Bunch in case a simple dictionary was provided
    files = Bunch(**files)
    for key in ['microarray', 'probes', 'annotation', 'pacall', 'ontology']:
        if key not in files:
            raise KeyError('Provided `files` dictionary is missing {}. '
                           'Please check inputs.'.format(key))

    # load atlas_info, if provided
    atlas = check_niimg_3d(atlas)
    if atlas_info is not None:
        atlas_info = utils.check_atlas_info(atlas, atlas_info)

    # get combination functions
    metric = utils.check_metric(metric)

    # get some info on the number of subjects, labels in `atlas_img`
    num_subj = len(files.microarray)
    all_labels = utils.get_unique_labels(atlas)
    if not exact:
        centroids = utils.get_centroids(atlas, labels_of_interest=all_labels)

    probes = process.filter_probes(files.pacall, files.probes,
                                   threshold=ibf_threshold)
    probes = process.get_stable_probes(files.microarray, files.annotation,
                                       probes)

    expression, counts = [], np.zeros((len(all_labels) + 1, num_subj))
    missing = []
    for subj in range(num_subj):
        # get rid of samples whose coordinates don't match ontological profile
        annotation = process.drop_mismatch_samples(files.annotation[subj],
                                                   files.ontology[subj],
                                                   corrected=corrected_mni)

        # get representative probes + samples from microarray data
        microarray = io.read_microarray(files.microarray[subj])
        samples = microarray.loc[probes.index, annotation.index].T
        samples.columns = probes.gene_symbol

        # assign samples to regions and aggregate samples w/i the same region
        sample_labels = label_samples(annotation, atlas,
                                      atlas_info=atlas_info,
                                      tolerance=tolerance)
        expression += [group_by_label(samples, sample_labels,
                                      all_labels, metric=metric)]

        # get counts of samples collapsed into each ROI
        labs, num = np.unique(sample_labels, return_counts=True)
        counts[labs, subj] = num

        if not exact:
            coords = utils.xyz_to_ijk(annotation[['mni_x', 'mni_y', 'mni_z']],
                                      atlas.affine)
            empty = np.setdiff1d(all_labels, labs)
            closest, dist = utils.closest_centroid(coords, centroids[empty],
                                                   return_dist=True)
            closest = samples.loc[annotation.iloc[closest].index]
            closest.index = pd.Series(empty, name='label')
            missing += [(closest, dict(zip(empty, np.diag(dist))))]

    # check for missing ROIs and fill in, as needed
    if not exact:
        empty = reduce(set.intersection, [set(f.index) for f, d in missing])
        for roi in empty:
            ind = np.argmin([d.get(roi) for f, d in missing])
            expression[ind].loc[roi] = missing[ind][0].loc[roi]
            counts[roi, ind] += 1

    # normalize data and aggregate across donors
    expression = [process.normalize_expression(e) for e in expression]
    expression = pd.concat(expression).groupby('label').aggregate(metric)

    if return_counts:
        return expression, counts[1:]

    return expression
