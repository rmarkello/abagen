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

from . import datasets, io, probes, process, utils

import logging
logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)
lgr = logging.getLogger('abagen')
lgr_levels = dict(zip(range(3), [40, 20, 10]))


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
        Parcel labels for `S` samples, as returned by e.g., `label_samples()`
    labels : (L,) array_like, optional
        All possible labels for parcellation (to account for possibility that
        some parcels have NO expression data). Default: None
    metric : str or func, optional
        Mechanism by which to collapse across samples within a parcel. If a
        str, should be in ['mean', 'median']; if a function, should be able to
        accept an `N`-dimensional input and the `axis` keyword argument and
        return an `N-1`-dimensional output. Default: 'mean'

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
                                      left_index=True,
                                      right_index=True)
                               .groupby('label')
                               .aggregate(metric)
                               .append(labels)
                               .drop([0])
                               .sort_index()
                               .rename_axis('label'))

    return gene_by_label


def get_expression_data(atlas, atlas_info=None, *, exact=True,
                        tolerance=2, metric='mean', ibf_threshold=0.5,
                        probe_selection='diff_stability',
                        corrected_mni=True, reannotated=True,
                        return_counts=False, return_donors=False,
                        donors='all', data_dir=None, verbose=1):
    """
    Assigns microarray expression data to ROIs defined in `atlas`

    This function aims to provide a workflow for generating pre-processed,
    microarray expression data from the Allen Human Brain Atlas ([A2]_) for
    abitrary `atlas` designations. First, some basic filtering of genetic
    probes is performed, including:

        1. Intensity-based filtering of microarray probes to remove probes that
           do not exceed a certain level of background noise (specified via the
           `ibf_threshold` parameter), and
        2. Selection of a single, representative probe (or collapsing across
           probes) for each gene, specified via the `probe_selection`
           parameter.

    Tissue samples are then matched to parcels in the defined `atlas` for each
    donor. If `atlas_info` is provided then this matching is constrained by
    both hemisphere and tissue class designation (e.g., cortical samples from
    the left hemisphere are only matched to ROIs in the left cortex,
    subcortical samples from the right hemisphere are only matched to ROIs in
    the left subcortex); see the `atlas_info` parameter description for more
    information.

    Matching of microarray samples to parcels in `atlas` is done via a multi-
    step process:

        1. Determine if the sample falls directly within a parcel,
        2. Check to see if there are nearby parcels by slowly expanding the
           search space to include nearby voxels, up to a specified distance
           (specified via the `tolerance` parameter),
        3. If there are multiple nearby parcels, the sample is assigned to the
           closest parcel, as determined by the parcel centroid.

    If at any step a sample can be assigned to a parcel the matching process is
    terminated. If multiple sample are assigned to the same parcel they are
    aggregated with the metric specified via the `metric` parameter. More
    control over the sample matching can be obtained by setting the `exact`
    parameter; see the parameter description for more information.

    Once all samples have been matched to parcels for all supplied donors, the
    microarray expression data are normalized within-donor via a scaled robust
    sigmoid (SRS) procedure before being combined across donors via the
    supplied `metric`.

    Parameters
    ----------
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
        matched to that parcel. This is only considered if the sample is not
        directly within a parcel. Default: 2
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
    probe_selection : str, optional
        Selection method for subsetting (or collapsing across) probes that
        index the same gene. Must be one of 'average', 'max_intensity',
        'max_variance', 'pc_loading', 'corr_variance', 'corr_intensity', or
        'diff_stability'; see Notes for more information. Default:
        'diff_stability'
    corrected_mni : bool, optional
        Whether to use the "corrected" MNI coordinates shipped with the
        `alleninf` package instead of the coordinates provided with the AHBA
        data when matching tissue samples to anatomical regions. Default: True
    reannotated : bool, optional
        Whether to use reannotated probe information provided by [A1]_ instead
        of the default probe information from the AHBA dataset. Using
        reannotated information will discard probes that could not be reliably
        matched to genes. Default: True
    return_counts : bool, optional
        Whether to return how many samples were assigned to each parcel in
        `atlas` for each donor. Default: False
    return_donors : bool, optional
        Whether to return donor-level expression arrays instead of aggregating
        expression across donors with provided `metric`. Default: False
    donors : list, optional
        List of donors to use as sources of expression data. Can be either
        donor numbers or UID. If not specified will use all available donors.
        Default: 'all'
    data_dir : str, optional
        Directory where expression data should be downloaded (if it does not
        already exist) / loaded. If not specified will use the current
        directory. Default: None
    verbose : int, optional
        Specifies verbosity of status messages to display during workflow.
        Higher numbers increase verbosity of messages while zero suppresses all
        messages. Default: 1

    Returns
    -------
    expression : (R, G) :class:`pandas.DataFrame`
        Microarray expression for `R` regions in `atlas` for `G` genes,
        aggregated across donors, where the index corresponds to the unique
        integer IDs of `atlas` and the columns are gene names. If
        ``return_donors`` is set to ``True`` then this is a list of (R, G)
        dataframes, one for each donor.
    counts : (R, D) :class:`pandas.DataFrame`
        Number of samples assigned to each of `R` regions in `atlas` for each
        of `D` donors (if multiple donors were specified); only returned if
        ``return_counts`` is set to ``True``.

    Notes
    -----
    The following methods can be used for collapsing across probes when
    multiple probes are available for the same gene.

    1. ``probe_selection='average'``

    Takes the average of expression data across all probes indexing the same
    gene; ``probe_selection='mean'`` is identical.

    2. ``probe_selection='max_intensity'``

    Selects the probe with the maximum average expression across samples from
    all donors.

    3. ``probe_selection='max_variance'``

    Selects the probe with the maximum variance in expression across samples
    from all donors.

    4. ``probe_selection='pc_loading'``

    Selects the probe with the maximum loading along the first principal
    component of a decomposition performed across samples from all donors.

    5. ``probe_selection='corr_intensity'``

    Selects the probe with the maximum correlation to other probes from the
    same gene when >2 probes exist; otherwise, uses the same procedure as
    `max_intensity`.

    6. ``probe_selection='corr_variance'``

    Selects the probe with the maximum correlation to other probes from the
    same gene when >2 probes exist; otherwise, uses the same procedure as
    `max_varance`.

    7. ``probe_selection='diff_stability'``

    Selects the probe with the most consistent pattern of regional variation
    across donors (i.e., the highest average correlation across brain regions
    between all pairs of donors).

    References
    ----------
    .. [A1] Arnatkevic̆iūtė, A., Fulcher, B. D., & Fornito, A. (2019). A
       practical guide to linking brain-wide gene expression and neuroimaging
       data. NeuroImage, 189, 353-367.
    .. [A2] Hawrylycz, M.J. et al. (2012) An anatomically comprehensive atlas
       of the adult human transcriptome. Nature, 489, 391-399.
    """

    # set logging verbosity level
    lgr.setLevel(lgr_levels.get(int(verbose), 1))

    # load atlas_info, if provided
    atlas = check_niimg_3d(atlas)
    if atlas_info is not None:
        atlas_info = utils.check_atlas_info(atlas, atlas_info)

    # get combination functions
    metric = utils.check_metric(metric)

    # check probe_selection input
    if probe_selection not in probes.SELECTION_METHODS:
        raise ValueError('Provided probe_selection method is invalid, must be '
                         'one of {}. Received value: \'{}\''
                         .format(list(probes.SELECTION_METHODS),
                                 probe_selection))

    # fetch files (downloading if necessary)
    files = datasets.fetch_microarray(data_dir=data_dir, donors=donors,
                                      verbose=verbose)
    for key in ['microarray', 'probes', 'annotation', 'pacall', 'ontology']:
        if key not in files:
            raise KeyError('Provided `files` dictionary is missing {}. '
                           'Please check inputs.'.format(key))

    # get some info on the number of subjects, labels in `atlas_img`
    num_subj = len(files.microarray)
    all_labels = utils.get_unique_labels(atlas)
    if not exact:
        lgr.info('Inexact matching of samples requested; pre-calculating ROI '
                 'centroids for {} labels in provided atlas image'
                 .format(len(all_labels)))
        centroids = utils.get_centroids(atlas, labels=all_labels)

    # get dataframe of probe information (reannotated or otherwise)
    if reannotated:
        lgr.info('Reannotating microarray probes with information from '
                 'Arnatkevic̆iūtė et al., 2018, NeuroImage')
        probe_info = probes.reannotate_probes(files.probes[0])
    else:
        probe_info = io.read_probes(files.probes[0])

    # intensity-based filtering of probes
    probe_info = probes.filter_probes(files.pacall, probe_info,
                                      threshold=ibf_threshold)
    lgr.info('{} genes survive intensity-based filtering with threshold of {}'
             .format(len(np.unique(probe_info['gene_symbol'])), ibf_threshold))

    # get probe-reduced microarray expression data for all donors based on
    # selection method; this will be a list of gene x sample dataframes (one
    # for each donor)
    lgr.info('Reducing probes indexing same gene with provided method: {}'
             .format(probe_selection))
    microarray = probes.collapse_probes(files.microarray, files.annotation,
                                        probe_info, method=probe_selection)

    expression, missing = [], []
    counts = pd.DataFrame(np.zeros((len(all_labels) + 1, num_subj)),
                          index=np.append([0], all_labels))
    for subj in range(num_subj):
        # get rid of samples whose coordinates don't match ontological profile
        annotation = process.drop_mismatch_samples(files.annotation[subj],
                                                   files.ontology[subj],
                                                   corrected=corrected_mni)

        # subset representative probes + samples from microarray data
        samples = microarray[subj].loc[annotation.index]

        # assign samples to regions and aggregate samples w/i the same region
        sample_labels = label_samples(annotation, atlas,
                                      atlas_info=atlas_info,
                                      tolerance=tolerance)
        expression += [group_by_label(samples, sample_labels,
                                      all_labels, metric=metric)]
        lgr.info('{:>3} / {} samples matched to ROIs for donor #{}'
                 .format(np.sum(np.asarray(sample_labels) != 0),
                         len(annotation),
                         datasets.WELL_KNOWN_IDS.value_set('subj')[subj]))

        # get counts of samples collapsed into each ROI
        labs, num = np.unique(sample_labels, return_counts=True)
        counts.loc[labs, subj] = num

        # if we don't want to do exact matching then cache which parcels are
        # missing data and the expression data for the closest sample to that
        # parcel; we'll use this once we've iterated through all donors
        if not exact:
            coords = utils.xyz_to_ijk(annotation[['mni_x', 'mni_y', 'mni_z']],
                                      atlas.affine)
            empty = ~np.in1d(all_labels, labs)
            closest, dist = utils.closest_centroid(coords, centroids[empty],
                                                   return_dist=True)
            closest = samples.loc[annotation.iloc[closest].index]
            empty = all_labels[empty]
            closest.index = pd.Series(empty, name='label')
            missing += [(closest, dict(zip(empty, np.diag(dist))))]

    # check for missing ROIs and fill in, as needed
    if not exact:
        # find labels that are missing across all donors
        empty = reduce(set.intersection, [set(f.index) for f, d in missing])
        lgr.info('Matching {} ROIs with no microarray data to nearest samples'
                 .format(len(empty)))
        for roi in empty:
            # find donor with sample closest to centroid of empty parcel
            ind = np.argmin([d.get(roi) for f, d in missing])
            donor = datasets.WELL_KNOWN_IDS.value_set("subj")[ind]
            lgr.debug('Assigning sample from donor {} to ROI {}'
                      .format(donor, roi))
            # assign expression data from that sample and add to count
            expression[ind].loc[roi] = missing[ind][0].loc[roi]
            counts.loc[roi, ind] += 1

    # normalize data with SRS
    expression = [process.normalize_expression(e) for e in expression]

    # aggregate across donors if individual donor dataframes not requested
    if not return_donors:
        expression = process.aggregate_donors(expression, metric)

    # drop the "zero" label from the counts dataframe (this is background)
    if return_counts:
        return expression, counts.iloc[1:]

    return expression
