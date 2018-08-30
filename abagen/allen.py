# -*- coding: utf-8 -*-
"""
Functions for mapping AHBA microarray dataset to atlases and and parcellations
in MNI space
"""
from nilearn._utils import check_niimg_3d
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from sklearn.utils import Bunch
from abagen import io, processing, utils


def assign_sample(sample, label_image, tolerance=3):
    """
    Determines what ROI `sample` belongs to in `label_image`

    Parameters
    ----------
    sample : (1, 3) array_like
        Coordinates (ijk) of microarray probe in `label_image` space
    label_image : niimg-like object
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

    # pull relevant info from label_image
    label_data = check_niimg_3d(label_image).get_data()

    # expand provided coordinates to include those w/i `tolerance` of `coords`
    # set a hard euclidean distance limit to account for different voxel sizes
    coords = utils.expand_roi(sample, dilation=tolerance, return_array=True)
    coords = coords[cdist(sample, coords).squeeze() < tolerance]

    # grab non-zero labels for expanded coordinates
    possible_labels = label_data[coords[:, 0], coords[:, 1], coords[:, 2]]
    nz_labels = possible_labels[possible_labels.nonzero()]
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
    centroids = utils.get_centroids(label_image, labels)
    return labels[utils.closest_centroid(sample, centroids)]


def _check_label(label, sample, atlas_info):
    """
    Checks that `label` assigned to `sample` is coherent with `atlas_info`

    Parameters
    ----------
    label : int
        Tenative label for `sample`
    sample
    atlas_info

    Returns
    -------
    label : int
    """

    if label != 0:
        cols = ['hemisphere', 'structure']
        sample = sample[cols].values
        atlas_check = atlas_info[atlas_info.id == label][cols].values
        if not np.all(sample == atlas_check):
            label = 0

    return label


def label_samples(annotation, atlas, atlas_info, tolerance=3,
                  use_centroids=False):
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
    atlas_info : niimg-like object
        ROI image, where each ROI should be identified with a unique integer ID
    tolerance : int, optional
        Distance that sample must be from ROI for it to be considered within a
        ROI. This is only used if the sample is not directly inside a ROI.
        Default: 3mm
    use_centroids : bool, optional
        If no ROI is within `tolerance` of a sample, assign sample to ROI with
        closest centroid. Default: False

    Returns
    -------
    labels : (S, 1) pandas.DataFrame
        Dataframe with ROI labels
    """

    # read annotation file, if provided
    if isinstance(annotation, str):
        annotation = io.read_sampleannot(annotation)

    # get image data
    atlas = check_niimg_3d(atlas)
    label_data, affine_trans = atlas.get_data(), atlas.affine

    # if we're going to use centroids, calculate them ahead of time
    if use_centroids:
        all_labels = utils.get_unique_labels(atlas)
        centroids = utils.get_centroids(atlas, all_labels)

    # grab xyz coordinates for microarray samples and convert to ijk
    g_xyz = annotation[['mni_x', 'mni_y', 'mni_z']].get_values()
    g_ijk = np.floor(utils.xyz_to_ijk(g_xyz, affine_trans)).astype(int)

    # get labels for all ijk values
    labelled_samples = label_data[g_ijk[:, 0], g_ijk[:, 1], g_ijk[:, 2]]

    # if coordinates aren't within the parcel, check for neighboring parcels
    # and slowly increase the radius around parcel up to `tolerance` to try
    # and find nearby parcels. if still no nearby parcel then ignore probe
    for idx in np.where(labelled_samples == 0)[0]:
        label, tol = labelled_samples[idx], 1
        sample = annotation.iloc[idx]
        while label == 0 and tol <= tolerance:
            label = assign_sample(g_ijk[[idx]], atlas, atlas_info, sample,
                                  tolerance=tol)
            tol += 1
        if label == 0 and use_centroids:
            label = all_labels[utils.closest_centroid(g_ijk[[idx]], centroids)]
            label = _check_label(label, annotation.iloc[idx], atlas_info)
        labelled_samples[idx] = label

    # return DataFrame for ease of use
    return pd.DataFrame(labelled_samples, dtype=int,
                        columns=['label'],
                        index=annotation.index)


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
    genes : (L, G) pandas.DataFrame
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


def label_rois(annotation, label_image, tolerance=3):
    """
    Matches all ROIs in `label_image` to closest samples in `annotation`

    Attempts to place each sample provided in `annotation` into a labelled ROI
    in `label_image`, where the latter is a 3D niimg-like object that contains
    ROI data, each with a unique integer ID.

    The function tries to best match the microarray sample to the ROIs by:
    1. Determining which samples falls directly within a ROI.
    2. If there is no sample assigned to an ROI, find the closest sample to the
       centroid of the ROI

    Parameters
    ----------
    annotation : (S, 13) pandas.DataFrame
        Annotation information, where `S` is samples
    label_image : niimg-like object
        ROI image, where each ROI should be identified with a unique integer ID
    tolerance : int, optional
        Not uses in this implementation.

    Returns
    -------
    labels : (L, ...) pandas.DataFrame
        Dataframe with lists of sample labels for each ROI
    """

    # read annotation file, if provided
    annotation = io.read_sampleannot(annotation)

    # get image data
    label_image = check_niimg_3d(label_image)
    label_data, affine_trans = label_image.get_data(), label_image.affine

    # calculate the centroids
    all_labels = utils.get_unique_labels(label_image)
    centroids = utils.get_centroids(label_image, all_labels)

    # grab xyz coordinates for microarray samples and convert to ijk
    g_xyz = annotation[['mni_x', 'mni_y', 'mni_z']].get_values()
    g_ijk = np.floor(utils.xyz_to_ijk(g_xyz, affine_trans)).astype(int)

    # get labels for all ijk values
    labelled_samples = label_data[g_ijk[:, 0], g_ijk[:, 1], g_ijk[:, 2]]

    # make a list of samples within each ROI
    # if no sample is assigned to an ROI, find closest sample to ROI centroid
    labelled_rois = []
    for nn in range(len(all_labels)):
        indices = [i for i, x in enumerate(labelled_samples) if x == nn + 1]
        labelled_rois.append(indices)
        if not labelled_rois[nn]:
            emptyroi = np.reshape(centroids[:, nn], (-1, 1))
            label = utils.closest_centroid(emptyroi, g_ijk.T)
            labelled_rois[nn] = [label]

    # return DataFrame for ease of use
    return pd.DataFrame(labelled_rois,
                        index=pd.Series(all_labels, name='label'))


def group_by_roi(microarray, roi_labels, labels=None, metric='mean'):
    """
    Averages expression data in `microarray` over samples within same ROI

    Parameters
    ----------
    microarray : (S, G) pandas.DataFrame
        Microarray expression data, where `S` is samples and `G` is genes
    roi_labels : (L, ...) pandas.DataFrame
        Lists of sample labels with varying lengths for `L` ROIs, as returned
        by e.g., `label_rois()`
    labels : (L,) array_like, optional
        All possible labels for parcellation scheme (to account for possibility
        that some parcels have NO expression data). Default: None
    metric : str or func
        Mechanism by which to collapse across samples within an ROI. If str,
        should be in ['mean', 'median'].

    Returns
    -------
    genes : (L x G) pandas.DataFrame
        Microarray expression data
    """
    # get combination function
    metric = utils.check_metric(metric)

    # take median of samples within an ROI to avoid outliers
    gene_by_label = np.zeros((len(roi_labels), microarray.shape[1]))
    subj_microarray = np.array(microarray)
    for nn in range(len(roi_labels)):
        temp = subj_microarray[roi_labels.iloc[nn].dropna().astype(int)]
        gene_by_label[nn, :] = metric(temp, 0)

    return pd.DataFrame(gene_by_label, columns=microarray.columns,
                        index=roi_labels.index)


def get_expression_data(files, atlas, atlas_info=None,
                        metric='mean', tolerance=3, use_centroids=False,
                        return_counts=False, dense=False, ibf_threshold=0.5):
    """
    Assigns microarray expression data in `files` to ROIs defined in `atlas`

    Parameters
    ----------
    files : dict
        Should have keys `microarray`, `probes`, and `anotation` pointing
        to filename (or list of filenames) of relevant files from Allen Brain
        Institute. Optimally obtained by calling `abagen.fetch_microarray()`.
    atlas : niimg-like object
        Parcel image, where each parcel should be identified with a unique
        integer ID
    atlas_info : str or pandas.DataFrame, optional
        Filepath to or pre-loaded dataframe containing information about
        `atlas`. Must jhave _at least_ columns 'id', 'hemisphere', and
        'structure' containing information mapping atlas IDs to hemisphere and
        broad structural class (i.e., "cortex", "subcortex", "cerebellum").
         Default: None
    metric : str or func, optional
        Mechanism by which to collapse across donors, if input `files` provides
        multiple microarray expression datasets. If str, should be in ['mean',
        'median']. If func, should be able to accept an `N`-dimensional input
        and the `axis` keyword argument and return an `N-1`-dimensional output.
        Default: 'mean'
    tolerance : int, optional
        Distance (in mm) that a sample must be from a ROI's boundary for it to
        be considered within that ROI. This is only used if the sample is not
        directly within a ROI. Default: 3
    use_centroids : bool, optional
        If no ROI is within `tolerance` of a sample, assign sample to ROI with
        closest centroid. Default: False
    return_counts : bool, optional
        Whether to return counts of how many samples were collapsed into each
        ROI, for each donor. Default: False
    dense : bool, optional
        Whether to return a dense microarray expression matrix by matching all
        ROIs in `atlas` to closest samples in `anotation`, instead of
        matching all samples to closes ROIs. Default: False
    ibf_threshold : [0, 1] float, optional
        Threshold for intensity-based filtering specifying the percentage of
        samples for which a probe must have signal above background noise in
        order to be retained for further consideration. Default: 0.5

    Returns
    -------
    expression : pandas.DataFrame
        Microarray expression averaged across samples within a given parcel and
        across probes within a given gene family
    labels : pandas.DataFrame
        Number of samples averaged into each ROI label, by donor (if multiple
        donors provided)
    """

    # coerce to Bunch in case simple dictionary was provided
    files = Bunch(**files)
    for key in ['microarray', 'probes', 'annotation']:
        if key not in files:
            raise KeyError('Provided `files` dictionary is missing {}. '
                           'Please check inputs.'.format(key))

    # load atlas_info, if porivded
    if atlas_info is not None:
        atlas_info = utils.check_atlas_info(atlas, atlas_info)

    # get combination functions
    metric = utils.check_metric(metric)
    if dense:
        label_func, group_func = label_rois, group_by_roi
    else:
        label_func, group_func = label_samples, group_by_label

    # get some info on the number of subjects, labels in `atlas_img`
    num_subj = len(files.microarray)
    all_labels = utils.get_unique_labels(atlas)

    # do some intensity-based filtering and DS selection on probes
    probes = processing.filter_probes(files.pacall,
                                      files.probes,
                                      threshold=ibf_threshold)
    probes = processing.get_stable_probes(files.microarray,
                                          files.annotation,
                                          probes)

    expression, labels = [], np.zeros((len(all_labels) + 1, num_subj))
    for subj in range(num_subj):
        # generate parcel labels for each sample
        annotation = processing.drop_mismatch_samples(files.annotation[subj],
                                                      files.ontology[subj])
        sample_labels = label_func(annotation, atlas, tolerance=tolerance)

        # get representative probes + samples from microarray data
        microarray = io.read_microarray(files.microarray[subj])
        sample_by_genes = microarray.loc[probes.index, annotation.index].T
        sample_by_genes.columns = probes.gene_symbol

        # aggregate samples within the same region and normalize data
        non_normalized = group_func(sample_by_genes, sample_labels,
                                    all_labels, metric=metric)
        expression += [processing.normalize_expression(non_normalized)]

        # get counts of samples collapsed into each ROI
        if dense:
            labs, counts = all_labels, sample_labels.count(axis=1)
        else:
            labs, counts = np.unique(sample_labels, return_counts=True)
        labels[labs, subj] = counts

    # aggregate over ROI across donors, if needed
    expression = pd.concat(expression).groupby('label').aggregate(metric)

    if return_counts:
        return expression, labels[1:]

    return expression
