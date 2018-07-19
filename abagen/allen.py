# -*- coding: utf-8 -*-

from nilearn._utils import check_niimg_3d
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from sklearn.utils import Bunch

from abagen import io, utils


AGG_FUNCS = dict(mean=np.mean,
                 median=np.median)


def assign_sample(sample, label_image, tolerance=3):
    """
    Determines what ROI ``sample`` belongs to in ``label_image``

    Parameters
    ----------
    sample : (3 x 1) array_like
        Coordinates (ijk) of microarray probe in ``label_image`` space
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
        ROI label of ``sample``
    """

    # pull relevant info from label_image
    label_data = check_niimg_3d(label_image).get_data()

    # expand provided coordinates to include those w/i `tolerance` of `coords`
    # set a hard euclidean distance limit to account for different voxel sizes
    coords = utils.expand_roi(sample, dilation=tolerance, return_array=True)
    coords = coords[cdist(sample.T, coords).squeeze() < tolerance]

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
    indmax = np.where(counts == counts.max())[0]
    if indmax.size == 1:
        return labels[indmax[0]]

    # if two or more parcels tied for neighboring frequency, use ROI
    # with closest centroid to `coords`
    centroids = utils.get_centroids(label_image, labels)
    return labels[utils.closest_centroid(sample, centroids)]


def label_samples(annotation, label_image, tolerance=3, use_centroids=False):
    """
    Matches all samples in ``annotation`` to closest ROIs in ``label_image``

    Attempts to place each sample provided in ``annotation`` into a labelled
    ROI in ``label_image``, where the latter is a 3D niimg-like object that
    contains ROI data, each with a unique integer ID.

    The function tries to best match the microarray sample to the ROIs by:
    1. Determining if the sample falls directly within a ROI.
    2. Checking to see if there are nearby ROIs by slowly expanding radius
       around sample (up to radius = ``tolerance``); if there are multiple
       nearby ROIs, determining which ROI is closer (i.e., by centroid of
       the ROI).
    3. If ``user_centroids`` is set, find the ROI with the closest centroid

    If there is still no ROI for a given sample, it is assigned as "NA" and
    should likely be ignored in further analyses.

    Parameters
    ----------
    annotation : (S x 13) pd.core.frame.DataFrame
        Annotation information, where ``S`` is samples
    label_image : niimg-like object
        ROI image, where each ROI should be identified with a unique integer ID
    tolerance : int, optional
        Distance that sample must be from ROI for it to be considered within
        a ROI. This is only used if the sample is not directly inside a ROI.
        Default: 3mm
    use_centroids : bool, optional
        If no ROI is within ``tolerance`` of a sample, assign sample to ROI
        with closest centroid. Default: False

    Returns
    -------
    labels : (S x 1) pd.core.frame.DataFrame
        Dataframe with ROI labels
    """

    # read annotation file, if provided
    if isinstance(annotation, str):
        annotation = io.read_sampleannot(annotation)

    # get image data
    label_image = check_niimg_3d(label_image)
    label_data, affine_trans = label_image.get_data(), label_image.affine[:-1]

    # if we're going to use centroids, calculate them ahead of time
    if use_centroids:
        all_labels = utils.get_unique_labels(label_image)
        centroids = utils.get_centroids(label_image, all_labels)

    # grab xyz coordinates for microarray samples and convert to ijk
    g_xyz = annotation[['mni_x', 'mni_y', 'mni_z']].get_values()
    g_ijk = np.floor(utils.xyz_to_ijk(g_xyz, affine_trans)).T.astype(int)

    # get labels for all ijk values
    labelled_samples = label_data[g_ijk[:, 0], g_ijk[:, 1], g_ijk[:, 2]]

    # if coordinates aren't within the parcel, check for neighboring parcels
    # and slowly increase the radius around parcel up to ``tolerance`` to try
    # and find nearby parcels. if still no nearby parcel then ignore probe
    for n, label in enumerate(labelled_samples):
        sample, tol = g_ijk[n, None].T, 1
        while label == 0 and tol <= tolerance:
            label = assign_sample(sample, label_image, tolerance=tol)
            tol += 1
        if label == 0 and use_centroids:
            label = all_labels[utils.closest_centroid(sample, centroids)]

        labelled_samples[n] = label

    # return DataFrame for ease of use
    return pd.DataFrame(labelled_samples,  columns=['label'], dtype=int)


def group_by_gene(microarray, probes):
    """
    Average over probes in ``microarray``, grouping by genes

    Parameters
    ----------
    microarray : (P x S) pd.core.frame.DataFrame
        Micoarray expression data, where ``P`` is probes and ``S`` is samples
    probes : (P x 6) pd.core.frame.DataFrame
        Probe information, where ``P`` is probes

    Returns
    -------
    microarray_by_gene : (S x G) pd.core.frame.DataFrame
        ``microarray``, where ``G`` is the number of unique genes
    """

    microarray_by_gene = (microarray.merge(probes[['gene_symbol']],
                                           left_index=True,
                                           right_index=True)
                                    .groupby(['gene_symbol'])
                                    .mean()
                                    .drop(['na'])
                                    .T
                                    .reset_index(drop=True))

    return microarray_by_gene


def group_by_label(microarray, sample_labels, labels=None, metric='mean'):
    """
    Averages expression data in ``microarray`` over samples with same label

    Parameters
    ----------
    microarray : (S x G) pd.core.frame.DataFrame
        Microarray expression data, where ``S`` is samples and ``G`` is genes
    sample_labels : (S x 1) pd.core.frame.DataFrame
        ROI labels for ``S`` samples, as returned by e.g., ``label_samples()``
    labels : (L,) array_like, optional
        All possible labels for parcellation scheme (to account for possibility
        that some parcels have NO expression data). Default: None
    metric : str or func
        Mechanism by which to collapse across samples within an ROI. If str,
        should be in ['mean', 'median'].

    Returns
    -------
    genes : (L x G) pd.core.frame.DataFrame
        Microarray expression data
    """

    # get combination function
    if isinstance(metric, str):
        try:
            metric = AGG_FUNCS[metric]
        except KeyError:
            raise ValueError('Provided metric {0} is not valid. If supplied'
                             'as string, metric must be in {1}.'
                             .format(metric, list(AGG_FUNCS.keys())))

    # get missing labels
    if labels is not None:
        labels = pd.DataFrame(index=np.setdiff1d(labels, sample_labels))

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
    Matches all ROIs in ``label_image`` to closest samples in ``annotation``

    Attempts to place each sample provided in ``annotation`` into a labelled
    ROI in ``label_image``, where the latter is a 3D niimg-like object that
    contains ROI data, each with a unique integer ID.

    The function tries to best match the microarray sample to the ROIs by:
    1. Determining which samples falls directly within a ROI.
    2. If there is no sample assigned to an ROI, find the closest sample to the
       centroid of the ROI

    Parameters
    ----------
    annotation : (S x 13) pd.core.frame.DataFrame
        Annotation information, where ``S`` is samples
    label_image : niimg-like object
        ROI image, where each ROI should be identified with a unique integer ID
    tolerance : int, optional
        Not uses in this implementation.

    Returns
    -------
    labels : (L x ...) pd.core.frame.DataFrame
        Dataframe with lists of sample labels for each ROI
    """

    # read annotation file, if provided
    if isinstance(annotation, str):
        annotation = io.read_sampleannot(annotation)

    # get image data
    label_image = check_niimg_3d(label_image)
    label_data, affine_trans = label_image.get_data(), label_image.affine[:-1]

    # calculate the centroids
    all_labels = utils.get_unique_labels(label_image)
    centroids = utils.get_centroids(label_image, all_labels)

    # grab xyz coordinates for microarray samples and convert to ijk
    g_xyz = annotation[['mni_x', 'mni_y', 'mni_z']].get_values()
    g_ijk = np.floor(utils.xyz_to_ijk(g_xyz, affine_trans)).T.astype(int)

    # get labels for all ijk values
    labelled_samples = label_data[g_ijk[:, 0], g_ijk[:, 1], g_ijk[:, 2]]

    # make a list of samples within each ROI
    # if no sample is assigned to an ROI, find closest sample to ROI centroid
    listlabelled_samples = labelled_samples.tolist()
    labelled_rois = []
    for nn in range(len(all_labels)):
        indices = [i for i, x in enumerate(listlabelled_samples) if x == nn+1]
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
    Averages expression data in ``microarray`` over samples within same ROI

    Parameters
    ----------
    microarray : (S x G) pd.core.frame.DataFrame
        Microarray expression data, where ``S`` is samples and ``G`` is genes
    roi_labels : (L x ...) pd.core.frame.DataFrame
        Lists of sample labels with varying lengths for ``L`` ROIs, as returned
        by e.g., ``label_rois()``
    labels : (L,) array_like, optional
        All possible labels for parcellation scheme (to account for possibility
        that some parcels have NO expression data). Default: None
    metric : str or func
        Mechanism by which to collapse across samples within an ROI. If str,
        should be in ['mean', 'median'].

    Returns
    -------
    genes : (L x G) pd.core.frame.DataFrame
        Microarray expression data
    """
    # get combination function
    if isinstance(metric, str):
        try:
            metric = AGG_FUNCS[metric]
        except KeyError:
            raise ValueError('Provided metric {0} is not valid. If supplied'
                             'as string, metric must be in {1}.'
                             .format(metric, list(AGG_FUNCS.keys())))

    # take median of samples within an ROI to avoid outliers
    gene_by_label = np.zeros((len(roi_labels), microarray.shape[1]))
    subj_microarray = np.array(microarray)
    for nn in range(len(roi_labels)):
        temp = subj_microarray[roi_labels.iloc[nn].dropna().astype(int)]
        gene_by_label[nn, :] = metric(temp, 0)

    return pd.DataFrame(gene_by_label, columns=microarray.columns,
                        index=roi_labels.index)


def get_expression_data(files, label_image, metric='mean', tolerance=3,
                        use_centroids=False, return_counts=False, dense=False):
    """
    Assigns microarray expression data in ``files`` to ROIs in ``label_image``

    Parameters
    ----------
    files : dict
        Should have keys ``microarray``, ``probes``, and ``anotation`` pointing
        to filename (or list of filenames) of relevant files from Allen Brain
        Institute. Optimally obtained by calling ``abagen.fetch_microarray()``.
    label_image : niimg-like object
        Parcel image, where each parcel should be identified with a unique
        integer ID
    metric : str or func
        Mechanism by which to collapse across donors, if input ``files``
        provides multiple microarray expression datasets. If str, should be in
        ['mean', 'median']. Default: 'mean'
    tolerance : int, optional
        Distance (in mm) that a sample must be from a ROI's boundary for it to
        be considered within that ROI. This is only used if the sample is not
        directly within a ROI. Default: 3
    use_centroids : bool, optional
        If no ROI is within ``tolerance`` of a sample, assign sample to ROI
        with closest centroid. Default: False
    return_counts : bool, optional
        Whether to return counts of how many samples were collapsed into each
        ROI, for each donor. Default: False
    dense : bool, optional
        Whether to return a dense microarray expression matrix by matching all
        ROIs in ``label_image`` to closest samples in ``anotation``, instead of
        matching all samples to closes ROIs. Default: False

    Returns
    -------
    expression : pd.core.frame.DataFrame
        Microarray expression averaged across samples within a given parcel and
        across probes within a given gene family
    labels : pd.core.frame.DataFrame
        Number of samples averaged into each ROI label, by donor (if multiple
        donors provided)
    """

    # coerce to Bunch in case simple dictionary was provided
    files = Bunch(**files)
    for key in ['microarray', 'probes', 'annotation']:
        if key not in files:
            raise KeyError('Provided ``files`` dictionary is missing {}. '
                           'Please check inputs.'.format(key))

    # get combination function
    if isinstance(metric, str):
        try:
            metric = AGG_FUNCS[metric]
        except KeyError:
            raise ValueError('Provided metric {0} is not valid. If supplied'
                             'as string, metric must be in {1}.'
                             .format(metric, list(AGG_FUNCS.keys())))

    # get some info on the number of subjects, labels in `label_image`
    num_subj = len(files.microarray)
    all_labels = utils.get_unique_labels(label_image)

    # empty lists and arrays to hold the expression information
    expression, labels = [], np.zeros((len(all_labels) + 1, num_subj))
    for subj in range(num_subj):
        # average microarray expression across all probes w/i same gene
        microarray, probes = files.microarray[subj], files.probes[subj]
        sample_by_genes = group_by_gene(io.read_microarray(microarray),
                                        io.read_probes(probes))

        # generate parcel labels for each microarray sample
        annotation = io.read_sampleannot(files.annotation[subj])

        if dense:
            label_func, group_func = label_rois, group_by_roi
        else:
            label_func, group_func = label_samples, group_by_label

        sample_labels = label_func(annotation, label_image,
                                   tolerance=tolerance)
        expression += [group_func(sample_by_genes, sample_labels,
                                  all_labels, metric=metric)]

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
