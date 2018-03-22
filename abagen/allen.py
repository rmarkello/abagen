# -*- coding: utf-8 -*-

import glob
import os
import numpy as np
import pandas as pd

from abagen import utils


def group_by_gene(microarray, info):
    """
    Average over probes in ``microarray``, grouping by gene ``info``

    Parameters
    ----------
    microarray : pd.core.frame.DataFrame
        From AIBS `MicroarrayExpression.csv`
    info : pd.core.frame.DataFrame
        From AIBS `Probes.csv`

    Returns
    -------
    expression : pd.core.frame.DataFrame
        ``microarray`` with probes averaged into gene families
    """

    microarray_by_gene = (pd.merge(microarray,
                                   info[['probe_id', 'gene_symbol']],
                                   on='probe_id')
                            .groupby(['gene_symbol'])
                            .mean())
    microarray_by_gene = microarray_by_gene.drop(['probe_id'], axis=1)

    return microarray_by_gene.T.drop(['na'], axis=1)


def assign_roi(sample, label_image, tolerance=3):
    """
    Determines what parcel ``sample`` belongs to in ``label_image``

    Parameters
    ----------
    sample : (3 x 1) array_like
        Coordinates (i, j, k) of microarray probe in ``label_image`` space
    label_image : niimg-like object
        Parcel image, where each parcel should be identified with a unique
        integer ID
    tolerance : int, optional
        Distance that sample must be from parcel for it to be considered
        within that parcel. This is only used if the sample is not directly
        within a parcel. Default: 3mm

    Returns
    -------
    label : int
        "Best fit" label for ``sample`` in ``label_image``
    """

    # pull relevant info from label_image
    label_data = label_image.get_data()

    # expand provided coordinates to include those w/i `tolerance` of `coords`
    coords = utils.expand_roi(sample, dilation=tolerance, return_array=True)
    possible_labels = label_data[coords[:, 0], coords[:, 1], coords[:, 2]]
    nz_labels = possible_labels[possible_labels.nonzero()]

    # determine unique labels and counts of those labels
    labels, counts = np.unique(nz_labels, return_counts=True)

    # if there is still nothing in the vicinity, we'll have to discard probe
    if labels.size == 0:
        return 0
    # if there is only one parcel in the vicinity, use that
    elif labels.size == 1:
        return labels[0]
    else:
        # if more than one parcel in the vicinity, return the most frequent
        indmax = np.where(counts == counts.max())[0]
        if indmax.size == 1:
            return labels[indmax[0]]
        # if two or more parcels tied for neighboring frequency, use parcel
        # with closest centroid to `coords`
        else:
            centroids = utils.get_centroids(label_image,
                                            labels,
                                            image_space=True)
            return labels[utils.closest_centroid(sample, centroids)]


def match_sample_to_parcel(annotation, label_image, tolerance=3):
    """
    Matches sample from ``annotation`` to closest parcel in ``label_image``

    Attempts to place each sample provided in ``annotation`` from the Allen
    Institute microarray data into a labelled parcel in ``label_image``, where
    ``label_image`` is a 3D Niimg-like object that contains parcel data, each
    with a unique integer ID.

    The function tries to best match the microarray sample to the parcels by:
    1. Determining if the sample falls directly within a parcel.
    2. Checking to see if there are nearby parcels by slowly expanding radius
       around sample (up to radius = ``tolerance``); if there are multiple
       nearby parcels, determining which parcel is closer (i.e., by centroid of
       the parcel).

    If there is still no parcel for a given sample, it is assigned as "NA" and
    should likely be ignored in further analyses.

    Parameters
    ----------
    annotation : pd.core.frame.DataFrame
        DataFrame containing data from Annotation.csv for one donor from AIBS,
        as obtained via `pd.read_csv('Annotation.csv')`
    label_image : niimg-like object
        Parcel image, where each parcel should be identified with a unique
        integer ID
    tolerance : int, optional
        Distance (in mm) that sample must be from parcel to be considered
        within that parcel. This is only used if the sample is not directly
        within a parcel. Default: 3mm

    Returns
    -------
    pd.core.frame.DataFrame
        Dataframe with parcel label; shape (probes x 1)
    """

    label_data, affine_trans = label_image.get_data(), label_image.affine[:-1]

    # grab xyz coordinates for microarray samples and convert to ijk
    g_xyz = annotation[['mni_x', 'mni_y', 'mni_z']].get_values()
    g_ijk = np.floor(utils.xyz_to_ijk(g_xyz, affine_trans)).T.astype(int)
    # TODO:

    # get labels for all ijk values
    labelled_samples = label_data[g_ijk[:, 0], g_ijk[:, 1], g_ijk[:, 2]]

    # go through every sample location
    for n, label in enumerate(labelled_samples):
        if label != 0:
            continue
        # if coordinates not within parcel, check for neighboring parcels
        # slowly increase radius around parcel (by 1mm each iteration)
        # up to ``tolerance`` to try and find nearby parcels. if still no
        # nearby parcel at ``tolerance``, then consider probe to be "off-brain"
        for tol in range(1, tolerance + 1):
            label = assign_roi(g_ijk[n, None].T,
                               label_image,
                               tolerance=tol)
            if label != 0:
                break
        labelled_samples[n] = label

    # convert probe information into DataFrame for ease of use
    labelled_samples = pd.DataFrame(labelled_samples, columns=['label'])

    return labelled_samples


def group_by_label(microarray, sample_labels, labels):
    """
    Averages expression data in `microarray` over samples within same parcel

    Parameters
    ----------
    microarray : (N x G) pd.core.frame.DataFrame
        Microarray expression data, where ``N`` is probes and ``G`` is genes
    sample_labels : (N x L) pd.core.frame.DataFrame
        Sample labels as returned by `match_probe_labels()`, where ``N`` is
        probes and ``L`` is labels
    labels : (L,) array_like
        All possible labels for parcellation scheme (to account for possibility
        that some parcels have NO expression data)

    Returns
    -------
    genes : (L x G) pd.core.frame.DataFrame
        Microarray expression data
    """

    # create empty DataFrame to store parcelled microexpression data
    gene_by_parcel = pd.DataFrame([], columns=microarray.columns)

    # iterate through all possible labels
    for label in labels:
        # find all probes with current label
        indices = sample_labels.query(f'label == {label}').index
        # take median of microexpression data across all probes with label
        label_expression = microarray.iloc[indices].median()
        label_expression.name = label
        gene_by_parcel = gene_by_parcel.append(label_expression)

    return gene_by_parcel


def get_expression_from_donor(donor_path, label_image, tolerance=3):
    """
    Parameters
    ----------
    donor_path : str
        Filepath to unzipped directory from AIBS Microarray Data. Directory
        should contain "MicroarrayExpression.csv", "Probes.csv", and
        "SampleAnnot.csv", at a minimum
    label_image : niimg-like object
        Parcel image, where each parcel should be identified with a unique
        integer ID
    tolerance : int, optional
        Distance (in mm) that a sample must be from a parcel's boundary for it
        to be considered *within* that parcel. This is only used if the sample
        is not directly within a parcel. Default: 3

    Returns
    -------
    parcel_by_gene : pandas.core.frame.DataFrame
        Microarray expression averaged across samples within a given parcel and
        across probes within a given gene family
    sample_labels : pandas.core.frame.DataFrame
        Parcel labels for each sample in original microarray data; for
        debugging purposes and to determine how many samples were averaged into
        each parcel
    """

    # load in annotation, probe, and microarray data from AIBS dataset
    annotation = pd.read_csv(glob.glob(os.path.join(donor_path,
                                                    'SampleAnnot.csv'))[0])
    probes = pd.read_csv(glob.glob(os.path.join(donor_path, 'Probes.csv'))[0])
    microarray = pd.read_csv(glob.glob(os.path.join(donor_path,
                                       'MicroarrayExpression.csv'))[0],
                             names=np.append(['probe_id'],
                                             annotation['well_id'].get_values()
                                             )
                             )

    # determine all possible parcel labels (for this scale)
    possible_labels = utils.get_unique_labels(label_image)
    # average microarray expression across all probes w/i same gene
    samples_by_genes = group_by_gene(microarray, probes)
    # generate parcel labels for each microarray sample
    sample_labels = match_sample_to_parcel(label_image,
                                           annotation,
                                           tolerance=tolerance)
    # average microarray expression across all samples w/i same parcel
    parcel_by_gene = group_by_label(samples_by_genes,
                                    sample_labels,
                                    possible_labels)

    return parcel_by_gene, sample_labels.get_values()
