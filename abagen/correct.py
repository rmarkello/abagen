# -*- coding: utf-8 -*-
"""
Functions for post-processing region x gene expression data
"""

import itertools
from nilearn._utils import check_niimg_3d
import numpy as np
from scipy.spatial.distance import cdist
from sklearn.utils.validation import check_symmetric
from abagen import utils


def remove_distance(coexpression, atlas, atlas_info=None, labels=None):
    """
    Corrects for distance-dependent correlation effects in `coexpression`

    Regresses Euclidean distance between regions in `atlas` from correlated
    gene expression array `coexpression`. If `atlas_info` is provided different
    connection types (e.g., cortex-cortex, cortex-subcortex, subcortex-
    subcortex) will be residualized independently.

    Parameters
    ----------
    coexpression : (R x R) array_like
        Correlate gene expression array, where `R` is the number of regions, as
        generated with e.g., `numpy.corrcoef(expression)`.
    atlas : niimg-like object
        A parcellation image in MNI space, where each parcel is identified by a
        unique integer ID
    atlas_info : str or :class:`pandas.DataFrame`, optional
        Filepath to or pre-loaded dataframe containing information about
        `atlas`. Must have at least columns 'id', 'hemisphere', and 'structure'
        containing information mapping atlas IDs to hemisphere (i.e, "L", "R")
        and broad structural class (i.e., "cortex", "subcortex", "cerebellum").
        Default: None
    labels : (N,) array_like, optional
        If only a subset `N` of the ROIs in `atlas` were used to generate the
        `coexpression` array this array should specify which. Default: None

    Returns
    -------
    residualized : (R x R) :class:`numpy.ndarray`
        Provided `coexpression` data residualized against spatial distance
         between region pairs
    """

    # load atlas_info, if provided
    atlas = check_niimg_3d(atlas)
    if atlas_info is not None:
        atlas_info = utils.check_atlas_info(atlas, atlas_info,
                                            labels_of_interest=labels)

    # check that provided coexpression array is symmetric
    check_symmetric(coexpression, raise_exception=True)

    # we'll do basic Euclidean distance correction for now
    # TODO: implement gray matter volume / cortical surface path distance
    centroids = utils.get_centroids(atlas, labels_of_interest=labels)
    dist = cdist(centroids, centroids, metric='euclidean')

    corr_resid = np.zeros_like(coexpression)
    triu_inds = np.triu_indices_from(coexpression, k=1)
    # if no atlas_info, just residualize all correlations against distance
    if atlas_info is None:
        corr_resid[triu_inds] = _resid_dist(coexpression[triu_inds],
                                            dist[triu_inds])
    # otherwise, we can residualize the different connection types separately
    else:
        triu_inds = np.ravel_multi_index(triu_inds, corr_resid.shape)
        coexpression, dist = coexpression.ravel(), dist.ravel()
        types = ['cortex', 'subcortex']
        for src, tar in itertools.combinations_with_replacement(types, 2):
            # get indices of sources and targets
            sources = np.where(atlas_info.structure == src)[0]
            targets = np.where(atlas_info.structure == tar)[0]
            inds = np.ravel_multi_index(np.ix_(sources, targets),
                                        corr_resid.shape)
            if src != tar:  # e.g., cortex + subcortex
                rev = np.ravel_multi_index(np.ix_(targets, sources),
                                           corr_resid.shape)
                inds = np.append(inds.ravel(), rev.ravel())
            # find intersection of source / target indices + upper triangle
            inds = np.intersect1d(triu_inds, inds)
            back = np.unravel_index(inds, corr_resid.shape)
            # residualize
            corr_resid[back] = _resid_dist(coexpression[inds], dist[inds])

    corr_resid = (corr_resid + corr_resid.T + np.eye(len(corr_resid)))

    return corr_resid


def _resid_dist(dv, iv):
    """
    Calculates residuals of `dv` after controlling for `iv`

    Parameters
    ----------
    dv : array_like
        Dependent variable
    iv : array_like
        Independent variable; removed from `dv`

    Returns
    -------
    residuals : array_like
        Residuals of `dv` after controlling for `iv`
    """
    dv = dv.squeeze()
    distance = np.column_stack((iv, np.ones_like(iv)))
    betas, *rest = np.linalg.lstsq(distance, dv, rcond=None)
    residuals = dv - (distance @ betas)

    return residuals.squeeze()


def keep_stable_genes(expression, threshold=0.9, percentile=True, rank=True):
    """
    Removes genes in `expression` with differential stability < `threshold`

    Calculates the similarity of gene expression across brain regions for every
    pair of donors in `expression`. Similarity is averaged across donor pairs
    and genes whose mean similarity falls below `threshold` are removed.

    Parameters
    ----------
    expression : list of (R x G) :class:`pandas.DataFrame`
        Where each entry is the microarray expression of `R` regions across `G`
        genes for a given donor
    threshold : [0, 1] float, optional
        Minimum required average similarity (e.g, correlation) across donors
        for a gene to be retained. Default: 0.1
    percentile : bool, optional
        Whether to treat `threshold` as a percentile instead of an absolute
        cutoff. For example, `threshold=0.9` and `percentile=True` would
        retain only those genes with a differential stability in the top 10% of
        all genes, whereas `percentile=False` would retain only those genes
        with differential stability > 0.9. Default: True
    rank : bool, optional
        Whether to calculate similarity as Spearman correlation instead of
        Pearson correlation. Default: True

    Returns
    -------
    expression : list of (R x Gr) :class:`pandas.DataFrame`
        Microarray expression for `R` regions across `Gr` genes, where `Gr` is
        the number of retained genes
    """

    # get number of donors and number of genes
    num_subj = len(expression)
    num_gene = expression[0].shape[-1]

    # rank data, if necessary
    for_corr = expression if not rank else [e.rank() for e in expression]

    # get correlation of gene expression across regions for all donor pairs
    gene_corrs = np.zeros((num_gene, sum(range(num_subj))))
    for n, (s1, s2) in enumerate(itertools.combinations(range(num_subj), 2)):
        regions = np.intersect1d(for_corr[s1].dropna(axis=0, how='all').index,
                                 for_corr[s2].dropna(axis=0, how='all').index)
        gene_corrs[:, n] = utils.efficient_corr(for_corr[s1].loc[regions],
                                                for_corr[s2].loc[regions])

    # average similarity across donors
    gene_corrs = gene_corrs.mean(axis=1)
    # calculate absolute threshold if percentile is desired
    if percentile:
        threshold = np.percentile(gene_corrs, threshold * 100)
    keep_genes = gene_corrs > threshold
    expression = [e.iloc[:, keep_genes] for e in expression]

    return expression
