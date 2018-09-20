# -*- coding: utf-8 -*-
"""
Functions for post-processing region x gene expression data
"""

import itertools
from nilearn._utils import check_niimg_3d
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from abagen import utils


def remove_distance(expression, atlas, atlas_info=None):
    """
    Corrects for distance-dependent correlation effects in `expression`

    Regresses Euclidean distance between regions in `atlas` from correlated
    gene expression array generated from `expression`. If `atlas_info` is
    provided different connection types (e.g., cortex-cortex, cortex-subcortex,
    subcortex-subcortex) will be residualized independently.

    Parameters
    ----------
    expression : (R x G) :class:`pandas.DataFrame`
        Microarray expression for `R` regions in `atlas` for `G` genes,
        aggregated across donors.
    atlas : niimg-like object
        A parcellation image in MNI space, where each parcel is identified by a
        unique integer ID
    atlas_info : str or :class:`pandas.DataFrame`, optional
        Filepath to or pre-loaded dataframe containing information about
        `atlas`. Must have at least columns 'id', 'hemisphere', and 'structure'
        containing information mapping atlas IDs to hemisphere (i.e, "L", "R")
        and broad structural class (i.e., "cortex", "subcortex", "cerebellum").
        Default: None

    Returns
    -------
    corrgene : (R x R) :class:`numpy.ndarray`
        Correlated gene `expression` data for `R` regions in `atlas`,
        residualized against the spatial distances between region pairs
    """

    # load atlas_info, if provided
    atlas = check_niimg_3d(atlas)
    if atlas_info is not None:
        atlas_info = utils.check_atlas_info(atlas, atlas_info)

    # check expression data and make correlation matrix
    if not isinstance(expression, pd.DataFrame):
        raise TypeError('Provided `expression` data must be type pd.DataFrame '
                        'not {}'.format(type(expression)))
    genecorr = np.corrcoef(expression.get_values())

    # we'll do basic Euclidean distance correction for now
    # TODO: implement gray matter volume / cortical surface path distance
    centroids = utils.get_centroids(atlas, labels_of_interest=expression.index)
    dist = cdist(centroids, centroids, metric='euclidean')

    corr_resid = np.zeros_like(genecorr)
    triu_inds = np.triu_indices_from(genecorr, k=1)
    # if no atlas_info, just residualize all correlations against distance
    if atlas_info is None:
        corr_resid[triu_inds] = _resid_dist(genecorr[triu_inds],
                                            dist[triu_inds])
    # otherwise, we can residualize the different connection types separately
    else:
        triu_inds = np.ravel_multi_index(triu_inds, genecorr.shape)
        genecorr, dist = genecorr.ravel(), dist.ravel()
        # FIXME: residualizing cortex/subcortex and subcortex/cortex separately
        for src, tar in itertools.product(['cortex', 'subcortex'], repeat=2):
            # get indices of sources and targets
            sources = atlas_info.query('structure == "{}"'.format(src)).index
            targets = atlas_info.query('structure == "{}"'.format(tar)).index
            inds = np.ix_(np.where(expression.index.isin(sources))[0],
                          np.where(expression.index.isin(targets))[0])
            # find intersection of source / target indices + upper triangle
            inds = np.intersect1d(triu_inds,
                                  np.ravel_multi_index(inds, corr_resid.shape))
            back = np.unravel_index(inds, corr_resid.shape)
            corr_resid[back] = _resid_dist(genecorr[inds], dist[inds])

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
    distance = np.column_stack((dv, np.ones_like(dv)))
    betas, *rest = np.linalg.lstsq(distance, iv[:, np.newaxis], rcond=None)
    residuals = iv[:, np.newaxis] - (distance @ betas)

    return residuals.squeeze()


def gene_stability(expression, threshold=0.8):
    """
    Removes genes with differential stability < `threshold` across donors

    Parameters
    ----------
    expression : list of (R x G) :class:`pandas.DataFrame`
        Microarray expression for `R` regions in across `G` genes for each
        donor
    threshold : [0, 1] float, optional
        Minimum required correlation of gene expression across regions averaged
        across donors for a gene to be retained

    Returns
    -------
    expression : list of (R x Gr) :class:`pandas.DataFrame`
        Microarray expression for `R` regions across `Gr` genes, where `Gr` is
        the number of retained genes
    """

    # get number of donors and number of genes
    num_subj = len(expression)
    num_gene = expression[0].shape[-1]

    # rank data for all subjects
    ranked = [e.rank() for e in expression]

    # get correlation of gene expression across regions for all donor pairs
    gene_corrs = np.zeros((num_gene, sum(range(num_subj))))
    for n, (s1, s2) in enumerate(itertools.combinations(range(num_subj), 2)):
        gene_corrs[:, n] = utils.efficient_corr(ranked[s1], ranked[s2])

    # calculate mean correlation for each gene across donor pairs and threshold
    keep_genes = gene_corrs.mean(axis=1) > threshold
    expression = [e.iloc[:, keep_genes] for e in expression]

    return expression
