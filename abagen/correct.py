# -*- coding: utf-8 -*-
"""
Functions for processing and correcting gene expression data
"""

import functools
import itertools
import warnings

import numpy as np
import pandas as pd
import scipy.stats as sstats
from scipy.spatial.distance import cdist

from . import images, utils


def _unpack_tuple(var):
    """
    Returns first entry of var if there is only one entry

    Parameters
    ----------
    var : tuple

    Returns
    -------
    var
    """

    if len(var) == 1:
        return var[0]
    return var


def _batch_correct(data):
    """
    Performs batch correction on `data`

    Parameters
    ----------
    data : list of (S, G) pandas.DataFrame
        Microarray expression data, where `S` is samples (or regions) and `G`
        is genes

    Returns
    -------
    residualized : list of (S, G) pandas.DataFrame
        Batch-corrected microarray expression data
    """

    if len(data) < 2:
        raise ValueError('Cannot perform batch correction with one batch.')

    # need to drop NaNs or the lstsq fit will choke
    data = [np.asarray(f) for f in data]
    data = [f[np.logical_not(np.all(np.isnan(f), axis=1))] for f in data]

    # generate donor label / "batch" array
    n_samp = [len(d) for d in data]
    batch = np.repeat(range(len(data)), n_samp)
    batch = np.column_stack([batch == f for f in np.unique(batch)]).astype(int)
    # add intercept for fit (but we won't regress this out)
    batch = np.column_stack([batch, np.ones(len(batch), dtype=int)])

    # fit least squares and residualize
    raw = np.row_stack(data)
    betas = np.linalg.lstsq(batch, raw, rcond=None)[0]
    resid = raw - (batch[:, :-1] @ betas[:-1])
    residualized = np.split(resid, np.cumsum(n_samp)[:-1])

    return residualized


def _rescale(data, low=0, high=1, axis=0):
    """
    Rescales `data` to range [`low`, `high`]

    Parameters
    ----------
    data : array_like
        Input data array to be rescaled
    low : float, optional
        Lower bound for rescaling. Default: 0
    high : float, optional
        Upper bound for rescaling. Default: 1
    axis : int, optional
        Axis of `data` to scale along. Default: 0

    Returns
    -------
    rescaled : np.ndarray
        Rescaled data
    """

    data = np.asanyarray(data)

    dmin = data.min(axis=axis, keepdims=True)
    dmax = data.max(axis=axis, keepdims=True)
    rescaled = ((data - dmin) / (dmax - dmin)) * (high - low) + low

    return rescaled


def _sigmoid(data, axis=0):
    """
    Normalizes `data` with a standard sigmoid function

    Parameters
    ----------
    data : array_like
        Input data array to be transformed
    axis : int, optional
        Axis of `data` to be normalized

    Returns
    -------
    normed : array_like
        Normalized input `data`
    """

    data = np.asanyarray(data)

    mean = np.mean(data, axis=axis, keepdims=True)
    std = np.std(data, ddof=1, axis=axis, keepdims=True)
    normed = 1 / (1 + np.exp(-(np.asarray(data) - mean) / std))

    return normed


def _scaledsig(data, low=0, high=1, axis=0):
    """
    Normalizes `data` with a scaled sigmoid function

    Parameters
    ----------
    data : array_like
        Input data array to be transformed
    low : float, optional
        Lower bound for rescaling. Default: 0
    high : float, optional
        Upper bound for rescaling. Default: 1
    axis : int, optional
        Axis of `data` to be normalized

    Returns
    -------
    normed : array_like
        Normalized input `data`
    """

    sig = _sigmoid(np.asanyarray(data), axis=axis)  # sigmoid transform
    normed = _rescale(sig, low=low, high=high, axis=axis)  # scales data

    return normed


def _rs(data, axis=0):
    """
    Normalizes `data` with a robust sigmoid function

    Parameters
    ----------
    data : array_like
        Input data array to be transformed
    axis : int, optional
        Axis of `data` to be normalized

    Returns
    -------
    normed : array_like
        Normalized input `data`
    """

    data = np.asanyarray(data)

    # calculate sigmoid normalization
    med = np.median(data, axis=axis, keepdims=True)
    iqr = sstats.iqr(data, axis=axis, scale='normal', keepdims=True)
    normed = 1 / (1 + np.exp(-(data - med) / iqr))

    return normed


def _srs(data, low=0, high=1, axis=0):
    """
    Normalizes `data` with a scaled robust sigmoid function

    Parameters
    ----------
    data : array_like
        Input data array to be transformed
    low : float, optional
        Lower bound for rescaling. Default: 0
    high : float, optional
        Upper bound for rescaling. Default: 1
    axis : int, optional
        Axis of `data` to be normalized

    Returns
    -------
    normed : array_like
        Normalized input `data`
    """

    srs = _rs(np.asanyarray(data), axis=axis)  # robust sigmoid transform
    normed = _rescale(srs, low=low, high=high, axis=axis)  # scales data

    return normed


def _scaledsig_qnt(data, low=0, high=1, axis=0):
    """
    Normalizes `data` with capped-quantile scaled sigmoid function

    Caps `data` at 5th and 95th percentiles and then uses `_scaledsig()`

    Parameters
    ----------
    data : array_like
        Input data array to be transformed
    low : float, optional
        Lower bound for rescaling. Default: 0
    high : float, optional
        Upper bound for rescaling. Default: 1
    axis : int, optional
        Axis of `data` to be normalized

    Returns
    -------
    normed : array_like
        Normalized input `data`
    """

    data = np.asanyarray(data)
    lq, hq = np.percentile(data, [5, 95], axis=axis, keepdims=True)
    mdata = np.ma.masked_array(data, np.logical_or(data > hq, data < lq))
    normed = np.asarray(_scaledsig(mdata, low=low, high=high, axis=axis))

    return normed


def _mixedsig(data, low=0, high=1, axis=0):
    """
    Uses `_scaledsig()` if IQR is 0; otherwise, uses `_srs()`

    Parameters
    ----------
    data : array_like
        Input data array to be transformed
    low : float, optional
        Lower bound for rescaling. Default: 0
    high : float, optional
        Upper bound for rescaling. Default: 1
    axis : int, optional
        Axis of `data` to be normalized

    Returns
    -------
    normed : array_like
        Normalized input `data`
    """

    data = np.asanyarray(data)

    iqr = sstats.iqr(data, axis=axis, scale='normal')
    mask = iqr == 0

    normed = np.zeros_like(data)
    normed[:, mask] = _scaledsig(data[:, mask])
    normed[:, ~mask] = _srs(data[:, ~mask])

    # constant columns set to 0
    return np.nan_to_num(normed)


NORMALIZATION_METHODS = dict(
    sig=functools.partial(_sigmoid, axis=0),
    scaled_sig=functools.partial(_scaledsig, low=0, high=1, axis=0),
    scaled_sig_qnt=functools.partial(_scaledsig_qnt, low=0, high=1, axis=0),
    mixed_sig=functools.partial(_mixedsig, low=0, high=1, axis=0),
    rs=functools.partial(_rs, axis=0),
    srs=functools.partial(_srs, low=0, high=1, axis=0),
    center=lambda x: x - x.mean(axis=0, keepdims=True),
    zscore=functools.partial(sstats.zscore, axis=0, ddof=1),
    minmax=functools.partial(_rescale, low=0, high=1, axis=0)
)
# aliases
for alias, orig in [('demean', 'center'),
                    ('rsig', 'rs'),
                    ('robust_sigmoid', 'rs'),
                    ('scaled_rsig', 'srs'),
                    ('scaled_robust_sigmoid', 'srs'),
                    ('sigmoid', 'sig'),
                    ('scaled_sigmoid', 'scaled_sig'),
                    ('scaled_sigmoid_quantiles', 'scaled_sig_qnt'),
                    ('mixed_sigmoid', 'mixed_sig')]:
    NORMALIZATION_METHODS[alias] = NORMALIZATION_METHODS[orig]


def normalize_expression(expression, norm='srs', structures=None,
                         ignore_warn=False):
    """
    Performs normalization on `expression` data

    Parameters
    ----------
    expression : list of (S, G) pandas.DataFrame
        Microarray expression data to be normalized, where `S` is samples (or
        regions) and `G` is genes
    norm : str, optional
        Function with which to normalize expression data. See Notes for more
        information on options. Default: 'scaled_robust_sigmoid'
    structures : list of (S,) pandas.DataFrame
        Structural designations of `S` samples (or regions) in `expression`.
        Index of provided data frames should be identical to `expression` and
        must have at least column 'structure'. If provided, normalization will
        be performed separately for each distinct structural class. Default:
        None
    ignore_warn : bool, optional
        Whether to suppress potential warnings raised by normalization.
        Default: False

    Returns
    -------
    normalized : list of (S, G) pandas.DataFrame
        Data from `expression` normalized separately for each gene

    Notes
    -----
    The following methods can be used for normalizing gene expression values
    for each donor (adapted from [PC2]_):

    1. ``norm='center'``

    Removes the mean of data in each column. Aliased to 'demean'

    2. ``norm='zscore'``

    Applies a basic z-score (subtract mean, divide by standard deviation) to
    each column; uses degrees of freedom equal to one for standard deviation

    3. ``norm='minmax'``

    Scales data in each column to the unit normal (i.e., range 0-1)

    4. ``norm='sigmoid'``

    Applies a sigmoidal transform function to normalize data in each column.
    Aliased to 'sig'

    5. ``norm='scaled_sigmoid'``

    Combines 'sigmoid' and 'minmax'. Aliased to 'scaled_sig'

    6. ``norm='scaled_sigmoid_quantiles'``

    Caps input data at the 5th and 95th percentiles before performing the
    'scaled_sigmoid' transform. Aliased to 'scaled_sig_qnt'

    7. ``norm='robust_sigmoid'``

    Uses a robust sigmoid function ([PC1]_) to normalize data in each column.
    Aliased to 'rs' and 'rsig'

    8. ``norm='scaled_robust_sigmoid'``

    Combines 'robust_sigmoid' and 'minmax'. Aliased to 'srs' and 'scaled_rsig'

    9. ``norm='mixed_sigmoid'``

    Uses 'scaled_sigmoid' transform for columns where the IQR is 0; otherwise,
    uses the 'scaled_robust_sigmoid' transform. Aliased to 'mixed_sig'

    10. ``norm='batch'``

    Uses a linear model to remove donor effects from data. Differs from other
    methods in that all donors are simultaneously fit to the same model and
    data are residualized based on estimated betas. Linear model includes the
    intercept but it is not removed during residualization

    References
    ----------
    .. [PC1] Fulcher, B. D., & Fornito, A. (2016). A transcriptional signature
       of hub connectivity in the mouse connectome. Proceedings of the National
       Academy of Sciences, 113(5), 1435-1440.
    .. [PC2] Fulcher, B. D., Little, M. A., & Jones, N. S. (2013). Highly
       comparative time-series analysis: the empirical structure of time series
       and their methods. Journal of the Royal Society Interface, 10(83),
       20130048
    """

    try:
        normfunc = NORMALIZATION_METHODS[norm]
    except KeyError:
        raise ValueError('Provided value for `norm` not recognized. Must be '
                         f'one of {list(NORMALIZATION_METHODS)}. Received: '
                         f'{norm}')
    kwargs = dict(all='ignore') if ignore_warn else {}

    # FIXME: I hate having to do this...
    if isinstance(expression, pd.DataFrame):
        expression = [expression]

    if structures is not None:
        if isinstance(structures, pd.DataFrame):
            structures = [structures]
        if any(len(s) != len(e) for s, e in zip(structures, expression)):
            raise ValueError('Length of `structures` class designations '
                             'differs from provided `expression`.')
    else:
        structures = [
            pd.DataFrame(dict(structure='same'), index=exp.index)
            for exp in expression
        ]

    if norm == 'batch':
        corrected = _batch_correct(expression)

    normexp = []
    for n, exp in enumerate(expression):
        normalized = pd.DataFrame(np.nan, columns=exp.columns, index=exp.index)
        notna = np.logical_not(exp.isna().all(axis=1))
        gb = structures[n].groupby('structure')
        for grp, samples in gb.groups.items():
            idx = samples[np.asarray(notna.loc[samples])]
            data = np.asarray(exp.loc[idx])
            # normalize the data (however was specified)
            with np.errstate(**kwargs):
                normed = normfunc(data) if norm != 'batch' else corrected[n]
            normalized.loc[idx] = normed

        normexp.append(normalized)

    return _unpack_tuple(normexp)


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
        Correlated gene expression array, where `R` is the number of regions,
        as generated with e.g., `numpy.corrcoef(expression)`.
    atlas : niimg-like object
        A parcellation image in MNI space, where each parcel is identified by a
        unique integer ID
    atlas_info : str or pandas.DataFrame, optional
        Filepath to or pre-loaded dataframe containing information about
        `atlas`. Must have at least columns 'id', 'hemisphere', and 'structure'
        containing information mapping atlas IDs to hemisphere (i.e, "L", "R")
        and broad structural class (i.e., "cortex", "subcortex", "cerebellum").
        Default: None
    labels : (N,) array_like, optional
        If only a subset `N` of the ROIs in `atlas` were used to generate the
        `coexpression` array this array should specify which to consider. Not
        specifying this may cause a ValueError if `atlas` and `atlas_info` do
        not match. Default: None

    Returns
    -------
    residualized : (R, R) numpy.ndarray
        Provided `coexpression` data residualized against spatial distance
        between region pairs
    """

    atlas = images.check_atlas(atlas)

    # check atlas + coexpression make sense
    if labels is None:
        labels = atlas.labels
    if len(labels) != len(coexpression):
        raise ValueError(f'Provided IDs for {labels.shape} parcels are a '
                         'different length than provided coexpression matrix '
                         f'of size {coexpression.shape}. Please confirm '
                         'inputs and try again.')

    # load atlas_info, if provided
    if atlas_info is not None:
        atlas_info = images.check_atlas_info(atlas_info, labels)

    # check that provided coexpression array is symmetric
    if not np.allclose(coexpression, coexpression.T, atol=1e-10):
        raise ValueError('Provided coexpression matrix is not symmetric')

    # we'll do basic Euclidean distance correction for now
    # TODO: implement gray matter volume / cortical surface path distance
    if labels is None:
        labels = atlas.labels
    centroids = np.row_stack([atlas.centroids[lab] for lab in labels])
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


def keep_stable_genes(expression, threshold=0.9, percentile=True, rank=True,
                      return_stability=False):
    """
    Removes genes in `expression` with differential stability < `threshold`

    Calculates the similarity of gene expression across brain regions for every
    pair of donors in `expression`. Similarity is averaged across donor pairs
    and genes whose mean similarity falls below `threshold` are removed.

    Parameters
    ----------
    expression : list of (R, G) pandas.DataFrame
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
    return_stability : bool, optional
        Whether to return stability estimates for each gene in addition to
        expression data. Default: False

    Returns
    -------
    expression : list of (R, Gr) pandas.DataFrame
        Microarray expression for `R` regions across `Gr` genes, where `Gr` is
        the number of retained genes
    stability : (G,) numpy.ndarray
        Stability (average correlation) of each gene across pairs of donors.
        Only returned if ``return_stability=True``
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

    # average similarity across donors (ignore NaNs)
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=RuntimeWarning,
                                message='Mean of empty slice')
        gene_corrs = np.nan_to_num(np.nanmean(gene_corrs, axis=1))

    # calculate absolute threshold if percentile is desired
    if percentile:
        threshold = np.percentile(gene_corrs, threshold * 100)
    keep_genes = gene_corrs >= threshold
    expression = [e.iloc[:, keep_genes] for e in expression]

    if return_stability:
        return expression, gene_corrs

    return expression
