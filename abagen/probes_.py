# -*- coding: utf-8 -*-
"""
Functions for annotating / selecting microarray probes
"""

import functools
import gzip
from io import StringIO
import itertools
import logging
from pkg_resources import resource_filename

import numpy as np
import pandas as pd
from scipy import stats as sstats

from . import datasets, io, utils

LGR = logging.getLogger('abagen')


def reannotate_probes(probes):
    """
    Replaces gene symbols in `probes` with reannotated data

    Uses annotations from [PR18]_ to replace probe annotations shipped with
    AHBA data. Any probes that were unable to be matched to a gene in
    reannotation procedure are not retained.

    Parameters
    ----------
    probes : str or pandas.DataFrame
        Probe file or loaded probe dataframe from Allen Brain Institute
        containing information on microarray probes

    Returns
    -------
    reannotated : pandas.DataFrame
        Provided probe information with updated gene symbols and Entrez IDs

    References
    ----------
    .. [PR18] Arnatkevic̆iūtė, A., Fulcher, B. D., & Fornito, A. (2019). A
       practical guide to linking brain-wide gene expression and neuroimaging
       data. NeuroImage, 189, 353-367.
    """

    LGR.info('Reannotating probes with information from Arnatkevic̆iūtė '
             'et al., 2019, NeuroImage')

    # load in reannotated probes
    reannot = resource_filename('abagen', 'data/reannotated.csv.gz')
    with gzip.open(reannot, 'r') as src:
        reannot = pd.read_csv(StringIO(src.read().decode('utf-8')))

    # merge reannotated with original, keeping only reannotated
    probes = io.read_probes(probes).reset_index()
    merged = pd.merge(reannot[['probe_name', 'gene_symbol', 'entrez_id']],
                      probes[['probe_name', 'probe_id']],
                      on='probe_name', how='left')

    # reset index as probe_id and sort
    reannotated = merged.set_index('probe_id') \
                        .sort_index() \
                        .dropna(subset=['entrez_id'])
    reannotated.loc[:, 'entrez_id'] = reannotated['entrez_id'].astype('int')

    return reannotated


def filter_probes(pacall, annotation, probes, threshold=0.5):
    """
    Performs intensity based filtering (IBF) of expression probes

    Uses binary indicator for expression levels in `pacall` to determine which
    probes have expression levels above background noise in `threshold` of
    samples across donors.

    Parameters
    ----------
    pacall : dict
        Dictionary where keys are donor IDs and values are filepaths to (or
        dataframes of) PACall.csv files from Allen Brain Institute
    annotation : dict
        Dictionary where keys are donor IDs and values are filepaths to (or
        dataframes of) SampleAnnot.csv files from Allen Brain Institute
    probes : str or pandas.DataFrame
        Filepath to Probes.csv or dataframe containing information on
        microarray probes that should be considered in filtering (probes not in
        this will be ignored)
    threshold : (0, 1) float, optional
        Threshold for filtering probes. Specifies the proportion of samples for
        which a given probe must have expression levels above background noise.
        Default: 0.5

    Returns
    -------
    filtered : pandas.DataFrame
        Dataframe containing information on probes that should be retained
        according to intensity-based filtering
    """

    threshold = np.clip(threshold, 0.0, 1.0)

    LGR.info(f'Filtering probes with intensity-based threshold of {threshold}')

    probes = io.read_probes(probes)
    signal, n_samp = np.zeros(len(probes), dtype=int), 0
    for donor, pa in pacall.items():
        annot = io.read_annotation(annotation[donor]).index
        data = io.read_pacall(pa).loc[probes.index, annot]
        n_samp += data.shape[-1]
        # sum binary expression indicator across samples for current subject
        signal += np.asarray(data.sum(axis=1))

    # calculate proportion of signal to noise for given probe across samples
    keep = (signal / n_samp) >= threshold

    LGR.info(f'{keep.sum()} probes survive intensity-based filtering')

    return probes[keep]


def _groupby_structure_id(microarray, annotation):
    """
    Averages samples in `microarray` having identical structure IDs

    Parameters
    ----------
    microarray : (P, S) pandas.DataFrame
        Dataframe should have `P` rows representing probes and `S` columns
        representing distinct samples, with values indicating microarray
        expression levels
    annotation : (S, A) pandas.DataFrame
        Annotation dataframe, obtained by loading a SampleAnnot.csv file from
        Allen Brain Institute

    Returns
    -------
    expression : (P, R) pandas.DataFrame
        Input `microarray` dataframe but with `S` samples averaged into `R`
        regions
    """

    sid = io.read_annotation(annotation)['structure_id']

    return io.read_microarray(microarray).groupby(sid, axis=1).mean()


def _groupby_and_apply(expression, probes, info, applyfunc):
    """
    Subsets `expression` based on most representative probe

    Parameters
    ----------
    expression : dict of (P, S) pandas.DataFrame
        Dictionary where keys are donor IDs and values are dataframes with `P`
        rows representing probes and `S` columns representing distinct samples
    probes : pandas.DataFrame
        Dataframe containing information on probes that should be considered in
        representative analysis. Generally, intensity-based-filtering (i.e.,
        `filter_probes()`) should have been used to reduce this list to only
        those probes with good expression signal
    info : pandas.DataFrame
        Dataframe containing information on probe expression information. Index
        should be unique probe IDs and must have at least 'gene_symbol' column
    applyfunc : callable
        Function used to select representative probe ID from those indexing
        the same gene. Must accept a pandas dataframe as input and return a
        string (i.e., the chosen probe ID)

    Returns
    -------
    representative : dict of (S, G) pandas.DataFrame
        Dictionary where keys are donor IDs and values are dataframes with `S`
        rows representing distinct samples and `G` columns representing unique
        genes
    """

    # group probes by gene and get probe corresponding to relevant feature
    retained = info.groupby('gene_symbol').apply(applyfunc).dropna().squeeze()
    probes = probes.loc[sorted(retained)].sort_values('gene_symbol')

    # subset expression dataframes to retain only desired probes
    representative = {
        d: e.loc[probes.index].T
        for d, e in utils.check_dict(expression).items()
    }

    return representative


def _diff_stability(expression, probes, annotation, *args, **kwargs):
    """
    Picks one probe to represent `expression` data for each gene in `probes`

    If there are multiple probes with expression data for the same gene, this
    function will calculate the similarity of each probes' expression across
    donors and select the probe with the most consistent pattern of regional
    variation (i.e., "differential stability" or DS). Regions are defined by
    the "structure_id" column in `annotation`; similarity is calculated by the
    Spearman correlation coefficient

    Parameters
    ----------
    expression : dict of (P, S) pandas.DataFrame
        Dictionary where keys are donor IDs and values are dataframes with `P`
        rows representing probes and `S` columns representing distinct samples
    probes : pandas.DataFrame
        Dataframe containing information on probes that should be considered in
        representative analysis. Generally, intensity-based-filtering (i.e.,
        `filter_probes()`) should have been used to reduce this list to only
        those probes with good expression signal
    annotation : list of str
        List of filepaths to annotation files from Allen Brain Institute (i.e.,
        as obtained by calling :func:`abagen.fetch_microarray` and accessing
        the `annotation` attribute on the resulting object).

    Returns
    -------
    representative : dict of (S, G) pandas.DataFrame
        Dictionary where keys are donor IDs and values are dataframes with `S`
        rows representing distinct samples and `G` columns representing unique
        genes
    """

    # confirm inputs are expected dictionaries
    expression = utils.check_dict(expression)
    annotation = utils.check_dict(annotation)

    # collapse (i.e., average) expression across AHBA anatomical regions
    region_exp = [
        _groupby_structure_id(microarray, annotation[donor])
        for donor, microarray in expression.items()
    ]

    # get correlation of probe expression across samples for all donor pairs
    probe_exp = np.zeros((len(probes), sum(range(len(expression)))))
    for n, (exp1, exp2) in enumerate(itertools.combinations(region_exp, 2)):
        # samples that current donor pair have in common
        samples = np.intersect1d(exp1.columns, exp2.columns)

        # the ranking process can take a few seconds on each loop
        # unfortunately, we have to do it each time because `samples` changes
        # based on which anatomical regions the two subjects have in common
        exp1 = exp1.loc[:, samples].T.rank()
        exp2 = exp2.loc[:, samples].T.rank()

        probe_exp[:, n] = utils.efficient_corr(exp1, exp2)

    info = pd.DataFrame(dict(gene_symbol=np.asarray(probes.gene_symbol),
                             diff_stability=probe_exp.mean(axis=1)),
                        index=probes.index)
    applyfunc = functools.partial(_max_idx, column='diff_stability')

    return _groupby_and_apply(expression, probes, info, applyfunc)


def _rnaseq(expression, probes, annotation, *args, **kwargs):
    """
    Picks one probe to represent `expression` data for each gene in `probes`

    If there are multiple probes with expression data for the same gene, this
    function will calculate the similarity between each probes' microarray
    expression and RNAseq expression data of the relevant gene, selecting the
    probe with the greatest similarity to the RNAseq data. Regions are defined
    by the "structure_id" column in `annotation`; similarity is calculated by
    the Spearman correlation coefficient.

    Parameters
    ----------
    expression : dict of (P, S) pandas.DataFrame
        Dictionary where keys are donor IDs and values are dataframes with `P`
        rows representing probes and `S` columns representing distinct samples
    probes : pandas.DataFrame
        Dataframe containing information on probes that should be considered in
        representative analysis. Generally, intensity-based-filtering (i.e.,
        `filter_probes()`) should have been used to reduce this list to only
        those probes with good expression signal
    annotation : list of str
        List of filepaths to annotation files from Allen Brain Institute (i.e.,
        as obtained by calling :func:`abagen.fetch_microarray` and accessing
        the `annotation` attribute on the resulting object).

    Returns
    -------
    representative : dict of (S, G) pandas.DataFrame
        Dictionary where keys are donor IDs and values are dataframes with `S`
        rows representing distinct samples and `G` columns representing unique
        genes
    """
    # confirm inputs are expected dictionaries
    expression = utils.check_dict(expression)
    annotation = utils.check_dict(annotation)

    # fetch RNAseq data
    rnaseq = datasets.fetch_rnaseq(donors=expression.keys())

    probe_exp = np.ones((len(probes), len(rnaseq))) * np.nan
    for n, (donor, data) in enumerate(rnaseq.items()):
        # collapse (i.e., average) data across AHAB  anatomical regions
        micro = _groupby_structure_id(expression[donor], annotation[donor])
        rna = _groupby_structure_id(io.read_tpm(data['tpm']),
                                    data['annotation'])

        # get rid of "constant" RNAseq genes
        rna = rna[np.logical_not(np.isclose(rna.std(axis=1, ddof=1), 0))]

        # get matching genes + strcutres between microarray + RNAseq
        regions = np.intersect1d(micro.columns, rna.columns)
        mask = np.isin(np.asarray(probes.gene_symbol),
                       np.intersect1d(probes.gene_symbol, rna.index))
        genes = np.asarray(probes.loc[mask, 'gene_symbol'])
        micro, rna = micro.loc[mask, regions].T, rna.loc[genes, regions].T

        # correlate expression values across regions for each gene
        probe_exp[mask, n] = utils.efficient_corr(micro.rank(), rna.rank())

    mask = np.sum(np.isnan(probe_exp), axis=1) < len(rnaseq)
    info = pd.DataFrame(dict(gene_symbol=np.asarray(probes[mask].gene_symbol),
                             rna_corr=np.nanmean(probe_exp[mask], axis=1)),
                        index=probes.index[mask])
    applyfunc = functools.partial(_max_idx, column='rna_corr')

    return _groupby_and_apply(expression, probes, info, applyfunc)


def _max_idx(df, column):
    """
    Returns probe ID with max index in `df`

    Parameters
    ----------
    df : (P, 1) pandas.DataFrame
        Dataframe with `P` rows indicating distinct probes and one column
        containing summary statistic of probe expression
    column : str, optional
        Column name from which to extract the max index. If not specified uses
        the first numerical column.

    Returns
    -------
    probe_id : str
        ID of probe selected as representative for given gene
    """

    return df.idxmax(numeric_only=True)[column]


def _max_loading(df):
    """
    Returns probe ID with max loading along first principal component of `df`

    Parameters
    ----------
    df : (P, S) pandas.DataFrame
        Dataframe with `P` rows indicating distinct probe expression values
        across `S` samples

    Returns
    -------
    probe_id : str
        ID of probe selected as representative for given gene
    """

    if len(df) == 1:
        return df.index[0]

    data = np.asarray(df)
    data = data - data.mean(axis=0, keepdims=True)

    # svd() is faster than eig() here because we don't need to construct the
    # covariance matrix
    u, s, v = np.linalg.svd(data, full_matrices=False)

    # use sign flip based on right singular vectors (as we would with eig())
    v *= np.sign(v[range(len(v)), np.argmax(np.abs(v), axis=1)])[:, np.newaxis]

    return df.index[(data @ v.T)[:, 0].argmax()]


def _correlate(df, method):
    """
    Returns probe ID with max avg correlation (>2 probes) or `method` (<2)

    Parameters
    ----------
    df : (P, S) pandas.DataFrame
        Dataframe with `P` rows indicating distinct probe expression values
        across `S` samples
    method : {'variance', 'intensity'}
        Method for selecting representative probe when only two probes index
        the same gene (>2 probes uses correlation method)

    Returns
    -------
    probe_id : str
        ID of probe selected as representative for given gene
    """

    data = np.asarray(df)

    if len(data) > 2:
        xmax = np.mean((np.corrcoef(data) + 1) / 2, axis=1)
    elif len(data) == 2:
        if method == 'variance':
            xmax = np.std(data, axis=1)
        elif method == 'intensity':
            xmax = np.mean(data, axis=1)
    else:
        return df.index[0]

    return df.index[xmax.argmax()]


def _average(expression, probes, *args, **kwargs):
    """
    Averages expression data for probes representing the same gene

    Parameters
    ----------
    expression : dict of (P, S) pandas.DataFrame
        Dictionary where keys are donor IDs and values are dataframes with `P`
        rows representing probes and `S` columns representing distinct samples
    probes : pandas.DataFrame
        Dataframe containing information on microarray probes that should be
        considered in representative analysis. Generally intensity-based
        filtering (i.e., via :func:`filter_probes()`) should have been used to
        reduce this list to only those probes with good expression signal

    Returns
    -------
    representative : dict of (S, G) pandas.DataFrame
        Dictionary where keys are donor IDs and values are dataframes with `S`
        rows representing distinct samples and `G` columns representing unique
        genes
    """

    def _avg(df):
        return df.rename(probes['gene_symbol'].to_dict()) \
                 .rename_axis('gene_symbol') \
                 .groupby('gene_symbol') \
                 .mean().T

    return {d: _avg(exp) for d, exp in utils.check_dict(expression).items()}


def _collapse(expression, probes, *args, method='max_variance', **kwargs):
    """
    Selects one representative probe per gene using provided `method`

    Parameters
    ----------
    expression : dict of (P, S) pandas.DataFrame
        Dictionary where keys are donor IDs and values are dataframes with `P`
        rows representing probes and `S` columns representing distinct samples
    probes : pandas.DataFrame
        Dataframe containing information on microarray probes that should be
        considered in representative analysis. Generally intensity-based
        filtering (i.e., via :func:`filter_probes()`) should have been used to
        reduce this list to only those probes with good expression signal
    method : str, optional
        Method by which to select represenative probes for each gene. Must be
        one of ['max_variance', 'max_intensity', 'pc_loading', 'corr_variance',
        'corr_intensity']

    Returns
    -------
    representative : dict of (S, G) pandas.DataFrame
        Dictionary where keys are donor IDs and values are dataframes with `S`
        rows representing distinct samples and `G` columns representing unique
        genes
    """

    # concatenate all donors into giant probe x sample expression dataframe
    expression = utils.check_dict(expression)
    probe_exp = pd.concat(expression.values(), axis=1)

    # determine aggregation function based on provided method; also reduce
    # probe expression if required (i.e., max_variance, max_intensity)
    if method == 'max_variance':
        probe_exp = pd.DataFrame(probe_exp.std(axis=1), columns=[method])
        agg = functools.partial(_max_idx, column=method)
    elif method == 'max_intensity':
        probe_exp = pd.DataFrame(probe_exp.mean(axis=1), columns=[method])
        # probe_exp.name = method
        agg = functools.partial(_max_idx, column=method)
    elif method == 'pc_loading':
        agg = _max_loading
    elif method in ['corr_variance', 'corr_intensity']:
        agg = functools.partial(_correlate, method=method[5:])
    else:
        raise ValueError(f'Provided method {method} is invalid. Please check '
                         'inputs and try again.')

    info = pd.merge(probes[['gene_symbol']], probe_exp, on='probe_id')
    return _groupby_and_apply(expression, probes, info, agg)


_max_variance = functools.partial(_collapse, method='max_variance')
_max_intensity = functools.partial(_collapse, method='max_intensity')
_pc_loading = functools.partial(_collapse, method='pc_loading')
_corr_variance = functools.partial(_collapse, method='corr_variance')
_corr_intensity = functools.partial(_collapse, method='corr_intensity')


SELECTION_METHODS = dict(
    mean=_average,
    average=_average,
    max_intensity=_max_intensity,
    max_variance=_max_variance,
    pc_loading=_pc_loading,
    corr_variance=_corr_variance,
    corr_intensity=_corr_intensity,
    diff_stability=_diff_stability,
    rnaseq=_rnaseq
)
AGG_METHODS = [  # can only be used with `donor_probes='aggregate'`
    'mean', 'average', 'diff_stability', 'rnaseq'
]
COLLAPSE_METHODS = [  # methods that don't SELECT but COLLAPSE ACROSS probes
    'mean', 'average'
]


def collapse_probes(microarray, annotation, probes, method='diff_stability',
                    donor_probes='aggregate'):
    """
    Reduces `microarray` to a sample x gene expression dataframe

    Using provided `method`, reduces `microarray` expression data by either (1)
    selecting a representative probe amongst all probes indexing the same gene,
    or (2) collapsing across all probes indexing the same gene. See Notes for
    more information on different methods available.

    Parameters
    ----------
    microarray : dict of str or pandas.DataFrame
        Dictionary where keys are donor IDs and values are filepaths to (or
        dataframes of) MicroarrayExpression.csv files from Allen Brain
        Institute
    annotation : dict of str or pandas.DataFrame
        Dictionary where keys are donor IDs and values are filepaths to (or
        dataframes of) SampleAnnot.csv files from Allen Brain Institute. Only
        used if `method='diff_stability'`
    probes : str or pandas.DataFrame
        Filepath to Probes.csv or dataframe containing information on
        microarray probes that should be considered in representative analysis.
        Generally intensity-based filtering (i.e., via :func:`filter_probes()`)
        should have been used to reduce this list to only those probes with
        good expression signal
    method : str, optional
        Selection method for subsetting (or collapsing across) probes from the
        same gene. Must be one of 'average', 'max_intensity', 'max_variance',
        'pc_loading', 'corr_intensity', 'corr_variance', 'diff_stability', or
        'rnaseq'; see Notes for more information. Default: 'diff_stability'
    donor_probes : str, optional
        Whether specified `probe_selection` method should be performed with
        microarray data from all donors ('aggregate'), independently for each
        donor ('independent'), or based on the most common selected probe
        across donors ('common'). Not all combinations of `probe_selection`
        and `donor_probes` methods are viable. Default: 'aggregate'

    Returns
    -------
    expression : dict of (S, G) pandas.DataFrame
        Dictionary where keys are donor IDs and values are dataframes with `S`
        rows representing distinct samples and `G` columns representing unique
        genes. Entries of dataframe indicate microarray expression levels for
        each combination of sample + gene. Columns will be identical across all
        dataframes, but `S` will vary by donor.

    Notes
    -----
    The following methods can be used for collapsing across probes when
    multiple probes are available for the same gene.

    1. ``method='average'``

    Uses the average of expression data across probes indexing the same gene as
    in [PR3]_, [PR5]_, [PR9]_, [PR10]_, [PR15]_, [PR16]_, and [PR17]_. Using
    `method='mean'` will do the same thing.

    2. ``method='max_intensity'``

    Selects probe with maximum average expression across samples from all
    donors as in [PR14]_.

    3. ``method='max_variance'``

    Selects probe with maximum variance in expression across samples from all
    donors as in [PR12]_.

    4. ``method='pc_loading'``

    Selects probe with maximum loading on first principal component of
    decomposition performed across samples from all donors as in [PR13]_.

    5. ``method='corr_intensity'``

    Selects probe with maximum correlation to other probes from same gene when
    >2 probes exist; otherwise, uses same procedure as `method=max_intensity`.
    Used in [PR1]_ and [PR11]_.

    6. ``method='corr_variance'``

    Selects probe with maximum correlation to other probes from same gene when
    >2 probes exist; otherwise, uses same procedure as `method=max_variance`.
    Used in [PR2]_, [PR4]_, and [PR6]_.

    7. ``method='diff_stability'``

    Selects probe with the most consistent pattern of regional variation across
    donors (i.e., highest average correlation across brain regions between all
    pairs of donors) as in [PR7]_ and [PR8]_.

    8. ``method='rnaseq'``

    Selects probes with most consistent pattern of regional variation to RNAseq
    data (across the two donors with RNAseq data).

    References
    ----------
    .. [PR1] Anderson, K. M., Krienen, F. M., Choi, E. Y., Reinen, J. M., Yeo,
       B. T., & Holmes, A. J. (2018). Gene expression links functional networks
       across cortex and striatum. Nature Communications, 9(1), 1428.

    .. [PR2] Burt, J. B., Demirtaş, M., Eckner, W. J., Navejar, N. M., Ji, J.
       L., Martin, W. J., ... & Murray, J. D. (2018). Hierarchy of
       transcriptomic specialization across human cortex captured by structural
       neuroimaging topography. Nature Neuroscience, 21(9), 1251.

    .. [PR3] Eising, E., Huisman, S. M., Mahfouz, A., Vijfhuizen, L. S.,
       Anttila, V., Winsvold, B. S., ... & Boomsma, D. I. (2016). Gene
       co-expression analysis identifies brain regions and cell types involved
       in migraine pathophysiology: a GWAS-based study using the Allen Human
       Brain Atlas. Human Genetics, 135(4), 425-439.

    .. [PR4] Forest, M., Iturria‐Medina, Y., Goldman, J. S., Kleinman, C. L.,
       Lovato, A., Oros Klein, K., ... & Greenwood, C. M. (2017). Gene networks
       show associations with seed region connectivity. Human Brain Mapping,
       38(6), 3126-3140.

    .. [PR5] French, L., & Paus, T. (2015). A FreeSurfer view of the cortical
       transcriptome generated from the Allen Human Brain Atlas. Frontiers in
       Neuroscience, 9, 323.

    .. [PR6] Hawrylycz, M. J., Lein, E. S., Guillozet-Bongaarts, A. L., Shen,
       E. H., Ng, L., Miller, J. A., ... & Abajian, C. (2012). An anatomically
       comprehensive atlas of the adult human brain transcriptome. Nature,
       489(7416), 391.

    .. [PR7] Hawrylycz, M., Miller, J. A., Menon, V., Feng, D., Dolbeare, T.,
       Guillozet-Bongaarts, A. L., ... & Glasser, M. F. (2015). Canonical
       genetic signatures of the adult human brain. Nature Neuroscience,
       18(12), 1832.

    .. [PR8] Kirsch, L., & Chechik, G. (2016). On expression patterns and
       developmental origin of human brain regions. PLoS Computational Biology,
       12(8), e1005064.

    .. [PR9] Krienen, F. M., Yeo, B. T., Ge, T., Buckner, R. L., & Sherwood, C.
       C. (2016). Transcriptional profiles of supragranular-enriched genes
       associate with corticocortical network architecture in the human brain.
       Proceedings of the National Academy of Sciences, 113(4), E469-E478.

    .. [PR10] McColgan, P., Gregory, S., Seunarine, K. K., Razi, A., Papoutsi,
       M., Johnson, E., ... & Scahill, R. I. (2018). Brain regions showing
       white matter loss in Huntington’s disease are enriched for synaptic and
       metabolic genes. Biological Psychiatry, 83(5), 456-465.

    .. [PR11] Myers, E. M., Bartlett, C. W., Machiraju, R., & Bohland, J. W.
       (2015). An integrative analysis of regional gene expression profiles in
       the human brain. Methods, 73, 54-70.S

    .. [PR12] Negi, S. K., & Guda, C. (2017). Global gene expression profiling
       of healthy human brain and its application in studying neurological
       disorders. Scientific Reports, 7(1), 897.

    .. [PR13] Parkes, L., Fulcher, B. D., Yücel, M., & Fornito, A. (2017).
       Transcriptional signatures of connectomic subregions of the human
       striatum. Genes, Brain and Behavior, 16(7), 647-663.

    .. [PR14] Romero-Garcia, R., Whitaker, K. J., Váša, F., Seidlitz, J.,
       Shinn, M., Fonagy, P., ... & Vértes, P. E. (2018). Structural covariance
       networks are coupled to expression of genes enriched in supragranular
       layers of the human cortex. NeuroImage, 171, 256-267.

    .. [PR15] Tan, P. P. C., French, L., & Pavlidis, P. (2013). Neuron-enriched
       gene expression patterns are regionally anti-correlated with
       oligodendrocyte-enriched patterns in the adult mouse and human brain.
       Frontiers in Neuroscience, 7, 5.

    .. [PR16] Vértes, P. E., Rittman, T., Whitaker, K. J., Romero-Garcia, R.,
       Váša, F., Kitzbichler, M. G., ... & Goodyer, I. M. (2016). Gene
       transcription profiles associated with inter-modular hubs and connection
       distance in human functional magnetic resonance imaging networks.
       Philosophical Transactions of the Royal Society B: Biological Sciences,
       371(1705), 20150362.

    .. [PR17] Whitaker, K. J., Vértes, P. E., Romero-Garcia, R., Váša, F.,
       Moutoussis, M., Prabhu, G., ... & Tait, R. (2016). Adolescence is
       associated with genomically patterned consolidation of the hubs of the
       human brain connectome. Proceedings of the National Academy of Sciences,
       113(32), 9105-9110.
    """

    try:
        collfunc = SELECTION_METHODS[method]
    except KeyError:
        raise ValueError(f'Provided `method` "{method}" is invalid; must be '
                         f'one of {list(SELECTION_METHODS)}')

    valid_probes = ['aggregate', 'common', 'independent']
    if donor_probes not in valid_probes:
        raise ValueError(f'Provided `donor_probes` "{donor_probes}" is '
                         f'invalid; must be one of {valid_probes}')

    LGR.info(f'Reducing probes indexing same gene with method: {method}')
    # subset microarray data for pre-selected probes + samples
    # this will also left/right mirror samples, if previously requested
    probes = io.read_probes(probes)
    microarray = utils.check_dict(microarray)
    annotation = utils.check_dict(annotation)
    for donor, micro in microarray.items():
        samp = io.read_annotation(annotation[donor]).index
        microarray[donor] = io.read_microarray(micro).loc[probes.index, samp]

    # now, "collect" the probes based on the provided `method`
    if method in AGG_METHODS or donor_probes == 'aggregate':
        # perform the collection function for all donors, together
        microarray = collfunc(microarray, probes, annotation)
    elif donor_probes == 'independent':
        # perform the collection function for each donor separately
        for donor in microarray:
            microarray.update(collfunc({donor: microarray[donor]}, probes,
                                       {donor: annotation[donor]}))
    elif donor_probes == 'common':
        # perform collection function for each donor separately and retain ONLY
        # the chose probe IDs
        probe_ids = [
            collfunc(
                {donor: microarray[donor]}, probes, {donor: annotation[donor]}
            )[donor].columns
            for donor in microarray
        ]
        # find the mode of the probe IDs chosen across donors and then subset
        # those probes from the original microarray dataframes for all donors
        probe_ids = np.squeeze(sstats.mode(probe_ids, axis=0)[0])
        for donor in microarray:
            microarray[donor] = microarray[donor].loc[probe_ids].T

    # convert probe IDs as column names to gene symbols
    if method not in COLLAPSE_METHODS:
        for donor, micro in microarray.items():
            symbols = probes.loc[micro.columns, 'gene_symbol']
            micro = micro.set_axis(symbols, axis=1, copy=True)
            microarray[donor] = micro.sort_index(axis=1)

    n_genes = utils.first_entry(microarray).shape[-1]
    LGR.info(f'{n_genes} genes remain after probe filtering + selection')

    return microarray
