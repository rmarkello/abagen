# -*- coding: utf-8 -*-
"""
Functions for cleaning and processing the AHBA microarray dataset
"""
import itertools
from nibabel.volumeutils import Recoder
import numpy as np
import pandas as pd
from abagen import io, utils

# AHBA structure IDs corresponding to different brain parts
ONTOLOGY = Recoder(
    (('4008', 'cerebral cortex', 'cortex'),
     ('4275', 'cerebral nuclei', 'subcortex'),
     ('4391', 'diencephalon', 'subcortex'),
     ('9001', 'mesencephalon', 'subcortex'),
     ('4696', 'cerebellum', 'cerebellum'),
     ('4833', 'metencephalon', 'cerebellum'),
     ('9512', 'myelencephalon', 'cerebellum')),
    fields=('id', 'name', 'structure')
)


def filter_probes(pacall, probes, threshold=0.5):
    """
    Performs intensity based filtering (IBF) of expression probes

    Uses binary indicator for expression levels in `pacall` to determine which
    probes have expression levels above background noise in `threshold` of
    samples across donors.

    Parameters
    ----------
    pacall : list
        List of PACall files from Allen Brain Institute. Optimally obtained by
        calling `abagen.fetch_microarray()` and accessing the `pacall`
        attribute on the resulting object
    probes : list
        List of probes files from Allen Brain Institute. Optimally obtained by
        calling `abagen.fetch_microarray()` and accessing the `probes`
        attribute on the resulting object
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

    signal, samples = [], 0
    for fname in pacall:
        data = io.read_pacall(fname)
        samples += data.shape[-1]
        # sum binary expression indicator across samples for current subject
        signal.append(data.sum(axis=1).values)

    # calculate proportion of signal to noise for given probe across samples
    keep = (np.sum(signal, axis=0) / samples) > threshold

    # read in probe file (they're all the same) and drop "bad" probes
    filtered = io.read_probes(probes[0])[keep]

    return filtered


def get_stable_probes(microarray, annotation, probes):
    """
    Picks one probe to represent `microarray` data for each gene in `probes`

    If there are multiple probes with expression data for the same gene, this
    function will calculate the similarity of each probes' expression across
    donors and select the probe with the most consistent pattern of regional
    variation (i.e., "differential stability" or DS). Regions are
    operationalized by the "structure_id" column in `annotation`.

    Parameters
    ----------
    microarray : list of str
        List of microarray expression files from Allen Brain Institute.
        Optimally obtained by calling `abagen.fetch_microarray()` and accessing
        the `microarray` attribute on the resulting object.
    annotation : list of str
        List of annotation files from Allen Brain Institute. Optimally obtained
        by calling `abagen.fetch_microarray()` and accessing the `annotation`
        attribute on the resulting object.
    probes : pandas.DataFrame
        Dataframe containing information on probes that should be considered in
        representative analysis. Generally, intensity-based-filtering (i.e.,
        `probe_ibf()`) should have been used to reduce this list to only those
        probes with good expression signal

    Returns
    -------
    representative : pandas.DataFrame
        Dataframe containing information on probes that are most representative
        of their genes based on differential stability analysis

    References
    ----------
    .. [1] Hawrylycz, M., Miller, J. A., Menon, V., Feng, D., Dolbeare, T.,
       Guillozet-Bongaarts, A. L., ... & Lein, E. (2015). Canonical genetic
       signatures of the adult human brain. Nature Neuroscience, 18(12), 1832.
    """

    # read in microarray data for all subjects
    num_subj = len(microarray)

    # this is a relatively slow procedure (i.e., takes a couple of seconds)
    micro = [_reduce_micro(microarray[n], annotation[n], probes)
             for n in range(num_subj)]

    # get correlation of probe expression across samples for all donor pairs
    probe_corrs = np.zeros((len(probes), sum(range(num_subj))))
    for n, (s1, s2) in enumerate(itertools.combinations(range(num_subj), 2)):

        # find samples that current donor pair have in common
        samples = np.intersect1d(micro[s1].columns, micro[s2].columns)

        # the ranking process can take a few seconds on each loop
        # unfortunately, we have to do it each time because `samples` changes
        probe_corrs[:, n] = utils.efficient_corr(micro[s1][samples].T.rank(),
                                                 micro[s2][samples].T.rank())

    # group probes by gene and get probe corresponding to max correlation
    df = pd.DataFrame(dict(gene_symbol=probes.gene_symbol.values,
                           corrs=probe_corrs.mean(axis=1),
                           probe_id=probes.index))
    retained = (df.groupby('gene_symbol')
                  .apply(lambda x: x.loc[x.corrs.idxmax(), 'probe_id']))

    return probes[probes.index.isin(retained.values)]


def _reduce_micro(micro, annot, probes):
    """
    Collapse expression data in `micro` taken from same anatomical structure

    Parameters
    ----------
    micro : str
        Filepath to microarray file for single donor
    annot : str
        Filepath to annotation file for single donor; should be same donor as
        `micro`
    probes : pandas.DataFrame
        Probes that should be considered for analysis (i.e., those that made it
        through intensity-based filtering)

    Returns
    -------
    micro : (P x S) pd.core.frame.DataFrame
        Microarray expression data, where `P` is probes and `S` is anatomical
        structures
    """
    # read in microarray data and set columns as anatomical structure ID
    micro = io.read_microarray(micro)
    micro.columns = io.read_annotation(annot).structure_id

    # average across same structures
    micro = micro.groupby(micro.columns, axis=1).mean()

    # only keep good probes
    return micro[micro.index.isin(probes.index)]


def drop_mismatch_samples(annotation, ontology):
    """
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

    Returns
    -------
    keep : pandas.DataFrame
        Annotation data
    """

    # read in data files
    annot = io.read_annotation(annotation)
    ont = io.read_ontology(ontology)

    # add hemisphere + brain "structure" designation to annotation data
    hemi = dict(zip(ont.id, ont.hemisphere))
    stru = dict(zip(ont.id, ont.structure_id_path.apply(_get_path_structure)))
    annot = annot.assign(hemisphere=annot.structure_id.replace(hemi),
                         structure=annot.structure_id.replace(stru))

    # only keep samples with consistent hemisphere + MNI coordinate designation
    keep = annot.query('(hemisphere == "L" & mni_x < 0) | '
                       '(hemisphere == "R" & mni_x > 0)')

    return keep


def _get_path_structure(structure_path):
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
    ids = list(ONTOLOGY.value_set('id').intersection(structure_path))

    try:
        return ONTOLOGY.structure[ids[0]]
    except IndexError:
        return


def normalize_expression(expression):
    """
    Performs scaled robust sigmoid (SRS) normalization on `expression` data

    Parameters
    ----------
    expression : (S, G) pandas.DataFrame
        Microarray expression data, where `S` is samples (or regions) and `G`
        is genes

    Returns
    -------
    normalized : (S, G) pandas.DataFrame
        Data from `expression` normalized separately for each gene
    """

    # calculate sigmoid normalization
    data = expression.get_values()
    norm = (data - np.median(data, axis=0)) / np.std(data, axis=0)
    srs = 1 / (1 + np.exp(-norm))

    # rescale normalized values to a unit interval
    srs_min = srs.min(axis=0, keepdims=True)
    srs_max = srs.max(axis=0, keepdims=True)
    scaled = (srs - srs_min) / (srs_max - srs_min)

    # recreate dataframe
    normalized = pd.DataFrame(scaled,
                              columns=expression.columns,
                              index=expression.index)

    return normalized
