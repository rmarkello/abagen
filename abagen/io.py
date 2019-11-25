# -*- coding: utf-8 -*-
"""
Functions for loading the various files associated with the AHBA microarray
dataset.

This also contains functionality for optionally converting the downloaded CSV
files to parquet format, which provides much faster I/O access / quicker load
times.
"""

import os.path as op
import pandas as pd
try:
    eng = pd.io.parquet.get_engine('fastparquet')
    assert 'SNAPPY' in eng.api.compression.compressions
    use_parq = True
# pandas version too low OR don't have fastparquet installed
except (AttributeError, ImportError, AssertionError):
    use_parq = False


def _make_parquet(fname, convert_only=False):
    """
    Loads `fname`, converting to parquet file if it does not already exist

    Parameters
    ----------
    fname : str
        Path to data file
    convert_only : bool, optional
        Check if parquet version of `fname` exists and convert if it doesn't;
        if it does, just return. Default: False

    Returns
    -------
    data : pandas.DataFrame
        Data loaded from `fname`
    """

    # get ideal parquet filename
    parqname = fname.rpartition('.csv')[0] + '.parq'

    # if it exists, load it as parquet
    if op.exists(parqname):
        if convert_only:
            return
        data = pd.read_parquet(parqname, engine='fastparquet')
    # otherwise, load CSV and save to parquet
    else:
        data = pd.read_csv(fname, header=None)
        data.columns = data.columns.astype(str)
        data.to_parquet(parqname, engine='fastparquet')
        if convert_only:
            return

    return data


def read_microarray(fname, copy=False, parquet=True):
    """
    Loads MicroarrayExpression.csv file found at `fname`

    Microarray files contain raw expression data for all the tissue samples
    taken from a single donor across all genetic probes.

    Parameters
    ----------
    fname : str
        Path to MicroarrayExpression.csv file
    copy : bool, optional
        Whether to return a copy if `fname` is a pre-loaded pandas.Dataframe.
        Default: False
    parquet : bool, optional
        Whether to load data from parquet file instead of CSV. If a parquet
        file does not already exist then one will be created for faster loading
        in the future. Only available if ``fastparquet`` and ``python-snappy``
        module are installed. Default: True

    Returns
    -------
    microarray : (P, S) pandas.DataFrame
        Dataframe containing microarray expression data, where `P` is probes
        and `S` is samples. The row index is the unique probe ID assigned
        during processing, which can be used to match data to the information
        obtained with :func:`read_probes`. The column index is the unique
        sample ID (integer, beginning at 0) which can be used to match data to
        the information obtained with :func:`read_annotation`.
    """

    try:
        if use_parq and parquet:
            data = _make_parquet(fname, convert_only=False)
            data = data.set_index('0')
        else:
            data = pd.read_csv(fname, header=None, index_col=0)
        data.index.name = 'probe_id'
        data.columns = pd.Series(range(len(data.columns)), name='sample_id')
    except (AttributeError, ValueError):
        if not isinstance(fname, pd.DataFrame):
            raise TypeError('Provided fname must be filepath to Microarray'
                            'Expression.csv file from Allen Human Brain '
                            'Atlas.')
        data = fname.copy() if copy else fname

    return data


def read_ontology(fname, copy=False):
    """
    Loads Ontology.csv file found at `fname`

    Ontology files contain information on the anatomical delineations used by
    the Allen Institute when obtaining samples from donor brains, and are used
    in their online Brain Viewer to colorize regions. These files should be the
    same for every donors.

    This information can be used to ensure that microarray samples are
    appropriately matched to anatomical regions.

    Parameters
    ----------
    fname : str
        Path to Ontology.csv file
    copy : bool, optional
        Whether to return a copy if `fname` is a pre-loaded pandas.Dataframe.
        Default: False

    Returns
    -------
    ontology : (R, 8) pandas.DataFrame
        Dataframe containing ontology information for `R` anatomical regions
        used by the Allen Institute. Columns include: 'id', 'acronym', 'name',
        'parent_structure_id', 'hemisphere', 'graph_order',
        'structure_id_path', and 'color_hex_triplet'.
    """

    try:
        data = pd.read_csv(fname)
    except ValueError:
        if not isinstance(fname, pd.DataFrame):
            raise TypeError('Provided fname must be filepath to Ontology.csv '
                            'file from Allen Human Brain Atlas.')
        data = fname.copy() if copy else fname

    return data


def read_pacall(fname, copy=False, parquet=True):
    """
    Loads PACall.csv file found at `fname`

    PA files contain a present/absent flag indicating whether the corresponding
    probe's expression is above background noise. It is set to 1 when both of
    the following conditions are met:

        1. The mean signal of the probe's expression is significantly different
           from the corresponding background, as assessed by a 2-sided t-test
           where p < 0.01, and
        2. The difference between the background subtracted signal and the
           background is significant (> 2.6 * background standard deviation).

    This information can be used to discard "noisy" probes that might not be
    contributing high-quality expression information.

    Parameters
    ----------
    fname : str
        Path to PACall.csv file
    copy : bool, optional
        Whether to return a copy if `fname` is a pre-loaded pandas.Dataframe.
        Default: False
    parquet : bool, optional
        Whether to load data from parquet file instead of CSV. If a parquet
        file does not already exist then one will be created for faster loading
        in the future. Only available if ``fastparquet`` and ``python-snappy``
        module are installed. Default: True

    Returns
    -------
    pacall : (P, S) pandas.DataFrame
        Dataframe containing a binary indicator determining whether expression
        information for each probe exceeded background noise in a given sample,
        where `P` is probes and `S` is samples. The row index is the unique
        probe ID assigned during processing, which can be used to match data to
        the information obtained with :func:`read_probes`. The column index is
        the unique sample ID (integer, beginning at 1) which can be used to
        match data to the information obtained with :func:`read_annotation`.
    """

    try:
        if use_parq and parquet:
            data = _make_parquet(fname, convert_only=False)
            data = data.set_index('0')
        else:
            data = pd.read_csv(fname, header=None, index_col=0)
        data.index.name = 'probe_id'
        data.columns = pd.Series(range(len(data.columns)), name='sample_id')
    except (AttributeError, ValueError):
        if not isinstance(fname, pd.DataFrame):
            raise TypeError('Provided fname must be filepath to PACall.csv '
                            'file from Allen Human Brain Atlas.')
        data = fname.copy() if copy else fname

    return data


def read_probes(fname, copy=False):
    """
    Loads Probes.csv file found at `fname`

    Probe files contain metadata on all genetic probes used in the AHBA data.
    These files should be the same for every donor.

    This information can be used to e.g., query expression data for certain
    genes, collapse data across probes from the same gene, etc.

    Parameters
    ----------
    fname : str
        Path to Probes.csv file
    copy : bool, optional
        Whether to return a copy if `fname` is a pre-loaded pandas.Dataframe.
        Default: False

    Returns
    -------
    probes : (P, 6) pandas.DataFrame
        Dataframe containing information for `P` genetic probes. The row index
        is the unique probe ID assigned during processing, which can be used to
        match metadata to information obtained with :func:`read_microarray` and
        :func:`read_pacall`. Columns include 'probe_name', 'gene_id',
        'gene_symbol', 'gene_name', 'entrez_id', and 'chromosome'.
    """

    try:
        data = pd.read_csv(fname, index_col=0)
    except ValueError:
        if not isinstance(fname, pd.DataFrame):
            raise TypeError('Provided fname must be filepath to Probes.csv '
                            'file from Allen Human Brain Atlas.')
        data = fname.copy() if copy else fname

    return data


def read_annotation(fname, copy=False):
    """
    Loads SampleAnnot.csv file found at `fname`

    Sample annotation files contain metadata on all the tissue samples taken
    from a single donor brain, including the spatial location of the samples.

    This information can be used to combine samples within the same anatomical
    region across donors.

    Parameters
    ----------
    fname : str
        Path to SampleAnnot.csv file
    copy : bool, optional
        Whether to return a copy if `fname` is a pre-loaded pandas.Dataframe.
        Default: False

    Returns
    -------
    annotation : (S, 13) pandas.DataFrame
        Dataframe containing structural information on `S` samples. The row
        index is the unique sample ID (integer, beginning with 1) which can be
        used to match data to the information obtained with
        :func:`read_microarray` and :func:`read_pacall`. Columns include
        'structure_id', 'slab_num', 'well_id', 'slab_type',
        'structure_acronym', 'structure_name', 'polygon_id', 'mri_voxel_x',
        'mri_voxel_y', 'mri_voxel_z', 'mni_x', 'mni_y', 'mni_z'.
    """

    try:
        data = pd.read_csv(fname)
        data.index.name = 'sample_id'
    except ValueError:
        if not isinstance(fname, pd.DataFrame):
            raise TypeError('Provided fname must be filepath to Annotation'
                            '.csv file from Allen Human Brain Atlas.')
        data = fname.copy() if copy else fname

    return data
