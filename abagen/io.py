# -*- coding: utf-8 -*-
"""
Functions for loading the various files associated with the AHBA microarray
dataset

This also contains functionality for optionally converting the downloaded CSV
files to the parquet format, which provides much faster I/O access and quicker
loading times
"""

import os.path as op
import pandas as pd
try:
    pd.io.parquet.get_engine('fastparquet')
    use_parq = True
except ImportError:
    use_parq = False


def _make_parquet(fname, convert_only=False):
    """
    Loads `fname`, converting to parquet file if it does not already exist

    Parameters
    ----------
    fname : str
        Path to large data file
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

    # do some cleaning up of the data
    data = data.set_index('0')
    data.index.name = 'probe_id'
    data.columns = pd.Series(range(len(data.columns)), name='sample_id')

    return data


def read_microarray(fname, parquet=True):
    """
    Reads in Allen Brain Institute microarray file found at `fname`

    Parameters
    ----------
    fname : str
        Path to microarray expression file
    parquet : bool, optional
        Whether, if available, to load data from parquet file instead of CSV.
        Data will be saved as a parquet file for faster loading in the future
        if such a file does not already exist

    Returns
    -------
    microarray : (P, S) pandas.DataFrame
        Dataframe containing microarray expression data, where `P` is probes
        and `S` is samples. Therow index is the unique probe ID assigned during
        processing
    """
    if not isinstance(fname, str):
        if isinstance(fname, pd.DataFrame):
            return fname.copy()
        else:
            raise TypeError('Provided fname {} must be a filepath.'
                            .format(fname))

    if use_parq and parquet:
        return _make_parquet(fname)

    return pd.read_csv(fname, header=None, index_col=0)


def read_ontology(fname, parquet=True):
    """
    Reads in Allen Brain Institute ontology file found at `fname`

    Parameters
    ----------
    fname : str
        Path to ontology file
    parquet : bool, optional
        Does nothing; here simply for compatibility with other io functionality

    Returns
    -------
    ontology : (R, 8) pandas.DataFrame
        Dataframe containing ontology information for all labelled brain
        structures used during sample collection
    """
    if not isinstance(fname, str):
        if isinstance(fname, pd.DataFrame):
            return fname.copy()
        else:
            raise TypeError('Provided fname {} must be a filepath.'
                            .format(fname))

    return pd.read_csv(fname)


def read_pacall(fname, parquet=True):
    """
    Reads in Allen Brain Institute PA call file found at `fname`

    Parameters
    ----------
    fname : str
        Path to PA call file
    parquet : bool, optional
        Whether, if available, to load data from parquet file instead of CSV.
        Data will be saved as a parquet file for faster loading in the future
        if such a file does not already exist

    Returns
    -------
    pacall : (P, S) pandas.DataFrame
        Dataframe containing a binary indicator determining whether expression
        information for each probe exceeded background noise in a given sample,
        where `P` is probes and `S` is samples
    """
    if not isinstance(fname, str):
        if isinstance(fname, pd.DataFrame):
            return fname.copy()
        else:
            raise TypeError('Provided fname {} must be a filepath.'
                            .format(fname))

    if use_parq and parquet:
        return _make_parquet(fname)

    return pd.read_csv(fname, header=None, index_col=0)


def read_probes(fname, parquet=True):
    """
    Reads in Allen Brain Institute probes file found at `fname`

    Parameters
    ----------
    fname : str
        Path to probes file
    parquet : bool, optional
        Does nothing; here simply for compatibility with other io functionality

    Returns
    -------
    probes : (P, 6) pandas.DataFrame
        Dataframe containing genetic information for `P` probes
    """
    if not isinstance(fname, str):
        if isinstance(fname, pd.DataFrame):
            return fname.copy()
        else:
            raise TypeError('Provided fname {} must be a filepath.'
                            .format(fname))

    return pd.read_csv(fname, index_col=0)


def read_annotation(fname, parquet=True):
    """
    Reads in Allen Brain Institute annotation file found at `fname`

    Parameters
    ----------
    fname : str
        Path to annotation file
    parquet : bool, optional
        Does nothing; here simply for compatibility with other io functionality

    Returns
    -------
    annotation : (S, 13) pandas.DataFrame
        Dataframe containing structural information on `S` samples
    """
    if not isinstance(fname, str):
        if isinstance(fname, pd.DataFrame):
            return fname.copy()
        else:
            raise TypeError('Provided fname {} must be a filepath.'
                            .format(fname))

    annotation = pd.read_csv(fname)
    annotation.index.name = 'sample_id'

    return annotation
