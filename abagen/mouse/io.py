# -*- coding: utf-8 -*-
"""
Functions to fetch (and load) mouse gene and structure lists from Allen API
"""

import os.path as op
from pkg_resources import resource_filename

import pandas as pd

from .utils import _make_api_query


def fetch_allenref_genes(entry_type=None, cache=True, verbose=True):
    """
    Loads all genes from Allen Reference database

    Parameters
    ----------
    entry_type: {'id', 'acronym', 'name'}, optional
        The type of gene identifier to load. Specifying 'id' returns a list of
        numerical gene IDs, 'acronym' returns a list of short-form gene
        acronyms, and 'name' returns full gene names. If not specified, returns
        a dataframe with all information. Default: None
    cache : bool, optional
        Whether to use cached gene information (if it exists). Setting to False
        will overwrite cache. Default: True
    verbose : bool, optional
        Whether to print status message. Default: True

    Returns
    -------
    genes : list or :obj:`pandas.DataFrame`
        Genes in Allen Reference database

    Notes
    -----
    May require internet access to make query to the Allen API (which will take
    some time); after query is made once the results are cached.
    """

    # check that provided entry_type is valid
    entries = ['id', 'acronym', 'name']
    if entry_type is not None and entry_type not in entries:
        raise ValueError('Provided entry_type {} is not valid. Specified '
                         'entry_type must be one of {}.'
                         .format(entry_type, entries))

    # if file doesn't exist or we want to overwrite the cache for some reason
    # download the data from the Allen API
    fname = resource_filename('abagen', 'data/allen_reference_genes.csv')
    if not op.isfile(fname) or not cache:
        if verbose:
            print('Gene information not available locally; querying '
                  'Allen API for information. This may take some time...')
        out = _make_api_query('Gene', criteria='products[id$eq1]',
                              attributes=entries, suffix='num_rows=all')
        genes = pd.DataFrame(out)
        # sort entries by gene ID
        genes = genes.sort_values('id').reset_index(drop=True)
        # save information to disk
        genes.to_csv(fname, index=False)
    else:
        genes = pd.read_csv(fname)[entries]

    # extract only relevant entry_type, if desired
    if entry_type is not None:
        genes = genes[entry_type].tolist()

    return genes


def fetch_allenref_structures(entry_type=None, cache=True, verbose=True):
    """
    Loads all anatomical structures in the Allen Reference Atlas

    Parameters
    ----------
    entry_type: {'id', 'acronym', 'name'}, optional
        The type of structural identifier to load. Specifying 'id' returns a
        list of numerical structure IDs, 'acronym' returns a list of short-form
        structure acronyms, and 'name' returns full structure names. If not
        specified, returns a dataframe with all information. Default: None
    cache : bool, optional
        Whether to use cached structure information (if it exists). Setting to
        False will overwrite cache. Default: True
    verbose : bool, optional
        Whether to print status message. Default: True

    Returns
    -------
    structures : list or :obj:`pandas.DataFrame`
        Anatomical structures in Allen Reference Atlas

    Notes
    -----
    May require internet access to make query to the Allen API (which will take
    some time); after query is made once the results are cached.
    """

    entries = ['id', 'acronym', 'name']
    if entry_type is not None and entry_type not in entries:
        raise ValueError('Provided entry_type {} is not valid. Specified '
                         'entry_type must be one of {}.'
                         .format(entry_type, entries))

    fname = resource_filename('abagen', 'data/allen_reference_atlas.csv')
    if not op.isfile(fname) or not cache:
        if verbose:
            print('Structure information not available locally; querying '
                  'Allen API for information...')
        out = _make_api_query('Structure', criteria='ontology[id$eq1]',
                              attributes=entries, suffix='num_rows=all')
        structures = pd.DataFrame(out)
        # sort entries by structure ID
        structures = structures.sort_values('id').reset_index(drop=True)
        # save information to disk
        structures.to_csv(fname, index=False)
    else:
        structures = pd.read_csv(fname)[entries]

    # extract only relevant entry_type, if desired
    if entry_type is not None:
        structures = structures[entry_type].tolist()

    return structures


def fetch_rubinov2015_structures(entry_type=None):
    """
    Loads subset of anatomical structures in Allen Reference Atlas used in [1]_

    Parameters
    ----------
    entry_type: {'id', 'acronym', 'name'}, optional
        The type of structural identifier to load. Specifying 'id' returns a
        list of numerical structure IDs, 'acronym' returns a list of short-form
        structure acronyms, and 'name' returns full structure names. If not
        specified, returns a dataframe with all information. Default: None

    Returns
    -------
    structures : list or :obj:`pandas.DataFrame`
        Anatomical structures in Allen Reference Atlas from [1]_

    References
    ----------
    .. [1] Rubinov, M., Ypma, R. J., Watson, C., & Bullmore, E. T. (2015).
       Wiring cost and topological participation of the mouse brain connectome.
       Proceedings of the National Academy of Sciences, 112(32), 10032-10037.
    """

    entries = ['id', 'acronym', 'name']
    if entry_type is not None and entry_type not in entries:
        raise ValueError('Provided entry_type {} is not valid. Specified '
                         'entry_type must be one of {}.'
                         .format(entry_type, entries))

    fname = resource_filename('abagen', 'data/rubinov2015_pnas.csv')
    structures = pd.read_csv(fname)[entries]

    if entry_type is not None:
        structures = structures[entry_type].tolist()

    return structures
