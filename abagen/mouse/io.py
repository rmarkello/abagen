# -*- coding: utf-8 -*-
"""
Functions to fetch (and load) mouse gene and structure lists from Allen API
"""

import json
import os.path as op
from pkg_resources import resource_filename

import pandas as pd
import requests


def fetch_allenref_genes(entry_type=None, cache=True, verbose=False):
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
        Whether to print status message. Default: False

    Returns
    -------
    genes : list or :obj:`pandas.DataFrame`
        Genes in Allen Reference database

    Notes
    -----
    May require internet access to make query to the Allen API (which will take
    some time); after query is made once the results are cached.
    """

    entries = ['id', 'acronym', 'name']
    url = 'http://api.brain-map.org/api/v2/data/Gene/query.json' \
          '?criteria=products[id$in1,12]&only=id,acronym,name' \
          '&start_row={}&num_rows=1000'
    fname = resource_filename('abagen', 'data/allen_reference_genes.csv')

    # check that provided entry_type is valid
    if entry_type is not None and entry_type not in entries:
        raise ValueError('Provided entry_type {} is not valid. Specified '
                         'entry_type must be one of {}.'
                         .format(entry_type, entries))

    # if file doesn't exist or we want to overwrite the cache for some reason
    # download the data from the Allen API
    if not op.isfile(fname) or not cache:
        genes = pd.DataFrame(columns=entries)
        if verbose:
            print('Gene information not available locally; querying '
                  'Allen API for information. This may take some time...')
        row = 0
        while True:
            root = json.loads(requests.get(url.format(row)).content)
            if not root['success']:
                raise ValueError('Unable to query Allen API at this time; '
                                 'please try again later.')
            genes = genes.append(pd.DataFrame(root['msg'])[entries])
            row += root['num_rows']
            if row == root['total_rows']:
                break
        genes = genes.sort_values('id').reset_index(drop=True)
        # save information to disk
        genes.to_csv(fname, index=False)
    else:
        genes = pd.read_csv(fname)

    # extract only relevant entry_type, if desired
    if entry_type is not None:
        genes = genes[entry_type].tolist()

    return genes


def fetch_allenref_structures(entry_type=None, cache=True, verbose=False):
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
        Whether to print status message. Default: False

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
    url = 'http://api.brain-map.org/api/v2/data/Structure/query.json' \
          '?criteria=ontology[id$in1,12]&num_rows=all&only=id,acronym,name'
    fname = resource_filename('abagen', 'data/allen_reference_atlas.csv')

    if entry_type is not None and entry_type not in entries:
        raise ValueError('Provided entry_type {} is not valid. Specified '
                         'entry_type must be one of {}.'
                         .format(entry_type, entries))

    if not op.isfile(fname) or not cache:
        if verbose:
            print('Structure information not available locally; querying '
                  'Allen API for information...')
        root = json.loads(requests.get(url).content)
        if not root['success']:
            raise ValueError('Unable to query Allen API at this time; please '
                             'try again later.')
        structures = pd.DataFrame(root['msg'])
        # sort entries by structure ID
        structures = structures.sort_values('id').reset_index(drop=True)
        # save information to disk
        structures.to_csv(fname, index=False)
    else:
        structures = pd.read_csv(fname)

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
