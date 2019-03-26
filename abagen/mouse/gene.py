# -*- coding: utf-8 -*-
"""
Functions to make mouse gene queries and manipulations
"""

import json
import requests

# available attributes of gene query
_GENE_ATTRIBUTES = [
    'acronym',
    'alias-tags',
    'chromosome-id',
    'ensembl-id',
    'entrez-id',
    'genomic-reference-update-id',
    'homologene-id',
    'id',
    'legacy-ensembl-gene-id',
    'name',
    'organism-id',
    'original-name',
    'original-symbol',
    'reference-genome-id',
    'sphinx-id',
    'version-status'
]


def _make_gene_query(gene_id=None, gene_acronym=None, gene_name=None,
                     attributes=None, suffix=None, verbose=False):
    """
    Check if a provided gene is valid and has records in the Allen database

    Requires internet access to make query to the Allen API

    Parameters
    ----------
    gene_id : int, optional
        Numerical gene ID
    gene_acronym : str, optional
        Short-form gene acronym (case sensitive)
    gene_name : str, optional
        Full gene name (case sensitive)
    attributes : list, optional
        List of attributes for which to query gene information
    verbose : bool, optional
        Whether to print query URL. Default: False

    Returns
    -------
    response : list
        Response from Allen API; empty list if no response

    Raises
    ------
    TypeError
        When none of the required parameters are specified

    Example
    -------
    Check if gene ID '18376' is valid:
    >>> out = _check_gene_validity(gene_id=18376)
    >>> len(out)
    1

    Check if gene acronym 'Pdyn' is valid:
    >>> out = _check_gene_validity(gene_acronym='Pdyn')
    >>> len(out)
    1
    """

    url = 'http://api.brain-map.org/api/v2/data/Gene/query.json' \
          '?criteria=[{}${}{}],products[id$eq1]'

    if attributes is not None:
        url += '&only={}'.format(','.join(attributes))
    if suffix is not None:
        url += suffix

    # preferred order of inputs: id > acronym > name
    query, params = 'eq', [gene_id, gene_acronym, gene_name]
    for criteria, param in zip(['id', 'acronym', 'name'], params):
        if param is not None:
            if isinstance(param, list):
                query, param = 'in', ','.join([str(f) for f in param])
            break
    if param is None:
        params = ['gene_id', 'gene_acronym', 'gene_name']
        raise TypeError('At least one of {} must be specified.'.format(params))

    # formt URL with desired information and make query to API
    query_url = url.format(criteria, query, param)
    if verbose:
        print("Querying {}...".format(query_url))
    response = json.loads(requests.get(query_url).content)

    return response.get('msg', [])


def available_gene_info():
    """ Lists available attributes for :func:`abagen.mouse.get_gene_info`
    """

    return _GENE_ATTRIBUTES.copy()


def get_gene_info(gene_id=None, gene_acronym=None, gene_name=None,
                  attributes=None, verbose=False):
    """
    Queries Allen API for information about given gene

    One of `gene_id`, `gene_acronym`, or `gene_name` must be provided.

    Parameters
    ----------
    gene_id : int, optional
        Numerical gene ID
    gene_acronym : str, optional
        Short-form gene acronym (case sensitive)
    gene_name : str, optional
        Full gene name (case sensitive)
    attributes : str or list, optional
        Which attributes / information to obtain for the provided gene. See
        :func:`abagen.mouse.available_gene_info` for list of available
        attributes to request. If not specified all available attributes will
        be returned. Default: None
    verbose : bool, optional
        Whether to print status messages. Default: False

    Returns
    -------
    info : int, str, or dict
        If `attributes` is a str, returns an int or str depending on specified
        attribute. If `attributes` is a list, return a dict where keys are
        attributes and values are str or int.

    Raises
    ------
    ValueError
        The provided gene is invalid
    AttributeError
        Only one attribute is specified but it is invalid

    Examples
    --------
    Get gene name corresponding to gene acronym 'Pdyn':
    >>> get_gene_info(gene_acronym='Pdyn', attributes='name')
    'prodynorphin'

    Get gene acronym corresponding to gene id 18376:
    >>> get_gene_info(gene_id=18376, attributes='acronym')
    'Pdyn'
    """

    # determine which attributes to request; if we don't have to request all
    # of them then we can speed up the API call
    if attributes is None:
        attributes = _GENE_ATTRIBUTES
    elif isinstance(attributes, str):
        attributes = [attributes]
    attributes = [attr.replace('-', '_') for attr in attributes]

    # query API and check if we received a valid output
    info = _make_gene_query(gene_id, gene_acronym, gene_name,
                            attributes=attributes, verbose=verbose)
    if len(info) == 0:
        provided = [gene_id, gene_acronym, gene_name]
        raise ValueError('Gene {} is invalid.'
                         .format(*[i for i in provided if i is not None]))

    if len(attributes) == 1:
        info = [f[attributes[0]] for f in info]
        if len(info) == 1:
            return info[0]

    return info
