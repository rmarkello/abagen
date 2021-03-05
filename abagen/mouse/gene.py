# -*- coding: utf-8 -*-
"""
Functions to make mouse gene queries and manipulations
"""

import re
import pandas as pd
from .utils import _coerce_inputs, _make_api_query


# available attributes of gene query
_GENE_ATTRIBUTES = [
    'acronym',
    'alias_tags',
    'chromosome_id',
    'ensembl_id',
    'entrez_id',
    'genomic_reference_update_id',
    'homologene_id',
    'id',
    'legacy_ensembl_gene_id',
    'name',
    'organism_id',
    'original_name',
    'original_symbol',
    'reference_genome_id',
    'sphinx_id',
    'version_status'
]


def available_gene_info():
    """ Lists available attributes for :func:`abagen.mouse.get_gene_info`
    """

    return _GENE_ATTRIBUTES.copy()


def get_gene_info(id=None, acronym=None, name=None, attributes=None,
                  verbose=False):
    """
    Queries Allen API for information about given gene

    One of `id`, `acronym`, or `name` must be provided.

    Parameters
    ----------
    id : int, optional
        Numerical gene ID
    acronym : str, optional
        Short-form gene acronym (case sensitive)
    name : str, optional
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
    info : pandas.DataFrame
        If `attributes` is a str, returns an int or str depending on specified
        attribute. If `attributes` is a list, return a dict where keys are
        attributes and values are str or int.

    Raises
    ------
    ValueError
        The provided gene is invalid

    Examples
    --------
    Get gene ID and name corresponding to gene acronym 'Pdyn':

    >>> from abagen import mouse
    >>> mouse.get_gene_info(acronym='Pdyn',
    ...                     attributes=['id', 'name'])  # doctest: +NORMALIZE_WHITESPACE
                id          name
    acronym
    Pdyn     18376  prodynorphin

    You can also supply multiple genes to the query:

    >>> mouse.get_gene_info(acronym=['Ace', 'Cd99'],
    ...                     attributes=['id', 'name'])  # doctest: +NORMALIZE_WHITESPACE
                 id                                               name
    acronym
    Ace       11210  angiotensin I converting enzyme (peptidyl-dipe...
    Cd99     163028                                       CD99 antigen
    """  # noqa

    criteria = [
        _coerce_inputs(id=id, acronym=acronym, name=name),
        'products[id$eq1]'
    ]
    provided = re.search(r'\[(\S+)\$', criteria[0]).group(1)

    # determine which attributes to request; if we don't have to request all
    # of them then we can speed up the API call
    if attributes is None:
        attributes = _GENE_ATTRIBUTES
    elif isinstance(attributes, str):
        attributes = [attributes]
    attributes = [a for a in attributes if a not in provided]

    for attr in attributes:
        if attr not in _GENE_ATTRIBUTES:
            raise ValueError('Provided attribute "{}" is invalid; please '
                             'check valid attributes with '
                             'abagen.mouse.available_gene_info().'
                             .format(attr))

    info = _make_api_query('Gene', criteria=criteria,
                           attributes=attributes + [provided], verbose=verbose)
    info = pd.DataFrame(info).set_index(provided)[attributes]

    return info.sort_index()
