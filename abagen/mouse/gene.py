# -*- coding: utf-8 -*-
"""
Functions to make mouse gene expression queries and manipulations
"""
import requests
from xml.etree import ElementTree as ET

# the attributes of gene query
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


def _check_gene_validity(gene_id=None, gene_acronym=None, gene_name=None,
                         verbose=False):
    """
    Check if a provided gene is valid and has records in the Allen database

    Parameters
    ----------
    gene_id : int, optional
        Numerical gene ID
    gene_acronym : str, optional
        Short-form gene acronym (case sensitive)
    gene_name : str, optional
        Full gene name (case sensitive)
    verbose : bool, optional
        Whether to print query URL. Default: False

    Returns
    -------
    validity : boolean
        Whether the provided gene has records in the Allen database
    root : :obj:`xml.etree.ElementTree`
        XML response from Allen API; empty if `validity` is False

    Raises
    ------
    TypeError
        If missing parameters

    Example
    -------
    >>> # check if gene ID 18376 is valid
    >>> validity, _ = check_structure_validity(gene_id=18376)
    >>> validity
    True
    >>> # check if structure Pdyn is valid
    >>> validity, root = check_structure_validity(gene_acronym='Pdyn')
    """

    url = 'http://api.brain-map.org/api/v2/data/SectionDataSet/query.xml' \
          '?criteria=products[id$eq1],genes[{}$eq{}]' \
          '&include=genes,plane_of_section'

    # if gene ID is given
    # preferred: id > acronym > name
    if gene_id is not None:
        query_url = url.format('id', gene_id)
    elif gene_acronym is not None:
        query_url = url.format('acronym', gene_acronym)
    elif gene_name is not None:
        query_url = url.format('name', gene_name)
    else:
        raise TypeError('At least one of [\'gene_id\', \'gene_acronym\','
                        '\'gene_name\'] must be specified.')

    # make the query
    if verbose:
        print("Querying {}...".format(query_url))
    root = ET.fromstring(requests.get(query_url).content)

    if root.attrib['total_rows'] != '0':  # successful
        return True, root
    else:
        return False, root


def available_gene_info():
    """ Lists available attributes for :py:func:`abagen.mouse.get_gene_info`
    """

    return _GENE_ATTRIBUTES


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
        :py:func:`abagen.mouse.available_gene_info` for list of available
        attributes to request. If not specified all available attributes will
        be returned. Default: None
    verbose : bool, optional
        Whether to print status messages. Default: False

    Returns
    -------
    gene_info : int, str, or dict
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
    Get gene name according to gene name 'Pdyn':
    >>> get_gene_info(gene_acronym='Pdyn', attributes='name')
    'prodynorphin'

    Get gene acronym according to gene id 1will hve t8376:
    >>> gene_acronym = get_gene_info(gene_id=18376, attributes='acronym')
    'Pdyn'
    """

    validity, root = _check_gene_validity(gene_id=gene_id,
                                          gene_acronym=gene_acronym,
                                          gene_name=gene_name)
    if validity is False:
        provided = [gene_id, gene_acronym, gene_name]
        raise ValueError('Gene {} is invalid. Try another gene.'
                         .format(*[i for i in provided if i is not None]))

    # if the query was successful
    attr_list = _GENE_ATTRIBUTES if attributes is None else attributes

    if isinstance(attr_list, list):
        gene_info = dict()
        for attr in attr_list:
            try:
                # extract the info of attr
                gene_info[attr] = _get_single_gene_attribute(root, attr)
            except AttributeError:
                print('Invalid attribute {}; skipping...'.format(attr))
                continue
    else:  # single attribute is given
        # return the single attr value, or raise AttributeError
        return _get_single_gene_attribute(root, attributes)

    return gene_info


def _get_single_gene_attribute(root, attr):
    """
    Finds attribute `attr` in `root` XML tree

    Parameters
    ----------
    root : :obj:`xml.etree.ElementTree`
        XML response for Allen API
    attr : str
        Attribute to query from XML response `root`

    Returns
    -------
    value: int, str, or None
        Value of `attr`

    Raises
    ------
    AttributeError
        If `attr` does not exist in `root`
    """

    item = root.findall('section-data-sets/section-data-set/genes/gene/{}'
                        .format(attr))

    # check if attr is valid (if any information is found)
    if len(item) == 0:
        raise AttributeError('There is no gene attribute called {}'
                             .format(attr))

    # check data type
    try:
        return int(item[0].text)
    except ValueError:
        return item[0].text
    except TypeError:
        # the attribute exists, but has no value
        return None
