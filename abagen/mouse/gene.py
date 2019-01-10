# -*- coding: utf-8 -*-
"""
functions to make mouse gene expression queries and manipulations
"""
import requests
from xml.etree import ElementTree as ET

# specify RMA query model, restrain the results to mouse
URL_PREFIX = "http://api.brain-map.org/api/v2/data/" \
             "SectionDataSet/query.xml?" \
             "criteria=products[id$eq1],"
# information to include
URL_INCLUDE = "&include=genes"

# an alternative to make the query
#URL_PREFIX = "http://api.brain-map.org/api/v2/data/" \
#             "query.xml?include=model::Gene"
# restrain the queries to mouse (products ID=1)
#URL_INCLUDE = ",products[id$eq1]"

# the attributes of gene query
GENE_ATTRIBUTES = [
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


def check_gene_validity(gene, entry_type='id'):
    """
    check if a gene is valid.
    :param gene: str or int
        gene identifier, either an int (id, entrez-id, chromosome-id,
        etc.), or a str (acronym, name, alias-tags, etc.), specified
        in entry_type.
        if acronym is given, it should be capitalized.
        recommend to use gene ID (int) or gene acronym (str).
    :param entry_type: str, optional
        the type of gene identifier. default: 'id'
        supported:
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

    :return: boolean
        whether the gene is valid
    """
    if entry_type not in GENE_ATTRIBUTES:
        raise ValueError('entry_type {} is invalid'
                         .format(entry_type))

    # URL
    if isinstance(gene, int):
        query_url = URL_PREFIX + \
            "genes[{0}$eq{1}]".format(entry_type, gene) + \
            URL_INCLUDE
    elif isinstance(gene, str):
        query_url = URL_PREFIX + \
            "genes[{0}$eq'{1}']".format(entry_type, gene) + \
            URL_INCLUDE
    else:
        raise ValueError('{} is an invalid gene identifier'
                         .format(gene))

    print('accessing {}...'.format(query_url))

    # make the query
    r = requests.get(query_url)
    root = ET.fromstring(r.content)

    return (True, root) if root.text else False


def get_gene_info(gene, entry_type='id', attributes='all'):
    """
    get gene attributes according to the identifier
    (e.g., gene id or acronym)

    :param gene: int or str
        gene identifier, either an ID (id, entrez-id, chromosome-id,
        etc.), or a name (acronym, name, alias-tags, etc.), specified
        in entry_type.
        if acronym is given, it should be capitalized.
    :param entry_type: str, optional
        the type of gene identifier. default: 'id'
        supported:
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

    :param attributes: str, or list, optional
        the attributes of the gene to return.
        default:'all', returning the value of all gene attributes.
        supported:
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

    :return: int, str or dict
        can be integer (id, entrez-id, chromosome-id, etc.), or
        str (name, acronym, etc.), if only one attribute is specified
        or dict {attribute:value} if multiple attributes are given
    """
    if check_gene_validity(gene, entry_type) is False:
        raise ValueError('{0} {1} is invalid'
                         .format(entry_type, gene))

    else:
        _, root = check_gene_validity(gene, entry_type)

    # if the query was successful
    if attributes == 'all':
        attr_list = GENE_ATTRIBUTES
    elif isinstance(attributes, list):
        attr_list = attributes

    gene_info = dict()
    if attributes == 'all' or isinstance(attributes, list):
        for attr in attr_list:
            try:
                # extract the info of attr
                gene_info[attr] = _get_single_gene_attribute(
                    root, attr
                )
            except AttributeError:
                print('There is no attribute called {}. '
                      'Skipped.'.format(attr))
                continue

    else:  # single attribute is given
        # return the single attr value, or AttributeError
        return _get_single_gene_attribute(root, attributes)

    return gene_info


def _get_single_gene_attribute(root, attr):
    """
    return the value of a single gene attribute
    :param root: object ElementTree 'Response'
    :param attr: str,
        attribute to return
    :return: str or int,
        the value of the attribute
    """
    item = root.find(
        'ncbi-genes/ncbi-gene/{}'.format(attr)
    )
    # check if attr is valid (if any information is found)
    if not item:
        raise AttributeError(
            'There is no attribute called {}'.format(attr)
        )

    # check data type
    # attr is an integer
    if 'type' in item.attrib \
            and item.attrib['type'] == 'integer':
        return int(item.text)
    else:
        # attr is a string
        return item.text
