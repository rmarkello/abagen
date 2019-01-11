# -*- coding: utf-8 -*-
"""
functions to make mouse gene expression queries and manipulations
"""
import requests
from xml.etree import ElementTree as ET

# specify RMA query model, restrain the results to mouse
# and to the genes that have section data sets available
URL_PREFIX = "http://api.brain-map.org/api/v2/data/" \
             "SectionDataSet/query.xml?" \
             "criteria=products[id$eq1],"
# information to include
URL_INCLUDE = "&include=genes,plane_of_section"

# #an alternative to make the query, but has no type info
# URL_PREFIX = "http://api.brain-map.org/api/v2/data/" \
#             "query.xml?include=model::Gene"
# #restrain the queries to mouse (products ID=1)
# URL_INCLUDE = ",products[id$eq1]"

GENE_ENTRY_TYPES = [
    'acronym',
    'chromosome-id',
    'ensembl-id',
    'entrez-id',
    'homologene-id',
    'id',
    'legacy-ensembl-gene-id',
    'name',
    'original-name',
    'original-symbol',
    'sphinx-id',
]

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


def check_gene_validity(gene_id=None, acronym=None, name=None):
    """
    check if a structure is valid or has records in the database
    Example
    -------
    # check if gene ID 18376 is valid
    validity, root = check_structure_validity(gene_id=18376)
    # check if structure Pdyn is valid
    validity, root = check_structure_validity(acronym='Pdyn')

    Parameters
    ----------
    gene_id: int, optional
        gene ID
    acronym: str, optional
        gene acronym (case sensitive)
    name: str, optional
        gene name (case sensitive)

    Returns
    -------
    validity: boolean
        if the gene has records in the database
    root: obj,
        Element 'Response' object, empty if query fails

    Raises
    ------
    TypeError:
        missing parameters
    """
    # if gene ID is given
    # preferred: id > acronym > name
    if gene_id is not None:
        query_url = URL_PREFIX + \
            "genes[id$eq{}]".format(gene_id) + \
            URL_INCLUDE
    elif acronym is not None:
        query_url = URL_PREFIX + \
            "genes[acronym$eq'{}']".format(acronym) + \
            URL_INCLUDE
    elif name is not None:
        query_url = URL_PREFIX + \
            "genes[name$eq'{}']".format(name) + \
            URL_INCLUDE
    else:
        raise TypeError(
            "at least one gene identifier should be specified"
        )
    # make the query
    print("access {}...".format(query_url))
    r = requests.get(query_url)
    root = ET.fromstring(r.content)

    if root.attrib['total_rows'] != '0':  # successful
        return True, root
    else:
        return False, root


def get_gene_info(gene_id=None, acronym=None, name=None, attributes='all'):
    """
    get attributes of a gene

    Examples
    --------
    # get gene name according to gene name 'Pdyn'
    gene_name = get_gene_info(acronym='Pdyn', attributes='name')
    # get gene acronym according to gene id 18376
    gene_acronym = get_gene_info(gene_id=18376, attributes='acronym')

    # multiple attributes
    gene_info = get_gene_info(
        gene_id=18376, attributes=['acronym', 'name']
    )
    # gene name
    gene_info['name']

    Parameters
    ----------
    gene_id: int, optional
        gene ID
    acronym: str, optional
        gene acronym (case sensitive)
    name: str, optional
        gene name (case sensitive)
    attributes: str or list, optional
        a single attribute or a list of attributes
        default: 'all', returning all the attributes
        available attributes:
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

    Returns
    -------
    gene_info: int, str or dict
        if a single attribute is given, return an int or str
        if multiple attributes are given, return a dict
        {attr:value}. attr is str (attribute) and value is str or int

    Raises
    ------
    ValueError:
        the gene given is invalid
    AttributeError:
        only one attribute is given, and it is invalid
    """
    validity, root = check_gene_validity(
        gene_id=gene_id, acronym=acronym, name=name
    )
    if validity is False:
        raise ValueError('The gene is invalid. '
                         'Try another gene.')

    # if the query was successful
    if attributes == 'all':
        attr_list = GENE_ATTRIBUTES
    else:
        attr_list = attributes

    if isinstance(attr_list, list):
        gene_info = dict()
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
        # return the single attr value, or raise AttributeError
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
    item = root.findall(
        'section-data-sets/section-data-set/'
        'genes/gene/{}'.format(attr)
    )
    # check if attr is valid (if any information is found)
    if len(item) == 0:
        raise AttributeError(
            'There is no gene attribute called {}'.format(attr)
        )

    # check data type
    # attr is an integer
    try:
        return int(item[0].text)
    except ValueError:
        return item[0].text
    except TypeError:
        # the attribute exists, but has no value
        return None


