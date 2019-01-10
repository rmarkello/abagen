# -*- coding: utf-8 -*-

import requests
from xml.etree import ElementTree as ET
import numpy as np
import pandas as pd
from .gene import (check_gene_validity, 
                   get_gene_info, 
                   GENE_ENTRY_TYPES, 
                   GENE_ATTRIBUTES)
from .io import read_all_structures

API_QUERY_STRING = 'http://api.brain-map.org/api/v2/data/' \
                   'SectionDataSet/query.xml?'


UNIONIZATION_ATTRIBUTES = [
    'expression-density',
    'expression-energy',
    'id',
    'section-data-set-id',
    'structure-id',
    'sum-expressing-pixel-intensity',
    'sum-expressing-pixels',
    'sum-pixel-intensity',
    'sum-pixels',
    'voxel-energy-cv',
    'voxel-energy-mean'
]


def get_unionization_from_experiment(
        experiment_id,
        structure_list=None,
        attributes='all'
):
    """
    fetches mouse Unionization data of a single experiment,
    either saggital or coronal, according to the experiment id

    Parameters
    ----------
    experiment_id: int, specifying the experiment id
    structure_list: a list of strings, optional
        the list of structures in the form of acronyms.
        default: structures as documented in Rubinov et al, 2015
    attributes: list, optional
        specify the unionization data attributes to include
        default: 'all'
        available attributes:
            'expression-density',
            'expression-energy', (gene expression)
            'id',
            'section-data-set-id',
            'structure-id',
            'sum-expressing-pixel-intensity',
            'sum-expressing-pixels',
            'sum-pixel-intensity',
            'sum-pixels',
            'voxel-energy-cv',
            'voxel-energy-mean'

    Returns
    -------
    unionization: dict
        unionization data attributes and the values {attribute:value}

    Raises
    ------
    ValueError:
        If experiment_id is invalid
    """

    if structure_list is None:
        # read default structure list
        # see Rubinov et al 2015 for more information
        structure_list = read_all_structures()

    # make the query
    api_query_criteria = 'id={}'.format(experiment_id)
    api_query_include = '&include=structure_unionizes%28structure%29'
    query_url = API_QUERY_STRING + \
        api_query_criteria + \
        api_query_include
    r = requests.get(query_url)
    root = ET.fromstring(r.content)

    if not root.text:  # check if any expressions are found
        raise ValueError(
            'No gene expression values are found '
            'associated experiment {}. '
            'Try another valid experiment ID'
            .format(experiment_id)
        )

    unionization = dict()
    path_prefix = 'section-data-sets/section-data-set/' \
                  'structure-unionizes/structure-unionize/'

    if attributes == 'all':  # find all expression-relevant values
        for attr in GENE_ATTRIBUTES:
            # find all items included in structure_list

            unionization[attr] = _get_single_unionization_attribute(
                root, attr, path_prefix, structure_list
            )
            
    # only one attribute is specified
    elif isinstance(attributes, str):
        return _get_single_unionization_attribute(
                root, attributes, path_prefix, structure_list
            )
    else:  # if multiple attributes are specified
        for attr in attributes:
            try:
                unionization[attr] = _get_single_unionization_attribute(
                root, attr, path_prefix, structure_list
            )

            except AttributeError:
                print('There is no attribute called {}. '
                      'Ignored.'.format(attr))
                continue

    return unionization


def _get_single_unionization_attribute(root, attr, path_prefix, structure_list):
    """
    return values of a single unionization attribute
    :param root:
        Element 'Response' to parse
    :param attr: str,
        the attribute to return
    :param path_prefix: str,
        specifying the path to find the attribute
    :param structure_list:
        list of structures to include
    :return: numpy_array
        the values of attr corresponding to the
        structures in structure_list
    """
    roi_count = len(structure_list)
    # if the attribute exists
    if not root.findall(path_prefix + attr):
        raise AttributeError(
            'There is no attribute called {}.'.format(attr)
        )
    all_items = [
        (float(val_item.text), structure_item.text)
        for val_item, structure_item in zip(
            root.findall(path_prefix + attr),
            root.findall(path_prefix + 'structure/acronym')
        ) if structure_item.text in structure_list
    ]

    structures_in_the_list = [item[1] for item in all_items]
    vals_in_the_list = np.array(
        [item[0] for item in all_items],
        dtype=np.float
    )

    vals = np.empty((roi_count,))
    vals[:] = np.nan
    # average duplicate expressions
    for k in range(roi_count):
        index = [
            idx
            for idx, item
            in enumerate(structures_in_the_list)
            if item == structure_list[k]
        ]
        # if gene expressions are found in kth region
        if index:
            vals[k] = vals_in_the_list[index].mean()
        else:
            print('No {0} values '
                  'are found for region {1}. '
                  'Set to NaN.'
                  .format(attr, structure_list[k]))

    return vals


def get_experiment_id_from_gene(gene, slicing_direction='sagittal'):
    """
    fetches mouse gene expression data of a single experiment,
    either saggital or coronal, according to the experiment id

    Parameters
    ----------
    gene: str or int
        specifying the acronym (capitalized) or ID of the gene
    slicing_direction: str, optional
        slicing scheme of the samples, 'sagittal' or 'coronal'.
        default: 'sagittal'

    Returns
    -------
    experiment_id_list: list
        integers (experiment IDs) of the gene specified by gene_acronym
    """
    if slicing_direction not in ['sagittal', 'coronal']:
        raise ValueError('Slicing slicing_direction {} is invalid. '
                         'Try sagittal or coronal instead'
                         .format(slicing_direction))

    if check_gene_validity(gene) is False:
        raise ValueError('Gene {} is invalid'.format(gene))

    # find the experiment IDs associated with this gene
    if isinstance(gene, int):  # if gene id is provided
        api_query_criteria = 'criteria=products%5Bid$eq1%5D,' \
                             'genes%5Bid$eq%27{}%27%5D,' \
                             'plane_of_section[name$eq%27{}%27]'\
            .format(gene, slicing_direction)
    else:  # if gene acronym is provided
        api_query_criteria = 'criteria=products%5Bid$eq1%5D,' \
                             'genes%5Bacronym$eq%27{}%27%5D,' \
                             'plane_of_section[name$eq%27{}%27]' \
            .format(gene, slicing_direction)
    # extra conditions
    api_query_include = '&include=genes,section_images'
    
    # make the query
    query_url = API_QUERY_STRING + \
        api_query_criteria + \
        api_query_include
    r = requests.get(query_url)
    root = ET.fromstring(r.content)

    experiment_id_list = []
    for item in root.findall('section-data-sets/section-data-set/id'):
        experiment_id_list.append(int(item.text))  # append experiment id

    if not experiment_id_list:
        print('No {0} experiments are found for gene {1}, '
              'return an empty experiment ID list'
              .format(slicing_direction, gene))

    return experiment_id_list


def get_unionization_from_gene(
        gene,
        slicing_direction='sagittal',
        structure_list=None,
        attributes='all'
):
    """
    fetches mouse gene expression data
    averaged across all the associated experiments
    (either sagittal or coronal) of a specific gene

    Parameters
    ----------
    gene: string or int
        specifying the acronym (capitalized) or ID of the gene
    slicing_direction: string, optional
        slicing scheme of the samples, 'sagittal' or 'coronal'.
        default: 'sagittal'
    structure_list: a list of strings, optional
        the list of ROIs in the form of acronyms.
        default: ROIs as included in Rubinov et al, 2015
    attributes: list, optional
        specify the unionization data attributes to include
        default: 'all'
        available attributes:
            'expression-density',
            'expression-energy', (gene expression)
            'id',
            'section-data-set-id',
            'structure-id',
            'sum-expressing-pixel-intensity',
            'sum-expressing-pixels',
            'sum-pixel-intensity',
            'sum-pixels',
            'voxel-energy-cv',
            'voxel-energy-mean'

    Returns
    -------
    numpy_array or dict
        if a single attribute is given, return a (N, ) numpy_array
        corresponding to the structures in structure_list
        if a list of attributes is given, return a dict

    """
    if attributes == 'all':
        attr_list = UNIONIZATION_ATTRIBUTES
    elif isinstance(attributes, list):
        attr_list = attributes

    if slicing_direction not in ['sagittal', 'coronal']:
        raise ValueError('Slicing slicing_direction {} is invalid. '
                         'Try sagittal or coronal instead'
                         .format(slicing_direction))

    if structure_list is None:
        # read default ROI list
        # (see Rubinov et al, 2015 for the criteria of
        # choosing the ROIs)
        roilabels = pd.read_csv(
            "abagen/data/roilabels-rubinov2015pnas.csv"
        )
        structure_list = roilabels['roiacronyms'].values
    roi_count = len(structure_list)

    experiment_id_list = get_experiment_id_from_gene(
        gene, slicing_direction=slicing_direction
    )

    # initialize a dict to store attr values and valid experimentIDs
    unionization = dict()

    if attributes == 'all' \
            or isinstance(attributes, list):
        for attr in attr_list:
            attr_vals = np.empty((0, roi_count))
            for experiment_id in experiment_id_list:
                # values of a single attr
                try:
                    vals = get_unionization_from_experiment(
                            experiment_id,
                            structure_list,
                            attributes=attr
                        )
                except AttributeError:
                    print('There is no attribute called {}. Skipped.'
                          .format(attr))
                    continue
                attr_vals = np.append(
                    attr_vals,
                    vals.reshape((1, roi_count)),
                    axis=0
                )
                #valid_experiment_id_list.append(experiment_id)
            unionization[attr] = attr_vals

    else:  # a single attribute is given
        attr_vals = np.empty((0, roi_count))
        for experiment_id in experiment_id_list:
            vals = get_unionization_from_experiment(
                experiment_id,
                structure_list,
                attributes=attributes
            )
            attr_vals = np.append(
                attr_vals,
                vals.reshape((1, roi_count)),
                axis=0
            )
        return attr_vals

    return unionization


def check_gene_validity(gene):
    """
    check if a gene is valid.
    :param gene: string or int
        specifying the acronym or ID of the gene
    :return: boolean
    """
    if isinstance(gene, int):  # if gene id is provided
        api_query_criteria = 'criteria=products%5Bid$eq1%5D,' \
                             'genes%5Bid$eq%27{}%27%5D'.format(gene)
    else:  # if gene acronym is provided
        api_query_criteria = 'criteria=products%5Bid$eq1%5D,' \
                             'genes%5Bacronym$eq%27{}%27%5D'.format(gene)
    # extra conditions
    api_query_include = '&include=genes'

    # make the query
    query_url = API_QUERY_STRING + \
        api_query_criteria + \
        api_query_include
    r = requests.get(query_url)
    root = ET.fromstring(r.content)

    return True if root.text else False


def get_gene_info(gene, attributes='all'):
    """
    get attributes associated with a gene

    Parameters
    ----------
    gene: string or int
        specifying the acronym or ID of the gene
    attributes: list
        specifying which attribute of the gene to return
        default: 'all' (return all the attributes in a dict)
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
    gene_info: dict {attribute:value}

    """
    # how many attributes are specified
    if attributes == 'all':
        attr_list = GENE_ATTRIBUTES
    elif isinstance(attributes, list):
        attr_list = attributes

    if check_gene_validity(gene) is False:
        raise ValueError('Gene {} is invalid'.format(gene))

    if isinstance(gene, int):  # if gene id is provided
        api_query_criteria = 'criteria=products%5Bid$eq1%5D,' \
                             'genes%5Bid$eq%27{}%27%5D'.format(gene)
    else:  # if gene acronym is provided
        api_query_criteria = 'criteria=products%5Bid$eq1%5D,' \
                             'genes%5Bacronym$eq%27{}%27%5D'.format(gene)
    # extra conditions
    api_query_include = '&include=genes'
    # make the query
    query_url = API_QUERY_STRING + \
        api_query_criteria + \
        api_query_include
    r = requests.get(query_url)
    root = ET.fromstring(r.content)

    # initialize a dict to store the attributes
    gene_info = dict()
    if attributes == 'all' or isinstance(attributes, list):
        for attr in attr_list:
            try:
                gene_info[attr] = _get_single_gene_attribute(
                    root, attr
                )
            except AttributeError:
                print('There is no attribute called {}. '
                      'Skipped.'.format(attr))
                continue

    else:  # single attribute is given
        return _get_single_gene_attribute(root, attributes)

    return gene_info


def _get_single_gene_attribute(root, attr):
    """
    get the value of single gene attribute
    :param root: Element 'Response' to parse
    :param attr: str, attribute to return
    :return: str or int, the value of the attribute
    """
    item = root.find(
        'section-data-sets/section-data-set/'
        'genes/gene/{}'.format(attr)
    )
    # check if attr is valid
    if not item:
        raise AttributeError(
            'There is no attribute called {}'.format(attr)
        )

    # check data type
    if 'type' in item.attrib \
            and item.attrib['type'] == 'integer':
        return int(item.text)
    else:  # attribute is a string
        return item.text






