# -*- coding: utf-8 -*-

import requests
from xml.etree import ElementTree as ET
import numpy as np
import pandas as pd

API_QUERY_STRING = 'http://api.brain-map.org/' \
                   'api/v2/data/SectionDataSet/query.xml?'
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


def get_expression_from_experiment(
        experiment_id,
        roi_list=None,
):
    """
    fetches mouse gene expression data of a single experiment,
    either saggital or coronal, according to the experiment id

    Parameters
    ----------
    experiment_id: an integer, specifying the experiment id
    roi_list: a list of strings, optional
        the list of ROIs in the form of acronyms.
        default: ROIs as included in Rubinov et al, 2015

    Returns
    -------
    gene_epr: (N,) numpy.array
        regional gene expressions corresponding to
        the ROIs in structure_list

    Raises
    ------
    ValueError:
        If experiment_id is valid
    """
    if not isinstance(experiment_id, int):
        raise ValueError(
            '{} is not a valid experiment ID'.format(experiment_id)
        )

    if roi_list is None:
        # read default ROI list
        # (see Rubinov et al, 2015 for the criteria of
        # choosing the ROIs)
        roilabels = pd.read_csv(
            "abagen/data/roilabels-rubinov2015pnas.csv"
        )
        roi_list = roilabels['roiacronyms']

    roi_count = len(roi_list)

    all_structure = []
    all_epr = np.empty((0, 1))

    # make the query
    api_query_criteria = 'id=' + str(experiment_id)
    api_query_include = '&include=structure_unionizes%28structure%29'
    query_url = API_QUERY_STRING + \
        api_query_criteria + \
        api_query_include
    r = requests.get(query_url)
    root = ET.fromstring(r.content)

    if not root.text:  # check if any expressions are found
        raise ValueError(
            'No gene expression values are found in experiment {}. '
            'Try another experiment ID'
            .format(experiment_id)
        )

    for item in root.findall('section-data-sets/section-data-set/'
                             'structure-unionizes/structure-unionize'):
        # append new epr value
        for subitem in item.findall('expression-energy'):
            all_epr = np.append(all_epr, float(subitem.text))
        # append new structure label
        for subitem in item.findall('structure/acronym'):
            all_structure.append(subitem.text)

    # extract the regions in roi_list
    data_in_the_list = [
        (item, all_epr[k])
        for k, item
        in enumerate(all_structure)
        if item in roi_list
    ]
    structure_in_the_list = [item[0] for item in data_in_the_list]
    epr_in_the_list = np.array(
        [item[1] for item in data_in_the_list]
    )

    gene_epr = np.empty((roi_count, ))
    gene_epr[:] = np.nan
    # average duplicate expressions
    for k in range(roi_count):
        index = [
            idx
            for idx, item
            in enumerate(structure_in_the_list)
            if item == roi_list[k]
        ]
        if not index:  # if no gene expressions are found, leave it to NaN
            print('No gene expression values are found for region {}. '
                  'Set to NaN.'
                  .format(roi_list[k]))
            continue
        gene_epr[k] = epr_in_the_list[index].mean()

    return gene_epr


def get_experiment_id_from_gene(gene, direction='sagittal'):
    """
    fetches mouse gene expression data of a single experiment,
    either saggital or coronal, according to the experiment id

    Parameters
    ----------
    gene: str or int
        specifying the acronym or ID of the gene
    direction: str, optional
        slicing scheme of the samples, 'sagittal' or 'coronal'.
        default: 'sagittal'

    Returns
    -------
    experiment_id_list: list
        integers (experiment IDs) of the gene specified by gene_acronym
    """
    if direction not in ['sagittal', 'coronal']:
        raise ValueError('Slicing direction {} is invalid. '
                         'Try sagittal or coronal instead'
                         .format(direction))

    if check_gene_validity(gene) is False:
        raise ValueError('Gene {} is invalid'.format(gene))

    # find the experiment IDs associated with this gene
    if isinstance(gene, int):  # if gene id is provided
        api_query_criteria = 'criteria=products%5Bid$eq1%5D,' \
                         'genes%5Bid$eq%27' + str(gene)
    else:  # if gene acronym is provided
        api_query_criteria = 'criteria=products%5Bid$eq1%5D,' \
                         'genes%5Bacronym$eq%27' + gene
    # extra conditions
    api_query_include = '%27%5D,plane_of_section[name$eq%27' + \
                        direction + \
                        '%27]&include=genes,section_images'
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
              'return an empty list'
              .format(direction, gene))

    return experiment_id_list


def get_expression_from_gene(
        gene,
        direction='sagittal',
        roi_list=None,
):
    """
    fetches mouse gene expression data
    averaged across all the associated experiments
    (either sagittal or coronal) of a specific gene

    Parameters
    ----------
    gene: string or int
        specifying the acronym or ID of the gene
    direction: string, optional
        slicing scheme of the samples, 'sagittal' or 'coronal'.
        default: 'sagittal'
    roi_list: a list of strings, optional
        the list of ROIs in the form of acronyms.
        default: ROIs as included in Rubinov et al, 2015

    Returns
    -------
    gene_epr: (N,) numpy.array
        regional gene expressions corresponding to
        the ROIs in structure_list
    valid_experiment_id_list: list or str
        the list of integers (experiment ID)
        corresponding to the rows of gene_epr
    """
    if direction not in ['sagittal', 'coronal']:
        raise ValueError('Slicing direction {} is invalid. '
                         'Try sagittal or coronal instead'
                         .format(direction))

    if roi_list is None:
        # read default ROI list
        # (see Rubinov et al, 2015 for the criteria of
        # choosing the ROIs)
        roilabels = pd.read_csv(
            "abagen/data/roilabels-rubinov2015pnas.csv"
        )
        roi_list = roilabels['roiacronyms']
    roi_count = len(roi_list)

    experiment_id_list = get_experiment_id_from_gene(
        gene, direction=direction
    )

    gene_epr = np.empty((0, roi_count))
    valid_experiment_id_list = []
    for experiment_id in experiment_id_list:
        try:
            gene_epr = np.append(
                gene_epr,
                get_expression_from_experiment(
                    experiment_id,
                    roi_list,
                ),
                axis=0
            )
            valid_experiment_id_list.append(experiment_id)
        except ValueError:
            print('No gene expression values are found in experiment {}'
                  .format(experiment_id))
            continue

    return gene_epr, valid_experiment_id_list


def check_gene_validity(gene):
    """
    check if a gene is valid.
    :param gene: string or int
        specifying the acronym or ID of the gene
    :return: boolean
    """
    if isinstance(gene, int):  # if gene id is provided
        api_query_criteria = 'criteria=products%5Bid$eq1%5D,' \
                         'genes%5Bid$eq%27' + str(gene)
    else:  # if gene acronym is provided
        api_query_criteria = 'criteria=products%5Bid$eq1%5D,' \
                         'genes%5Bacronym$eq%27' + gene
    # extra conditions
    api_query_include = '%27%5D&include=genes'
    # make the query
    query_url = API_QUERY_STRING + \
        api_query_criteria + \
        api_query_include
    r = requests.get(query_url)
    root = ET.fromstring(r.content)

    return True if root.text else False


def get_gene_info(gene):
    """
    get attributes associated with a gene

    Parameters
    ----------
    gene: string or int
        specifying the acronym or ID of the gene

    Returns
    -------
    gene_info: dict

    """
    if check_gene_validity(gene) is False:
        raise ValueError('Gene {} is invalid'.format(gene))

    if isinstance(gene, int):  # if gene id is provided
        api_query_criteria = 'criteria=products%5Bid$eq1%5D,' \
                         'genes%5Bid$eq%27' + str(gene)
    else:  # if gene acronym is provided
        api_query_criteria = 'criteria=products%5Bid$eq1%5D,' \
                         'genes%5Bacronym$eq%27' + gene
    # extra conditions
    api_query_include = '%27%5D&include=genes'
    # make the query
    query_url = API_QUERY_STRING + \
        api_query_criteria + \
        api_query_include
    r = requests.get(query_url)
    root = ET.fromstring(r.content)

    gene_info = dict()

    for attr in GENE_ATTRIBUTES:
        item = root.find(
            'section-data-sets/section-data-set/genes/gene/{}'.format(attr)
        )
        # check data type
        if 'type' in item.attrib and item.attrib['type'] == 'integer':
            gene_info[attr] = int(item.text)
        else:  # attribute is a string
            gene_info[attr] = item.text

    return gene_info


