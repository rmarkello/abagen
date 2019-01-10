# -*- coding: utf-8 -*-
"""
functions to fetch unionization relevant data
"""
import requests
from xml.etree import ElementTree as ET
import numpy as np
import pandas as pd
from .gene import check_gene_validity
from .io import read_all_structures

URL_PREFIX = "http://api.brain-map.org/api/v2/data/" \
             "SectionDataSet/query.xml?"
URL_INCLUDE = "&include=structure_unionizes%28structure%29"

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
        the list of structures, either id (int) or acronym (str).
        default: structures as documented in Rubinov et al, 2015
    attributes: list, optional
        specify the unionization data attributes to include
        default: 'all'
        available attributes:
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

    Returns
    -------
    unionization: dict or numpy_array
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
    query_url = URL_PREFIX + \
        'id={}'.format(experiment_id) + \
        URL_INCLUDE
    print('accessing {}...'.format(query_url))
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

    # if no attribute is given
    if attributes == 'all':
        attr_list = UNIONIZATION_ATTRIBUTES
    else:
        attr_list = attributes
            
    # only one attribute is specified
    if isinstance(attr_list, str):
        return _get_single_unionization_attribute(
                root, attributes, path_prefix, structure_list
            )  # or raise AttributError
    else:  # if multiple attributes are given
        for attr in attr_list:
            print('fetching {} unionization data...'.format(attr))
            try:
                unionization[attr] = _get_single_unionization_attribute(
                    root, attr, path_prefix, structure_list
                )
            except AttributeError:
                print('There is no attribute called {}. '
                      'Skipped. '.format(attr))
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
    value_items = root.findall(path_prefix + attr)
    if not value_items:  # items is empty
        raise AttributeError(
            'There is no unionization attribute called {}.'
            .format(attr)
        )

    # if structure acronym is given
    if isinstance(structure_list[0], str):
        structure_items = root.findall(path_prefix + 'structure/acronym')
        # extract structure and values in structure_list
        all_items = [
            (float(val_item.text), structure_item.text)
            for val_item, structure_item in zip(
                value_items,
                structure_items
            ) if structure_item.text in structure_list
        ]
    elif isinstance(structure_list[0], int):
        structure_items = root.findall(path_prefix + 'structure/id')
        # extract structure and values in structure_list
        all_items = [
            (float(val_item.text), int(structure_item.text))
            for val_item, structure_item in zip(
                value_items,
                structure_items
            ) if int(structure_item.text) in structure_list
        ]
    else:
        raise ValueError('structures are invalid. '
                         'Please use ID (integers) or '
                         'acronym (strings) instead.')

    structures_in_the_list = [item[1] for item in all_items]
    vals_in_the_list = np.array(
        [item[0] for item in all_items],
        dtype=np.float
    )

    vals = np.empty((roi_count,))
    vals[:] = np.nan
    # average duplicate unionizations
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
                  'are found in region {1}. '
                  'Set to NaN.'
                  .format(attr, structure_list[k]))

    return vals


def get_experiment_id_from_gene(
        gene, entry_type, slicing_direction='sagittal'
):
    """
    fetches mouse unionization data of a single experiment,
    either saggital or coronal, according to the experiment id

    Parameters
    ----------
    gene: str or int
        specifying the acronym (capitalized) or ID of the gene
    entry_type: str
        the type of gene identifier.
        supported:
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

    if check_gene_validity(gene, entry_type) is False:
        raise ValueError('{0} {1} is invalid'.format(entry_type, gene))

    # find the experiment IDs associated with this gene
    if isinstance(gene, int):  # if integer is provided
        api_query_criteria = "criteria=products[id$eq1]," \
                             "genes[{0}$eq{1}]," \
                             "plane_of_section[name$eq'{2}']"\
                             .format(entry_type, gene, slicing_direction)
    else:  # if string is provided
        api_query_criteria = "criteria=products[id$eq1]," \
                             "genes[{0}$eq'{1}']," \
                             "plane_of_section[name$eq'{2}']" \
                             .format(entry_type, gene, slicing_direction)

    # extra conditions
    api_query_include = "&include=genes"
    
    # make the query
    query_url = URL_PREFIX + \
        api_query_criteria + \
        api_query_include
    print('accessing {}...'.format(query_url))
    r = requests.get(query_url)
    root = ET.fromstring(r.content)

    experiment_id_list = []
    for item in root.findall('section-data-sets/section-data-set/id'):
        experiment_id_list.append(int(item.text))  # append experiment id

    if not experiment_id_list:  # no experiments are found
        print('No {0} experiments are found for gene {1}, '
              'return an empty experiment ID list.'
              .format(slicing_direction, gene))

    return experiment_id_list


def get_unionization_from_gene(
        gene,
        entry_type,
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
    entry_type: str
        the type of gene identifier.
        supported:
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
        gene, entry_type, slicing_direction=slicing_direction
    )

    # initialize a dict to store attr values and valid experimentIDs
    unionization = dict()

    if attributes == 'all':
        attr_list = UNIONIZATION_ATTRIBUTES
    else:
        attr_list = attributes

    if isinstance(attributes, str):
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

    else:
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
            unionization[attr] = attr_vals

    return unionization
