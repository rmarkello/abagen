# -*- coding: utf-8 -*-
"""
functions to fetch unionization relevant data
"""
import requests
from xml.etree import ElementTree as ET
import numpy as np
from .gene import check_gene_validity
from .io import read_all_structures
import random

URL_PREFIX = "http://api.brain-map.org/api/v2/data/" \
    "SectionDataSet/query.xml?"
URL_INCLUDE = "&include=structure_unionizes%28structure%29"

PATH_PREFIX = "section-data-sets/section-data-set/" \
    "structure-unionizes/structure-unionize/"

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


def get_experiment_id_from_gene(
    gene_id=None,
    gene_acronym=None,
    gene_name=None,
    slicing_direction='sagittal'
):
    """
    fetches mouse unionization data of a single experiment,
    either saggital or coronal, according to the experiment id

    Parameters
    ----------
    gene_id : int, optional
    gene_acronym : str, optional
    gene_name : str, optional
        at least one gene identifier should be given.
    slicing_direction : str, optional
        slicing scheme of the samples, 'sagittal' or 'coronal'.
        default: 'sagittal'

    Returns
    -------
    experiment_id_list : list
        list of integers (experiment IDs in slicing_direction) of the gene

    Raises
    ------
    ValueError
        slicing_direction is invalid, or
        gene identifier is invalid

    Examples
    --------
    >>> get_experiment_id_from_gene(
    >>>     gene_id=84193, slicing_direction='coronal'
    >>> )
    [74047443]

    """
    if slicing_direction not in ['sagittal', 'coronal']:
        raise ValueError('Slicing slicing_direction {} is invalid. '
                         'Try sagittal or coronal instead'
                         .format(slicing_direction))

    validity, root = check_gene_validity(
        gene_id=gene_id,
        gene_acronym=gene_acronym,
        gene_name=gene_name
    )

    if validity is False:
        raise ValueError(
            'the gene {} is invalid'
            .format(
                [
                    item
                    for item in [gene_id, gene_acronym, gene_name]
                    if item is not None
                ][0]
            )
        )
    experiment_id_list = [
        int(item.text)
        for item, plane_of_section in zip(
            root.findall(
                'section-data-sets/section-data-set/id'
            ), root.findall(
                'section-data-sets/section-data-set/plane-of-section/name'
            )
        ) if plane_of_section.text == slicing_direction
    ]

    if len(experiment_id_list) == 0:  # no experiments are found
        print(
            'No {0} experiments are found for gene {1}, '
            'return an empty experiment ID list.'
            .format(
                slicing_direction, [
                    item for item in [gene_id, gene_acronym, gene_name]
                    if item is not None
                ]
            )
        )

    return experiment_id_list


# to be tested
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
    experiment_id : int
        specifying the experiment id
    structure_list : list, optional
        the list of structures, either id (int) or acronym (str).
        default: structures as documented in Rubinov et al, 2015
        we recommend to use ID against acronym
        as acronyms may not be unique
        (e.g., 'PH' have names 'pontine hindbrain' in mouse atlas
        and 'pontine hindbrain (pons proper)' in developing mouse atlas)
    attributes : list, optional
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
    unionization : dict or numpy_array
        unionization data attributes and the values {attribute:value}

    Raises
    ------
    ValueError
        If experiment_id is invalid
    """

    if structure_list is None:
        # read default structure list
        # see Rubinov et al 2015 for more information
        structure_list = read_all_structures(entry_type='id')

    # make the query
    query_url = URL_PREFIX + \
        'id={}'.format(experiment_id) + \
        URL_INCLUDE
    print('accessing {}...'.format(query_url))
    r = requests.get(query_url)
    root = ET.fromstring(r.content)

    if root.attrib['total_rows'] == '0':  # check if any expressions are found
        raise ValueError(
            'No unionization data is found found '
            'in experiment {}. '
            'Try another valid experiment ID'
            .format(experiment_id)
        )

    # if no attribute is given
    if attributes == 'all':
        attr_list = UNIONIZATION_ATTRIBUTES
    else:
        attr_list = attributes

    # only one attribute is specified
    if isinstance(attr_list, str):
        return _get_single_unionization_attribute(
                root, attr_list, structure_list
            )  # or raise AttributError
    else:  # if multiple attributes are given
        unionization = dict()
        for attr in attr_list:
            print('fetching {} unionization data...'.format(attr))
            try:
                unionization[attr] = _get_single_unionization_attribute(
                    root, attr, structure_list
                )
            except AttributeError:
                print('There is no attribute called {}. '
                      'Skipped. '.format(attr))
                continue

    return unionization


def _get_single_unionization_attribute(root, attr, structure_list):
    """
    Parameters
    ----------
    root: obj, ElementTree 'Response'
    attr: str
        the attribute to return
    structure_list: array_like
        the list of structures (ID or acronym)

    Returns
    -------
    None or numpy_array
        if attr exists but has no associated values, None is returned.
        if attr does not exist, raise AttributeError
        if attr exits and has values, return (N, ) numpy_array corresponding to
        the structures in structure_list
        if structure_list is in valid format, but certain structures have no
        unionization data, the values will be set to NaN

    Raises
    ------
    ValueError:
        structure_list is invalid
    AttributeError:
        Unionization attribute given is invalid
    """
    roi_count = len(structure_list)

    # if the attribute does not exist, raise AttributeError
    value_items = root.findall(PATH_PREFIX + attr)
    if len(value_items) == 0:  # val_items is empty
        raise AttributeError(
            'There is no unionization attribute called {}.'
            .format(attr)
        )

    try:
        # if structure acronym is given
        if isinstance(random.choice(structure_list), str):
            structure_items = root.findall(PATH_PREFIX + 'structure/acronym')
            # extract structure and values in structure_list
            all_items = [
                (float(val_item.text), structure_item.text)
                for val_item, structure_item in zip(
                    value_items,
                    structure_items
                ) if structure_item.text in structure_list
            ]
        elif isinstance(random.choice(structure_list), (np.integer, int)):
            structure_items = root.findall(PATH_PREFIX + 'structure/id')
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
    except TypeError:  # when float(val_item.text) return TypeError
        # namely, the attribute given is valid, but has no values
        # return an empty numpy_array
        return np.empty((roi_count, 0))

    vals_in_the_list = np.array(
        [item[0] for item in all_items], dtype=np.float
    )
    structures_in_the_list = [item[1] for item in all_items]

    vals = np.empty((roi_count,))
    vals[:] = np.nan  # initialize the values with NaNs

    # average duplicate unionizations,
    # especially when structure acronyms are given
    # unlikely to happen if structure id is given
    for k, structure in enumerate(structure_list):
        index = [
            idx
            for idx, item
            in enumerate(structures_in_the_list)
            if item == structure
        ]
        # if gene expressions are found in kth region
        if len(index) != 0:
            vals[k] = vals_in_the_list[index].mean()
        else:
            print('No {0} values '
                  'are found in region {1}. '
                  'Set to NaN.'
                  .format(attr, structure))

    return vals


def get_unionization_from_gene(
    gene_id=None,
    gene_acronym=None,
    gene_name=None,
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
    gene_id: int, optional
    gene_acronym: str, optional
    gene_name: str, optional
        at least one gene identifier should be given.
    slicing_direction: str, optional
        slicing scheme of the samples, 'sagittal' or 'coronal'.
        default: 'sagittal'
    structure_list: a list of int or str, optional
        the list of structure IDs or acronyms.
        default: structures as documented in Rubinov et al, 2015
    attributes: str or list, optional
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
        if single attribute is given, return a numpy_array
        unionization data attributes and the values {attribute:value}
    """
    if slicing_direction not in ['sagittal', 'coronal']:
        raise ValueError("Slicing slicing_direction {} is invalid. "
                         "Try 'sagittal' or 'coronal' instead."
                         .format(slicing_direction))

    if structure_list is None:
        # read default structure list
        # (see Rubinov et al, 2015 for the criteria of
        # choosing the ROIs)
        structure_list = read_all_structures(entry_type='id')
    roi_count = len(structure_list)

    experiment_id_list = get_experiment_id_from_gene(
        gene_id=gene_id,
        gene_acronym=gene_acronym,
        gene_name=gene_name,
        slicing_direction=slicing_direction,
    )

    if attributes == 'all':
        attr_list = UNIONIZATION_ATTRIBUTES
    else:
        attr_list = attributes

    if isinstance(attr_list, str):  # if single attribute is given
        attr_vals = np.empty((0, roi_count))
        for experiment_id in experiment_id_list:
            vals = get_unionization_from_experiment(
                experiment_id,
                structure_list,
                attributes=attr_list
            )  # can raise AttributeError, or return None
            attr_vals = np.append(
                attr_vals,
                vals.reshape(-1, roi_count),
                axis=0
            )
        return attr_vals
    else:  # multiple attributes
        # initialize a dict to store unionization data
        unionization = dict.fromkeys(
            attr_list, np.empty((0, roi_count))
        )
        for experiment_id in experiment_id_list:
            # get unionization dict of a single experiment
            # then concatenate each elements
            try:
                unionization_tmp = get_unionization_from_experiment(
                    experiment_id,
                    structure_list,
                    attributes=attr_list
                )
                for item in unionization_tmp:
                    unionization[item] = np.append(
                        unionization[item],
                        unionization_tmp[item].reshape(-1, roi_count),
                        axis=0
                    )
            except ValueError:
                print('No unionization data is found '
                      'in experiment {}. skipped'
                      .format(experiment_id))
        # remove the nan arrays
        unionization = {
            item: unionization[item]
            for item in unionization
            if unionization[item].size != 0
        }
    return unionization
