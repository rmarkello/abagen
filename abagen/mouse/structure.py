# -*- coding: utf-8 -*-
"""
functions to make mouse structure queries and manipulations
such as structure ID <-> structure acronym conversions,
fetch structure full names,
find structure coordinates, etc
"""
import requests
from xml.etree import ElementTree as ET

URL_PREFIX = "http://api.brain-map.org/api/v2/data/" \
    "query.xml?include=model::Structure"
# Allen mouse atlas
URL_INCLUDE_MOUSE = ",structure_centers," \
    "ontology[id$eq1]"
# Allen developing mouse atlas
URL_INCLUDE_DEV_MOUSE = ",structure_centers," \
    "ontology[id$eq12]"
# for more available attributes, see
# http://api.brain-map.org/doc/Structure.html
# more attributes to be incorporated

STRUCTURE_ENTRY_TYPES = [
    'acronym',
    'id',
    'name'
]

STRUCTURE_ATTRIBUTES = [
    'acronym',
    'atlas-id',
    'color-hex-triplet',
    'depth',
    'graph-id',
    'graph-order',
    'hemisphere-id',
    'id',
    'name',
    'neuro-name-structure-id',
    'neuro-name-structure-id-path',
    'ontology-id',
    'parent-structure-id',
    'safe-name',
    'sphinx-id',
    'structure-id-path',
    'weight'
]


def check_structure_validity(
    structure_id=None,
    structure_acronym=None,
    structure_name=None,
):
    """
    check if a structure is valid or has records in the database
    Example
    -------
    # check if structure ID 1018 is valid
    validity, _, _ = check_structure_validity(structure_id=1018)
    # check if structure SSp is valid
    validity, _, _ = check_structure_validity(structure_acronym='SSp')

    Parameters
    ----------
    structure_id: int, optional
        structure ID
    structure_acronym: str, optional
        structure acronym (case sensitive)
    structure_name: str, optional
        structure name (case sensitive)

    Returns
    -------
    validity: boolean
        if the structure has records in the database
    root_mouse: obj,
        Element 'Response' object, empty if query fails
    root_dev_mouse: obj
        Element 'Response' object, empty if query fails

    Raises
    ------
    TypeError:
        missing parameters
    """
    # if structure ID is given
    # preferred: id > acronym > name
    if structure_id is not None:
        query_url_mouse = URL_PREFIX + \
            "[id$eq{}]".format(structure_id) + \
            URL_INCLUDE_MOUSE
        query_url_dev_mouse = URL_PREFIX + \
            "[id$eq{}]".format(structure_id) + \
            URL_INCLUDE_DEV_MOUSE
    # then if acronym is given
    elif structure_acronym is not None:
        query_url_mouse = URL_PREFIX + \
            "[acronym$eq'{}']".format(structure_acronym) + \
            URL_INCLUDE_MOUSE
        query_url_dev_mouse = URL_PREFIX + \
            "[acronym$eq'{}']".format(structure_acronym) + \
            URL_INCLUDE_DEV_MOUSE
    # then if name is given
    elif structure_name is not None:
        query_url_mouse = URL_PREFIX + \
            "[name$eq'{}']".format(structure_name) + \
            URL_INCLUDE_MOUSE
        query_url_dev_mouse = URL_PREFIX + \
            "[name$eq'{}']".format(structure_name) + \
            URL_INCLUDE_DEV_MOUSE
    else:
        # if no input
        raise TypeError(
            "at least one structure identifier should be specified"
        )

    print('accessing {}...'.format(query_url_mouse))
    r_mouse = requests.get(query_url_mouse)
    root_mouse = ET.fromstring(r_mouse.content)

    print('accessing {}...'.format(query_url_dev_mouse))
    r_dev_mouse = requests.get(query_url_dev_mouse)
    root_dev_mouse = ET.fromstring(r_dev_mouse.content)

    # both are empty
    if root_mouse.attrib['total_rows'] == '0' \
            and root_dev_mouse.attrib['total_rows'] == '0':
        return False, root_mouse, root_dev_mouse
    else:  # at least one success
        return True, root_mouse, root_dev_mouse


def get_structure_info(
    structure_id=None,
    structure_acronym=None,
    structure_name=None,
    attributes='all'
):
    """
    get attributes of a structure
    Parameters
    ----------
    structure_id: int, optional
        structure ID
    structure_acronym: str, optional
        structure acronym (case sensitive)
    structure_name: str, optional
        structure name (case sensitive)
    attributes: str or list, optional
        a single attribute or a list of attributes
        default: 'all', returning all the attributes
        available attributes:
            'acronym',
            'atlas-id',
            'color-hex-triplet',
            'depth',
            'graph-id',
            'graph-order',
            'hemisphere-id',
            'id',
            'name',
            'neuro-name-structure-id',
            'neuro-name-structure-id-path',
            'ontology-id',
            'parent-structure-id',
            'safe-name',
            'sphinx-id',
            'structure-id-path',
            'weight'

    Returns
    -------
    structure_info: list or dict
        if a single attribute is given, return a list of values
        (int or str) of that attribute acquired from mouse atlas and
        developing mouse atlas
        if multiple attributes are given, return a dict
        {attr:values}. attr is str (attribute) and values is list

    Raises
    ------
    ValueError:
        the structure given is invalid
    """
    validity, root1, root2 = check_structure_validity(
        structure_id=structure_id,
        structure_acronym=structure_acronym,
        structure_name=structure_name
    )

    if validity is False:
        raise ValueError('the structure given is invalid')

    if attributes == 'all':
        attr_list = STRUCTURE_ATTRIBUTES
    else:
        attr_list = attributes

    if isinstance(attr_list, str):
        # attr_list is a str (single attribute)
        # return the single attr value, or AttributeError
        try:
            structure_info = [
                _get_single_structure_attribute(root, attr_list)
                for root in [root1, root2]
                if root.attrib['total_rows'] != '0'
            ]
            return structure_info
        except AttributeError:
            print('There is no attribute called {}. '
                  .format(attr_list))
            return []
    else:
        structure_info = dict()
        # iterate through the attributes
        for attr in attr_list:
            try:
                # extract the info of attr
                structure_info[attr] = [
                    _get_single_structure_attribute(root, attr)
                    for root in [root1, root2]
                    if root.attrib['total_rows'] != '0'
                ]
            except AttributeError:
                print('There is no attribute called {}. '
                      'Skipped.'.format(attr))
                continue

    return structure_info


def _get_single_structure_attribute(root, attr):
    """
    return single attribute of a structure.

    Parameters
    ----------
    root: obj, Element 'Response'
    attr: str
        the attribute to return

    Returns
    -------
    int or str, the value of structure's attr
    """
    item = root.findall(
        'structures/structure/{}'
        .format(attr)
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
        # has this attribute, but no value in the database
        return None


def get_structure_coordinates(
    structure_id=None,
    structure_acronym=None,
    structure_name=None
):
    """
    get structure coordinates in reference atlas.

    Examples
    --------
    # get the coordinates of structure ID 1018
    coor = get_structure_coordinates(structure_id=1018)
    # if the structure has records in mouse brain atlas
    # coordinates in space with reference id 10
    # (Mouse brain atlas left hemisphere)
    coordinates = coor[10]
    # coordinates in space with reference id 9
    # (Mouse brain atlas right hemisphere)
    coordinates = coor[9]
    # get the coordinates of structure acronym SSp
    coor = get_structure_coordinates(structure_acronym='CA')
    coor

    Parameters
    ----------
    structure_id: int, optional
        structure ID
    structure_acronym: str, optional
        structure acronym (case sensitive)
    structure_name: str, optional
        structure name (case sensitive)

    Returns
    -------
    coor: dict {int: list of (3, } tuple}
        {reference-space-id:[(x, y, z)]}
    """

    validity, root1, root2 = check_structure_validity(
        structure_id=structure_id,
        structure_acronym=structure_acronym,
        structure_name=structure_name
    )

    if validity is False:
        raise ValueError('the structure given is invalid')

    coor = dict()
    # find coor in mouse atlas and developing mouse atlas
    # in the form {reference-space-id:(x, y, z)}
    for root in [root1, root2]:
        if root.attrib['total_rows'] != '0':  # if root is not empty
            for item in root.findall(
                    'structures/structure/'
                    'structure-centers/structure-center'
            ):
                reference_atlas_id = int(
                    item.find('reference-space-id').text
                )
                if reference_atlas_id in coor:
                    # append a tuple (x, y, z)
                    coor[reference_atlas_id].append(
                        (
                            int(item.find('x').text),
                            int(item.find('y').text),
                            int(item.find('z').text)
                        )
                    )
                else:  # create a new list for reference_atlas_id
                    coor[reference_atlas_id] = [
                        (
                            int(item.find('x').text),
                            int(item.find('y').text),
                            int(item.find('z').text)
                        )
                    ]

    if len(coor) == 0:
        print('No coordinates information is found')

    return coor
