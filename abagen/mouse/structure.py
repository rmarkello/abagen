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

URL_INCLUDE = ",structure_centers," \
              "ontology[abbreviation$eq'Mouse']"
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
    acronym=None,
    name=None,
):
    """
    check if a structure is valid or has records in the database
    Example
    -------
    # check if structure ID 1018 is valid
    validity, root = check_structure_validity(id=1018)
    # check if structure SSp is valid
    validity, root = check_structure_validity(acronym='SSp')

    Parameters
    ----------
    structure_id: int, optional
        structure ID
    acronym: str, optional
        structure acronym (case sensitive)
    name: str, optional
        structure name (case sensitive)

    Returns
    -------
    validity: boolean
        if the structure has records in the database
    root: obj,
        ElementTree 'Response' object

    Raises
    ------
    TypeError:
        missing parameters
    """
    # if structure ID is given
    if structure_id is not None:
        query_url = URL_PREFIX + \
            "[id$eq{}]".format(structure_id) + \
            URL_INCLUDE
    elif acronym is not None:
        query_url = URL_PREFIX + \
            "[acronym$eq'{}']".format(acronym) + \
            URL_INCLUDE
    elif name is not None:
        query_url = URL_PREFIX + \
            "[name$eq'{}']".format(name) + \
            URL_INCLUDE
    else:
        # if no input
        raise TypeError(
            "at least one structure identifier should be specified"
        )

    print('accessing {}...'.format(query_url))

    r = requests.get(query_url)
    root = ET.fromstring(r.content)

    return True, root if root.text else False, root


def get_structure_info(
    structure_id=None,
    acronym=None,
    name=None,
    attributes='all'
):
    """
    get attributes of a structure
    Parameters
    ----------
    structure_id: int, optional
        structure ID
    acronym: str, optional
        structure acronym (case sensitive)
    name: str, optional
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
    structure_info: int, str, or dict
        if a single attribute is given, return int or str
        if multiple attributes are given, return a dict

    Raises
    ------
    ValueError:
        the structure given is invalid
    """
    validity, root = check_structure_validity(
        structure_id=structure_id, acronym=acronym, name=name
    )

    if validity is False:
        raise ValueError('the structure given is invalid')

    if attributes == 'all':
        attr_list = STRUCTURE_ATTRIBUTES
    else:
        attr_list = attributes

    structure_info = dict()

    if isinstance(attr_list, str):
        # attr_list is a str (single attribute)
        # return the single attr value, or AttributeError
        try:
            return _get_single_structure_attribute(root, attr_list)
        except AttributeError:
            print('There is no attribute called {}. '
                  .format(attr_list))
            # then return an empty dict...
    else:
        # iterate through the attributes
        for attr in attr_list:
            try:
                # extract the info of attr
                structure_info[attr] = _get_single_structure_attribute(
                    root, attr
                )
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
    root: obj, ElementTree obj
    attr: str
        the attribute to return

    Returns
    -------
    int or str, the value of structure's attr
    """
    item = root.find(
        'structures/structure/{}'
        .format(attr)
    )
    # check if attr is valid (if any information is found)
    if not item:
        raise AttributeError(
            'There is no gene attribute called {}'.format(attr)
        )

    # check data type
    # attr is an integer
    try:
        return int(item.text)
    except ValueError:
        return item.text


def get_structure_coordinates(
    structure_id=None,
    acronym=None,
    name=None
):
    """
    get structure coordinates in reference atlas.

    Examples
    --------
    # get the coordinates of structure ID 1018
    coor = get_structure_coordinates(structure_id=1018)
    # get the coordinates of structure acronym SSp
    coor = get_structure_coordinates(acronym='SSp')

    Parameters
    ----------
    structure_id: int, optional
        structure ID
    acronym: str, optional
        structure acronym (case sensitive)
    name: str, optional
        structure name (case sensitive)

    Returns
    -------
    coor: list
        a list of (3, ) tuples specifying (x, y, z) position
        if len(list) == 2, coor[0] is left hemisphere
        and coor[1] is right hemisphere

    """

    validity, root = check_structure = check_structure_validity(
        structure_id=structure_id,
        acronym=acronym,
        name=name
    )

    if validity is False:
        raise ValueError('the structure given is invalid')

    coor = []
    for item in root.findall(
            'structures/structure/'
            'structure-centers/structure-center'
    ):
        coor.append(
            (
                int(item.find('x').text),
                int(item.find('y').text),
                int(item.find('z').text)
            )
        )

    if len(coor) == 0:
        print('No coordinates information is found')

    return coor
