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
    'name',
    'safe-name',
    'sphinx-id',
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


def check_structure_validity(structure, entry_type='id'):
    """
    check if a structure id or acronym is valid
    :param structure: str or int
        structure identifier, either an int (id, sphinx-id),
        or a str (acronym, name, safe-name, etc.), specified
        in entry_type.
        recommend to use structure ID (int) or
        acronym (str, case sensitive).
    :param entry_type: str, optional
        the type of structure identifier. default: 'id'
        supported:
            'acronym',
            'id',
            'name',
            'safe-name',
            'sphinx-id'

    :return: tuple or boolean
        if the structure is invalid, return False
        else return a tuple (True, obj)
        where obj is the ET 'response' object to parse
    """
    if entry_type not in STRUCTURE_ENTRY_TYPES:
        raise ValueError('entry_type {} is invalid'
                         .format(entry_type))

    # URL
    if isinstance(structure, int):
        query_url = URL_PREFIX + \
                    "[{0}$eq{1}]".format(entry_type, structure) + \
                    URL_INCLUDE
    elif isinstance(structure, str):
        query_url = URL_PREFIX + \
                    "[{0}$eq'{1}']".format(entry_type, structure) + \
                    URL_INCLUDE
    else:
        raise ValueError('{} is an invalid structure identifier'
                         .format(structure))

    print('accessing {}...'.format(query_url))

    # make the query
    r = requests.get(query_url)
    root = ET.fromstring(r.content)

    return (True, root) if root.text else False


def get_structure_info(structure, entry_type='id', attributes='all'):
    """
    get structure attributes according to the identifier
    (e.g., gene id or acronym)
    :param structure: str or int
        structure identifier, either an int (id, sphinx-id),
        or a str (acronym, name, safe-name, etc.), specified
        in entry_type.
        recommend to use structure ID (int) or
        acronym (str, case sensitive).
    :param entry_type: str, optional
        the type of structure identifier. default: 'id'
        supported:
            'acronym',
            'id',
            'name',
            'safe-name',
            'sphinx-id'
    :param attributes: str, or list, optional
        the attributes of the structure to return.
        default:'all', returning the values of all structure attributes.
        supported:
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

    :return: int, str or dict
        can be integer (id, hemisphere-id, parent-structure-id,
        etc.), or str (name, acronym, etc.) if only one attribute
        is specified, or dict {attribute:value} if multiple attributes
        are given

    """
    check_structure = check_structure_validity(
        structure, entry_type
    )

    if len(check_structure) == 1:
        raise ValueError('{0} {1} is invalid'
                         .format(entry_type, structure))

    # the object to parse
    root = check_structure[1]

    if attributes == 'all':
        attr_list = STRUCTURE_ATTRIBUTES
    elif isinstance(attributes, list):
        attr_list = attributes

    structure_info = dict()

    if attributes == 'all' or isinstance(attributes, list):
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

    else:  # single attribute is given
        # return the single attr value, or AttributeError
        return _get_single_structure_attribute(root, attributes)

    return structure_info


def _get_single_structure_attribute(root, attr):
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


def get_structure_coordinates(structure, entry_type=id):
    """
    get structure coordinates in reference atlas
    :param structure: str or int
        structure identifier, either an int (id, sphinx-id),
        or a str (acronym, name, safe-name, etc.), specified
        in entry_type.
        recommend to use structure ID (int) or
        acronym (str, case sensitive).
    :param entry_type: str, optional
        the type of structure identifier. default: 'id'
        supported:
            'acronym',
            'id',
            'name',
            'safe-name',
            'sphinx-id'
    :return: coor: list of (3, ) tuples
        left and right hemisphere coordinates [(x, y, z), (x, y, z)]
    """
    check_structure = check_structure_validity(
        structure, entry_type
    )

    if len(check_structure) == 1:
        raise ValueError('{0} {1} is invalid'
                         .format(entry_type, structure))

    # the object to parse
    root = check_structure[1]

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

    return coor
