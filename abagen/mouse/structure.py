# -*- coding: utf-8 -*-
"""
Functions to make mouse structure queries and manipulations
"""

import json
import requests

# available attributes of structure query
_STRUCTURE_ATTRIBUTES = [
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


def _make_structure_query(structure_id=None, structure_acronym=None,
                          structure_name=None, use_developing=False,
                          attributes=None, suffix=None, verbose=False):
    """
    Check if provided structure is valid and has records in the Allen database

    Parameters
    ----------
    structure_id : int, optional
        Numerical structure ID
    structure_acronym : str, optional
        Short-form structure acronym (case sensitive)
    structure_name : str, optional
        Full structure name (case sensitive)
    verbose : bool, optional
        Whether to print query URL. Default: False

    Returns
    -------
    response : list
        Response from Allen API; empty list if no response

    Raises
    ------
    TypeError
        When none of the required parameters are specified

    Examples
    --------
    Check if structure ID 1018 is valid
    >>> validity, _, _ = check_structure_validity(structure_id=1018)
    >>> validity
    True

    Check if structure acronym SSp is valid
    >>> validity, _, _ = check_structure_validity(structure_acronym='SSp')
    >>> Validity
    True
    """

    url = 'http://api.brain-map.org/api/v2/data/Structure/query.json?' \
          'criteria=[{}${}{}],ontology[id$'
    url += 'in1,12]' if use_developing else 'eq1]'

    if attributes is not None:
        url += '&only={}'.format(','.join(attributes))
    if suffix is not None:
        url += suffix

    # preferred order of inputs: id > acronym > name
    query, params = 'eq', [structure_id, structure_acronym, structure_name]
    for criteria, param in zip(['id', 'acronym', 'name'], params):
        if param is not None:
            if isinstance(param, list):
                query, param = 'in', ','.join([str(f) for f in param])
            break
    if param is None:
        params = ['structure_id', 'structure_acronym', 'structure_name']
        raise TypeError('At least one of {} must be specified.'.format(params))

    # formt URL with desired information and make query to API
    query_url = url.format(criteria, query, param)
    if verbose:
        print("Querying {}...".format(query_url))
    root = json.loads(requests.get(query_url).content)

    return root.get('msg', [])


def available_structure_info():
    """ Lists available attributes for :func:`abagen.mouse.get_structure_info`
    """

    return _STRUCTURE_ATTRIBUTES.copy()


def get_structure_info(structure_id=None, structure_acronym=None,
                       structure_name=None, attributes=None, verbose=False):
    """
    get attributes of a structure

    Parameters
    ----------
    structure_id : int, optional
        Numerical structure ID
    structure_acronym : str, optional
        Short-form structure acronym (case sensitive)
    structure_name : str, optional
        Full structure name (case sensitive)
    attributes : str or list, optional
        Which attributes / information to obtain for the provided structure.
        See :func:`abagen.mouse.available_structure_info` for list of available
        attributes to request. If not specified all available attributes will
        be returned. Default: None
    verbose : bool, optional
        Whether to print status messages. Default: False

    Returns
    -------
    info : int, str, or dict
        If `attributes` is a str, returns an int or str depending on specified
        attribute. If `attributes` is a list, return a dict where keys are
        attributes and values are str or int.

    Raises
    ------
    ValueError
        the structure given is invalid

    Examples
    --------
    # get the full name of structure 1018
    >>> structure_name = get_structure_info(
    >>>     structure_id=1018, attributes='name'
    >>> )
    >>> structure_name
    ['Ventral auditory area']

    """

    # determine which attributes to request; if we don't have to request all
    # of them then we can speed up the API call
    if attributes is None:
        attributes = _STRUCTURE_ATTRIBUTES
    elif isinstance(attributes, str):
        attributes = [attributes]
    attributes = [attr.replace('-', '_') for attr in attributes]

    info = _make_structure_query(structure_id, structure_acronym,
                                 structure_name, attributes=attributes,
                                 verbose=verbose)
    if len(info) == 0:
        provided = [structure_id, structure_acronym, structure_name]
        raise ValueError('Structure {} is invalid.'
                         .format(*[i for i in provided if i is not None]))

    if len(attributes) == 1:
        info = [f[attributes[0]] for f in info]
        if len(info) == 1:
            return info[0]

    return info


def get_structure_coordinates(structure_id=None, structure_acronym=None,
                              structure_name=None, verbose=False):
    """
    get structure coordinates in reference atlas.

    Parameters
    ----------
    structure_id : int, optional
        structure ID
    structure_acronym : str, optional
        structure acronym (case sensitive)
    structure_name : str, optional
        structure name (case sensitive)

    Returns
    -------
    coor : dict
        {reference-space-id:[(x, y, z)]}

    Raises
    ------
    ValueError
        structure given is invalid

    Examples
    --------
    Get the coordinates of structure ID 1018
    >>> coor = get_structure_coordinates(structure_id=1018)
    >>> coor
    {10: [(7800, 3400, 1050)]}
    f the structure has records in mouse brain atlas
    >>> # coordinates in space with reference id 10
    >>> # (Mouse brain atlas left hemisphere)
    >>> coor[10]
    [(7800, 3400, 1050)]
    """

    suffix = ('&include={0}&only={0},{0}.structure_id,{0}.reference_space_id,'
              '{0}.x,{0}.y,{0}.z'.format('structure_centers'))

    info = _make_structure_query(structure_id, structure_acronym,
                                 structure_name, suffix=suffix,
                                 verbose=verbose)
    if len(info) == 0:
        provided = [structure_id, structure_acronym, structure_name]
        raise ValueError('Structure {} is invalid.'
                         .format(*[i for i in provided if i is not None]))

    coords = dict()
    for struct in info:
        for item in struct['structure_centers']:
            coords.setdefault(item['structure_id'], dict())
            reference_atlas_id = int(item['reference_space_id'])
            xyz = (int(item.get('x')), int(item.get('y')), int(item.get('z')))

            if reference_atlas_id in coords:
                coords[item['structure_id']][reference_atlas_id].append(xyz)
            else:
                coords[item['structure_id']][reference_atlas_id] = [xyz]

    # if only one structure was queried don't return a nasty nested dictionary
    if len(coords) == 1:
        return list(coords.values())[0]

    return coords
