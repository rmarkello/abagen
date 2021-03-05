# -*- coding: utf-8 -*-
"""
Functions to make mouse structure queries and manipulations
"""

import itertools
import re
import pandas as pd
from .utils import _coerce_inputs, _make_api_query

# available attributes of structure query
_STRUCTURE_ATTRIBUTES = [
    'acronym',
    'atlas_id',
    'color_hex_triplet',
    'depth',
    'graph_id',
    'graph_order',
    'hemisphere_id',
    'id',
    'name',
    'neuro_name_structure_id',
    'neuro_name_structure_id_path',
    'ontology_id',
    'parent_structure_id',
    'safe_name',
    'sphinx_id',
    'structure_id_path',
    'weight'
]


def available_structure_info():
    """ Lists available attributes for :func:`abagen.mouse.get_structure_info`
    """

    return _STRUCTURE_ATTRIBUTES.copy()


def get_structure_info(id=None, acronym=None, name=None, attributes=None,
                       verbose=False):
    """
    Queries Allen API for information about given gene

    One of `structure_id`, `structure_acronym`, or `structure_name` must be
    provided.

    Parameters
    ----------
    id : int, optional
        Numerical structure ID
    acronym : str, optional
        Short-form structure acronym (case sensitive)
    name : str, optional
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
    info : pandas.DataFrame
        Where columns are the requested attributes and index is the provided
        structural identifier type (e.g., 'id', 'acronym', 'name')

    Raises
    ------
    ValueError
        The provided structure is invalid

    Examples
    --------
    Get the full names of structures 22 and 1018:

    >>> from abagen import mouse
    >>> mouse.get_structure_info(id=[22, 1018],
    ...                          attributes=['acronym', 'name'])  # doctest: +NORMALIZE_WHITESPACE
         acronym                                  name
    id
    22      PTLp  Posterior parietal association areas
    1018    AUDv                 Ventral auditory area
    """  # noqa

    criteria = [
        _coerce_inputs(id=id, acronym=acronym, name=name),
        'ontology[id$eq1]'
    ]
    provided = re.search(r'\[(\S+)\$', criteria[0]).group(1)

    # determine which attributes to request; if we don't have to request all
    # of them then we can speed up the API call
    if attributes is None:
        attributes = _STRUCTURE_ATTRIBUTES
    elif isinstance(attributes, str):
        attributes = [attributes]
    attributes = [a for a in attributes if a not in provided]

    for attr in attributes:
        if attr not in _STRUCTURE_ATTRIBUTES:
            raise ValueError('Provided attribute "{}" is invalid; please '
                             'check valid attributes with '
                             'abagen.mouse.available_structure_info().'
                             .format(attr))

    info = _make_api_query('Structure', criteria=criteria,
                           attributes=attributes + [provided], verbose=verbose)
    info = pd.DataFrame(info).set_index(provided)[attributes]

    return info.sort_index()


def get_structure_coordinates(id=None, acronym=None, name=None,
                              reference_space='sagittal', verbose=False):
    """
    Finds xyz coordinates of provided structure(s) in `reference_space`

    Parameters
    ----------
    id : int, optional
        Numerical structure ID
    acronym : str, optional
        Short-form structure acronym (case sensitive)
    name : str, optional
        Full structure name (case sensitive)
    reference_space : {'sagittal', 'coronal'}, optional
        Reference space from which to extract coordinates. Default: 'sagittal'
    verbose : bool, optional
        Whether to print status messages. Default: False

    Returns
    -------
    coords : pandas.DataFrame
        With columns ['structure_id', 'x', 'y', 'z']

    Examples
    --------
    Get the coordinates of structure 1018:

    >>> from abagen import mouse
    >>> mouse.get_structure_coordinates(id=1018)  # doctest: +NORMALIZE_WHITESPACE
       structure_id     x     y     z
    0          1018  7800  3400  1050
    """  # noqa

    spaces = {'sagittal': 10, 'coronal': 9}
    if reference_space not in spaces:
        raise ValueError('Reference_space {} is invalid. Must be in {}.'
                         .format(reference_space, list(spaces.keys())))
    space = spaces[reference_space]

    includes = 'structure_centers'
    criteria = [
        _coerce_inputs(id=id, acronym=acronym, name=name),
        'ontology[id$eq1]',
        'structure_centers[reference_space_id$eq{}]'.format(space)
    ]
    attributes = [
        'structure_centers', 'structure_centers.structure_id',
        'structure_centers.x', 'structure_centers.y', 'structure_centers.z'
    ]

    info = _make_api_query('Structure', includes=includes, criteria=criteria,
                           attributes=attributes, verbose=verbose)
    info = [struct['structure_centers'] for struct in info]
    coords = pd.DataFrame(list(itertools.chain.from_iterable(info)))

    return coords
