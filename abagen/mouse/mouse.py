# -*- coding: utf-8 -*-
"""
Functions to fetch mouse unionization (i.e., expression) data
"""

import itertools

from nibabel.volumeutils import Recoder
import numpy as np
import pandas as pd

from .io import fetch_allenref_structures, fetch_rubinov2015_structures
from .utils import _coerce_inputs, _make_api_query

# available attributes of unionization query
_UNIONIZATION_ATTRIBUTES = [
    'expression_density',
    'expression_energy',
    'sum_expressing_pixel_intensity',
    'sum_expressing_pixels',
    'sum_pixel_intensity',
    'sum_pixels',
    'voxel_energy_cv',
    'voxel_energy_mean'
]


def _get_experiments_from_gene(id=None, acronym=None, name=None,
                               slicing_direction='sagittal', verbose=False):
    """
    Fetches experiment IDs associated with specified gene(s)

    One of `id`, `acronym`, or `name` must be provided.

    Parameters
    ----------
    id : int, optional
        Numerical gene ID
    acronym : str, optional
        Short-form gene acronym (case sensitive)
    name : str, optional
        Full gene name (case sensitive)
    slicing_direction : {'sagittal', 'coronal'}, optional
        Slicing direction of brain tissue
    verbose : bool, optional
        Whether to print status messages. Default: False

    Returns
    -------
    experiments : list of int
        List of experiment IDs that probed the specified gene(s)
    """

    directions = ['sagittal', 'coronal']
    if slicing_direction not in directions:
        raise ValueError('Slicing_direction {} is invalid. Must be in {}.'
                         .format(slicing_direction, directions))

    criteria = [
        '[failed$eqfalse]',
        'products[id$eq1]',
        'genes{}'.format(_coerce_inputs(id=id, acronym=acronym, name=name)),
        'plane_of_section[name$eq{}]'.format(slicing_direction)
    ]

    info = _make_api_query('SectionDataSet', criteria=criteria,
                           attributes='data_sets.id', verbose=verbose)

    return [exp['id'] for exp in info]


def _get_unionization_from_experiment(experiment_id, structures=None,
                                      attributes=None, average=True,
                                      verbose=False):
    """
    Gets unionization data for provided experiment(s) `experiment_id`

    Parameters
    ----------
    experiment_id : int or list
        Numerical experiment ID. If multiple experiments are provided the
        requested `attributes` will be averaged across experiments
    structures : list, optional
        List of structures (id, acronym, or name) for which to get unionization
        information associated with provided `experiment_id`. If not specified
        uses structures documented in [MI1]_. Specifying either the id or name
        is recommended as acronyms are not unique to structures. Default: None
    attributes : str or list, optional
        Which attributes / information to obtain for the provided structure.
        See :func:`abagen.mouse.available_unionization_info` for list of
        available attributes to request. If not specified all available
        attributes will be returned. Default: None
    average : bool, optional
        Whether to average across experiments if `experiment_id` is provided as
        a list. Only experiments probing the same gene will be considered for
        averaging. Default: True
    verbose : bool, optional
        Whether to print status messages. Default: False

    Returns
    -------
    unionization : pandas.DataFrame
        Where columns are unionization attributes and the index corresponds to
        gene ids and strucuture ids
    """

    if isinstance(experiment_id, (str, int)):
        experiment_id = [experiment_id]

    if structures is None:
        # read default structure list (from Rubinov et al., 2015, PNAS)
        structures = fetch_rubinov2015_structures(entry_type='id')
    elif isinstance(structures, (str, int)):
        structures = [structures]

    # we need to coerce all provided structures to be integer ids, NOT strings
    # so fetch all available structures then recode them to ids
    if any(isinstance(f, str) for f in structures):
        structs = np.asarray(fetch_allenref_structures(verbose=False))
        structs = Recoder(structs.tolist(), fields=['acronym', 'id', 'name'])
        structures = list(set(structs.id.get(f) for f in structures))

    # determine which attributes to request; if we don't have to request all
    # of them then we can speed up the API call
    if attributes is None:
        attributes = ['expression_density']
    elif attributes == 'all':
        attributes = _UNIONIZATION_ATTRIBUTES
    elif isinstance(attributes, str):
        attributes = [attributes]

    includes = [
        'structure_unionizes', 'genes'
    ]
    criteria = [
        '[id$in{}]'.format(
            ','.join([str(f) for f in experiment_id])
        ),
        'products[id$eq1]',
        'structure_unionizes[structure_id$in{}]'.format(
            ','.join([str(f) for f in structures])
        ),
    ]
    req_attributes = [
        'id', 'structure_unionizes', 'structure_unionizes.structure_id'
    ] + [
        'structure_unionizes.' + attr for attr in attributes
    ]

    info = _make_api_query('SectionDataSet', includes=includes,
                           criteria=criteria, attributes=req_attributes,
                           verbose=verbose)

    for n, exp in enumerate(info):
        keep = exp['structure_unionizes']
        for struc in keep:
            struc['gene_id'] = exp['genes'][0]['id']
            struc['experiment_id'] = exp['id']
        info[n] = keep

    # construct data frame from requested unionization info
    info = pd.DataFrame(list(itertools.chain.from_iterable(info)))
    if average:
        info = info.groupby(['gene_id', 'structure_id']).mean()
    else:
        info = info.set_index(['gene_id', 'experiment_id', 'structure_id'])

    return info[attributes]


def available_unionization_info():
    """ Lists attributes for :func:`abagen.mouse.get_unionization_from_gene`
    """

    return _UNIONIZATION_ATTRIBUTES


def get_unionization_from_gene(id=None, acronym=None, name=None,
                               slicing_direction='sagittal', structures=None,
                               attributes=None, average=True, verbose=False):
    """
    Gets unionization data for provided gene(s)

    One of `id`, `acronym`, or `name` must be provided.

    Parameters
    ----------
    id : int, optional
        Numerical gene ID
    acronym : str, optional
        Short-form gene acronym (case sensitive)
    name : str, optional
        Full gene name (case sensitive)
    slicing_direction : {'sagittal', 'coronal'}, optional
        Slicing direction of brain tissue
    structures : list, optional
        List of structures (id, acronym, or name) for which to get unionization
        information associated with provided `experiment_id`. If not specified
        uses structures documented in [MI1]_. Specifying either the id or name
        is recommended as acronyms are not unique to structures. Default: None
    attributes : str or list, optional
        Which attributes / information to obtain for the provided gene. See
        :func:`abagen.mouse.available_gene_info` for list of available
        attributes to request. If not specified then only 'expression_density'
        will be returned. Specifying 'all' will return all information.
        Default: None
    average : bool, optional
        Whether to average across experiments if there are multiple experiments
        corresponding to any provided gene(s). Only experiments probing the
        same gene will be considered for averaging, and distinct structures
        will be retained. Default: True
    verbose : bool, optional
        Whether to print status messages. Default: False

    Returns
    -------
    unionization : pandas.DataFrame
        Where columns are unionization attributes and the index corresponds to
        strucuture and gene ids (if `experiments` is provided as a list
        with multiple genes). If `average=False`, `experiments` will also be a
        level in index

    Examples
    --------
    >>> from abagen import mouse
    >>> mouse.get_unionization_from_gene(acronym='Pdyn',
    ...                                  structures=[22, 31])  # doctest: +NORMALIZE_WHITESPACE
                          expression_density
    gene_id structure_id
    18376   22                      0.024840
            31                      0.017199
    >>> mouse.get_unionization_from_gene(acronym=['Ace', 'Cd99'],
    ...                                  structures=[22, 31])  # doctest: +NORMALIZE_WHITESPACE
                          expression_density
    gene_id structure_id
    11210   22                      0.001283
            31                      0.001427
    163028  22                      0.067537
            31                      0.056442
    """  # noqa

    directions = ['sagittal', 'coronal']
    if slicing_direction not in directions:
        raise ValueError('Slicing_direction {} is invalid. Must be in {}.'
                         .format(slicing_direction, directions))

    if structures is None:
        # read default structure list (from Rubinov et al., 2015, PNAS)
        structures = fetch_rubinov2015_structures(entry_type='id')
    elif isinstance(structures, (str, int)):
        structures = [structures]

    # we need to coerce all provided structures to be integer ids, NOT strings
    # so fetch all available structures then recode them to ids
    if any(isinstance(f, str) for f in structures):
        structs = np.asarray(fetch_allenref_structures(verbose=False))
        structs = Recoder(structs.tolist(), fields=['acronym', 'id', 'name'])
        structures = list(set(structs.id.get(f) for f in structures))

    exp_ids = _get_experiments_from_gene(id=id, acronym=acronym, name=name,
                                         slicing_direction=slicing_direction,
                                         verbose=verbose)

    data = _get_unionization_from_experiment(exp_ids, structures=structures,
                                             attributes=attributes,
                                             average=average, verbose=verbose)

    return data
