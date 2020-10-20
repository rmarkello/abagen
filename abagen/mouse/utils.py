# -*- coding: utf-8 -*-

import json
from typing import Iterable
import urllib
from urllib.request import HTTPError, urlopen


def _coerce_inputs(id=None, acronym=None, name=None):
    """
    Checks that at least one of required inputs was supplied

    Parameters
    ----------
    id : int, optional
        Numerical ID
    acronym : str, optional
        Short-form acronym (case sensitive)
    name : str, optional
        Full name (case sensitive)

    Returns
    -------
    criteria : str
        Should be provided to `criteria` kwarg of :func:`_make_api_query`
    """

    # preferred order of inputs: id > acronym > name
    query, params = 'eq', [id, acronym, name]
    for criteria, param in zip(['id', 'acronym', 'name'], params):
        if param is not None:
            if isinstance(param, Iterable) and not isinstance(param, str):
                # if we're doing "in" queries with strings they must be quoted
                if criteria in ['acronym', 'name']:
                    param = ["'{}'".format(f) for f in param]
                query, param = 'in', ','.join([str(f) for f in param])
            break
    if param is None:
        params = ['id', 'acronym', 'name']
        raise TypeError('At least one of {} must be specified.'.format(params))

    return '[{}${}{}]'.format(criteria, query, param)


def _make_api_query(dtype, includes=None, criteria=None, attributes=None,
                    suffix=None, returns='msg', verbose=False):
    """
    """

    url = 'https://api.brain-map.org/api/v2/data/{}/query.json?'.format(dtype)

    params = [includes, criteria, attributes]
    for key, value in zip(['include', 'criteria', 'only'], params):
        if value is not None:
            if isinstance(value, list):
                value = ','.join(value)
            url += '{}={}&'.format(key, urllib.parse.quote_plus(value))

    if suffix is not None:
        url += suffix

    if verbose:
        print("Querying {}...".format(urllib.parse.unquote_plus(url)))
    response = urlopen(url)
    if response.status != 200:
        raise HTTPError('Failed to query API with code {}: {}'
                        .format(response.status, response.reason))

    info = json.loads(response.read().decode('utf-8'))

    if not info['success']:
        raise ValueError('Provided query {} is invalid. Please check '
                         'parameters and try again.'
                         .format(urllib.parse.unquote_plus(url)))
    elif info['total_rows'] == 0:
        raise ValueError('Provided query {} returned no results. Please '
                         'check parameters and try again.'
                         .format(urllib.parse.unquote_plus(url)))

    if returns is not None:
        info = info.get(returns, [])

    return info
