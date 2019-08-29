# -*- coding: utf-8 -*-
"""
Functions for cleaning and processing the AHBA microarray dataset
"""

import numpy as np
import pandas as pd

from . import utils


def normalize_expression(expression):
    """
    Performs scaled robust sigmoid (SRS) normalization on `expression` data

    Parameters
    ----------
    expression : (S, G) pandas.DataFrame
        Microarray expression data, where `S` is samples (or regions) and `G`
        is genes

    Returns
    -------
    normalized : (S, G) pandas.DataFrame
        Data from `expression` normalized separately for each gene
    """

    # get non-NaN values
    data = np.asarray(expression.dropna(axis=0, how='all'))

    # calculate sigmoid normalization
    norm = (data - np.median(data, axis=0)) / np.std(data, axis=0)
    srs = 1 / (1 + np.exp(-norm))

    # rescale normalized values to a unit interval
    srs_min = srs.min(axis=0, keepdims=True)
    srs_max = srs.max(axis=0, keepdims=True)
    scaled = (srs - srs_min) / (srs_max - srs_min)

    # recreate dataframe and fill non-NaN values
    normalized = pd.DataFrame(np.nan,
                              columns=expression.columns,
                              index=expression.index)
    inds = expression[expression.notna().all(axis=1)].index
    normalized.loc[inds] = scaled

    return normalized


def aggregate_donors(expression, metric='mean'):
    """
    Aggregates microarray `expression` across donors using `metric`

    Parameters
    ----------
    expression : list of (R, G) pandas.DataFrame
        Where each entry is the microarray expression of `R` regions across `G`
        genes for a given donor
    metric : {'mean', 'median'} or callable, optional
        Mechanism by which to reduce donor-level expression data into a single
        dataframe. If a callable, should be able to accept an `N`-dimensional
        input and the `axis` keyword argument and return an `N-1`-dimensional
        output. Default: 'mean'

    Returns
    -------
    expression : (R, G) pandas.DataFrame
        Microarray expression for `R` regions in `atlas` for `G` genes,
        aggregated across donors.
    """

    metric = utils.check_metric(metric)

    return pd.concat(expression).groupby('label').aggregate(np.mean)
