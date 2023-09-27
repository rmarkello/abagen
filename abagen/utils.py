# -*- coding: utf-8 -*-

import gzip

import nibabel as nib
from nibabel.filebasedimages import ImageFileError
import numpy as np
import pandas as pd
from scipy.stats import zscore

AGG_FUNCS = dict(
    mean=np.mean,
    median=np.median
)


def check_dict(dictionary):
    """
    Checks that `dictionary` is a dict and makes it one, if not

    Parameters
    ----------
    dictionary : dict
        Presumptive dictionary

    Returns
    -------
    dictionary : dict
        Actual dictionary. If it wasn't one before, keys are chronological int,
        starting at 0
    """

    if isinstance(dictionary, str):
        dictionary = [dictionary]
    if not isinstance(dictionary, dict):
        return dict(zip(range(len(dictionary)), dictionary))
    return dictionary


def flatten_dict(dictionary, subkey):
    """
    Extracts `subkey` from entries of `dictionary` and creates new dictionary

    Parameters
    ----------
    dictionary : dict
        Input dictionary, where values are sub-dictionaries
    subkey : str
        Key to extract from `dictionary` entries

    Returns
    -------
    flattened : dict
        Flattened input `dictionary`

    Examples
    --------
    >>> test = {'one': {'value': 1}, 'two': {'value': 2}}
    >>> flatten_dict(test, 'value')
    {'one': 1, 'two': 2}
    """

    return {k: v.get(subkey, None) for k, v in check_dict(dictionary).items()}


def first_entry(dictionary, subkey=None):
    """
    Extracts `subkey` from the first entry in `dictionary`

    Parameters
    ----------
     dictionary : dict
        Input dictionary, where values are sub-dictionaries
    subkey : str, optional
        Key to extract from `dictionary` entries. If not provided just returns
        value of first key in `dictionary`. Default: None

    Returns
    -------
    entry
        Extracted entry
    """

    dictionary = check_dict(dictionary)
    entry = dictionary[list(dictionary)[0]]
    if subkey is not None:
        entry = entry.get(subkey, None)

    return entry


def check_metric(metric):
    """
    Confirms `metric` is a valid aggregation metric

    Parameters
    ----------
    metric : {'mean', 'median'} or callable, optional
        Mechanism by which to reduce arrays. If a callable, should be able to
        accept an `N`-dimensional input and the `axis` keyword argument and
        return an `N-1`-dimensional output. Default: 'mean'

    Returns
    -------
    metric : callable
    """

    if not callable(metric):
        try:
            metric = AGG_FUNCS[metric]
        except KeyError:
            raise ValueError('Provided metric {0} is not valid. If supplied'
                             'as string, metric must be in {1}.'
                             .format(metric, list(AGG_FUNCS)))
    else:
        try:
            test = np.random.rand(10, 10, 10)
            assert metric(test, axis=0).shape == (10, 10)
            assert isinstance(metric(test), float)
        except (AssertionError, TypeError):
            raise TypeError('Provided metric {0} does not perform as '
                            'expected. Please ensure it accepts the `axis` '
                            'keyword argument and can reduce an `N`-'
                            'dimensional input to an `N-1`-dimensional output.'
                            .format(metric))

    return metric


def efficient_corr(x, y):
    """
    Computes correlation of matching columns in `x` and `y`

    Parameters
    ----------
    x, y : (N, M) array_like
        Input data arrays

    Returns
    -------
    corr : (M,) numpy.ndarray
        Correlations of columns in `x` and `y`
    """

    x, y = np.asarray(x), np.asarray(y)
    if len(x) != len(y):
        raise ValueError('Provided arrays do not have same length')

    if x.size == 0 or y.size == 0:
        return np.nan

    corr = np.sum(zscore(x, ddof=1) * zscore(y, ddof=1), axis=0) / (len(x) - 1)

    return corr


def labeltable_to_df(labels):
    """
    Converts GIFTI label table dictionaries to uniform pandas.DataFrame

    Parameters
    ----------
    labels : (2,) tuple of dict
        Where dictionaries have parcel IDs as keys and parcel labels as values

    Returns
    -------
    info : pandas.DataFrame or None
        If provided `labels` are not empty, returns a dataframe with info about
        corresponding atlas. Dataframe will have columns ['id', 'label',
        'hemisphere', 'structure'] containing information mapping atlas IDs to
        hemisphere (i.e., "L" or "R") and structure (i.e., "cortex")
    """

    info = pd.DataFrame(columns=['id', 'label', 'hemisphere', 'structure'])
    for table, hemi in zip(labels, ('L', 'R')):
        if len(table) == 0:
            continue
        ids, label = zip(*table.items())
        info = pd.concat([
            info,
            pd.DataFrame(dict(id=ids, label=label, hemisphere=hemi, structure='cortex'))
        ], ignore_index=True)
    info = info.set_index('id').drop([0], axis=0).sort_index()

    if len(info) != 0:
        return info


def load_gifti(img):
    """
    Loads gifti file `img`

    Will try to gunzip `img` if gzip is detected, and will pass pre-loaded
    GiftiImage object

    Parameters
    ----------
    img : os.PathLike or nib.GiftiImage object
        Image to be loaded

    Returns
    -------
    img : nib.GiftiImage
        Loaded GIFTI images
    """

    try:
        img = nib.load(img)
    except (ImageFileError, TypeError) as err:
        # it's gzipped, so read the gzip and pipe it in
        if isinstance(err, ImageFileError) and str(err).endswith('.gii.gz"'):
            with gzip.GzipFile(img) as gz:
                img = nib.GiftiImage.from_bytes(gz.read())
        # it's not a pre-loaded GiftiImage so error out
        elif (isinstance(err, TypeError)
              and not str(err) == 'stat: path should be string, bytes, os.'
                                  'PathLike or integer, not GiftiImage'):
            raise err

    return img
