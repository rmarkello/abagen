# -*- coding: utf-8 -*-
"""
Tests for abagen.utils module
"""

import nibabel as nib
import numpy as np
import pytest

from abagen import utils


@pytest.mark.parametrize('indict, outdict', [
    ('val', {0: 'val'}),
    (['val1', 'val2'], {0: 'val1', 1: 'val2'}),
    ({'key': 'val'}, {'key': 'val'})
])
def test_check_dict(indict, outdict):
    assert utils.check_dict(indict) == outdict


@pytest.mark.parametrize('indict, subkey, outdict', [
    ({'0': {'key': 0}, '1': {'key': 1}}, 'key', {'0': 0, '1': 1}),
    ({'0': {'key': 0}, '1': {}}, 'key', {'0': 0, '1': None})
])
def test_flatten_dict(indict, subkey, outdict):
    assert utils.flatten_dict(indict, subkey) == outdict


@pytest.mark.parametrize('indict, subkey, out', [
    ({'0': {'key': 'val0'}, '1': {'key': 'val1'}}, None, {'key': 'val0'}),
    ({'0': {'key': 'val0'}, '1': {'key': 'val1'}}, 'key', 'val0')
])
def test_first_entry(indict, subkey, out):
    assert utils.first_entry(indict, subkey) == out


@pytest.mark.parametrize('metric, check, confirm, kwargs', [
    ('mean', np.array([[1, 2], [3, 4], [5, 6]]),
     np.array([1.5, 3.5, 5.5]), {'axis': 1}),
    ('median', np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]),
     np.array([2, 5, 8]), {'axis': 1}),
])
def test_check_metric(metric, check, confirm, kwargs):
    metric = utils.check_metric(metric)
    assert np.allclose(metric(check, **kwargs), confirm)


def test_check_metric_errors():
    # not an appropriate str
    with pytest.raises(ValueError):
        utils.check_metric('blargh')

    # doesn't accept axis kwarg
    with pytest.raises(TypeError):
        utils.check_metric(lambda x: np.mean(x))


def test_efficient_corr():
    # valid inputs
    a, b = np.random.rand(2, 100, 10)
    corrs = utils.efficient_corr(a, b)
    assert len(corrs) == 10

    # known output
    a, b = np.arange(9).reshape(3, 3), np.arange(9).reshape(3, 3)[::-1]
    corrs = utils.efficient_corr(a, b)
    assert np.all(corrs == [-1, -1, -1])

    # empty input yields NaN
    assert np.isnan(utils.efficient_corr([], []))

    # different lengths
    with pytest.raises(ValueError):
        utils.efficient_corr(a[:2], b)


def test_load_gifti(atlas, surface):
    gii = utils.load_gifti(surface['image'][0])
    assert isinstance(gii, nib.GiftiImage)

    # providing a GiftiImage object will just pass that object back
    gii2 = utils.load_gifti(gii)
    assert gii is gii2

    # cannot load non-GiftiImage object
    with pytest.raises(TypeError):
        utils.load_gifti(nib.load(atlas['image']))
