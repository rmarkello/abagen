# -*- coding: utf-8 -*-
"""
Tests for abagen.utils module
"""

import numpy as np
import pytest

from abagen import utils


@pytest.mark.xfail
def test_check_dict():
    assert False


@pytest.mark.xfail
def test_flatten_dict():
    assert False


@pytest.mark.xfail
def test_first_entry():
    assert False


@pytest.mark.parametrize('metric, check, confirm, kwargs', [
    ('mean', np.array([[1, 2], [3, 4], [5, 6]]),
     np.array([1.5, 3.5, 5.5]), {'axis': 1}),
    ('median', np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]),
     np.array([2, 5, 8]), {'axis': 1}),
])
def test_check_metric(metric, check, confirm, kwargs):
    metric = utils.check_metric(metric)
    assert np.allclose(metric(check, **kwargs), confirm)

    # not an appropriate str
    with pytest.raises(ValueError):
        utils.check_metric('blargh')

    # doesn't accept axis kwarg
    with pytest.raises(TypeError):
        utils.check_metric(lambda x: np.mean(x))


@pytest.mark.xfail
def test_efficient_corr():
    assert False
