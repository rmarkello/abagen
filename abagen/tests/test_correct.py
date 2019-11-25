# -*- coding: utf-8 -*-
"""
Tests for abagen.correct module
"""

import itertools

import numpy as np
import pandas as pd
import pytest
import scipy.stats as sstats

from abagen import allen, correct, io


@pytest.fixture(scope='module')
def donor_expression(testfiles, atlas):
    return allen.get_expression_data(atlas['image'], atlas['info'],
                                     exact=False, return_donors=True,
                                     donors=['12876', '15496'])


def test__unpack_tuple():
    assert correct._unpack_tuple((3,)) == 3
    assert correct._unpack_tuple((3, 3)) == (3, 3)
    assert correct._unpack_tuple([2]) == 2
    assert correct._unpack_tuple([2, 4]) == [2, 4]
    assert correct._unpack_tuple(np.array([3])) == 3
    assert np.all(correct._unpack_tuple(np.array([3, 3])) == [3, 3])


def test__batch():
    rs = np.random.RandomState(1234)
    # p-values for ANOVA should all be ~0 (large group differences) before
    # batch correction
    y = [rs.normal(size=(100, 1000)) + f for f in [5, 0, 0]]
    assert np.allclose(sstats.f_oneway(*y)[1], 0)

    # F-values for ANOVA should all be ~0 (no group differences) after batch
    # correction; p-values returned here are sometimes NaN so not a good test
    out = correct._batch_correct(y)
    assert np.allclose(sstats.f_oneway(*out)[0], 0)

    # mean expressions after correction should be ~equal
    assert np.allclose([o.mean() for o in out], 1.24871965683026)

    with pytest.raises(ValueError):
        correct._batch_correct([y[0]])


def test__rescale():
    rs = np.random.RandomState(1234)
    y = rs.normal(size=(100, 1000)) + 10
    out = correct._rescale(y)

    # default max = 1, min =0
    assert np.allclose(out.max(axis=0), 1) and np.allclose(out.min(axis=0), 0)

    # can specify alternative min/max
    out = correct._rescale(y, low=5, high=6)
    assert np.allclose(out.max(axis=0), 6) and np.allclose(out.min(axis=0), 5)

    # different axis works, too!
    out = correct._rescale(y, axis=1)
    assert np.allclose(out.max(axis=1), 1) and np.allclose(out.min(axis=1), 0)


@pytest.mark.parametrize('a', [0, 1])
def test__rs(a):
    rs = np.random.RandomState(1234)

    # create an array with a pretty ridiculous outlier effect to try and fix
    y = rs.normal(size=(100, 1000))
    y[0] += 1000
    y[:, 0] += 1000
    out = correct._rs(y, axis=a)

    # max will always be less than one, min will always be greater than zero
    assert np.all(out.max(axis=a) <= 1) and np.all(out.min(axis=a) >= 0)

    # we should have reduced skewness / kurtosis compared to the original
    assert np.all(sstats.skew(out, axis=a) < sstats.skew(y, axis=a))
    assert np.all(sstats.kurtosis(out, axis=a) < sstats.kurtosis(y, axis=a))

    # this is a weird test; we're gonna bin the data at 0.2 intervals and make
    # sure no bins are empty. if one is something probably went wrong, right?
    for low in np.arange(0, 1, 0.2):
        hi = low + 0.2 + np.spacing(1)  # include 1
        assert np.all(np.sum(np.logical_and(out >= low, out < hi), axis=a) > 0)


@pytest.mark.parametrize('a', [0, 1])
def test__srs(a):
    rs = np.random.RandomState(1234)

    # create an array with a pretty ridiculous outlier effect to try and fix
    y = rs.normal(size=(100, 1000))
    y[0] += 1000
    y[:, 0] += 1000
    out = correct._srs(y, axis=a)

    # max will always be one, min will always be zero
    assert np.allclose(out.max(axis=a), 1) and np.allclose(out.min(axis=a), 0)

    # we should have reduced skewness / kurtosis compared to the original
    assert np.all(sstats.skew(out, axis=a) < sstats.skew(y, axis=a))
    assert np.all(sstats.kurtosis(out, axis=a) < sstats.kurtosis(y, axis=a))

    # this is a weird test; we're gonna bin the data at 0.2 intervals and make
    # sure no bins are empty. if one is something probably went wrong, right?
    for low in np.arange(0, 1, 0.2):
        hi = low + 0.2 + np.spacing(1)  # include 1
        assert np.all(np.sum(np.logical_and(out >= low, out < hi), axis=a) > 0)


def test_normalize_expression_real(testfiles):
    # load in data and add some NaN values for "realness"
    micro = [io.read_microarray(f).T for f in testfiles['microarray']]
    inds = [[5, 15, 25], [0, 10, 20]]
    for n, idx in enumerate(inds):
        micro[n].iloc[idx] = np.nan

    # min-max scaling
    out = correct.normalize_expression(micro, norm='minmax')
    for exp, idx in zip(out, inds):
        assert np.all(np.isnan(exp.iloc[idx]))
        exp = exp.dropna(axis=1, how='all')
        assert np.allclose(exp.max(axis=0), 1)
        assert np.allclose(exp.min(axis=0), 0)
    del out

    # robust sigmoid
    out = correct.normalize_expression(micro, norm='rs')
    for exp, idx in zip(out, inds):
        assert np.all(np.isnan(exp.iloc[idx]))
        exp = exp.dropna(axis=1, how='all')
        assert np.all(exp.max(axis=0) <= 1)
        assert np.all(exp.min(axis=0) >= 0)
    del out

    # scaled robust sigmoid
    out = correct.normalize_expression(micro, norm='srs')
    for exp, idx in zip(out, inds):
        assert np.all(np.isnan(exp.iloc[idx]))
        exp = exp.dropna(axis=1, how='all')
        assert np.allclose(exp.max(axis=0), 1)
        assert np.allclose(exp.min(axis=0), 0)
    del out

    # centering
    out = correct.normalize_expression(micro, norm='center')
    for exp, idx in zip(out, inds):
        assert np.all(np.isnan(exp.iloc[idx]))
        exp = exp.dropna(axis=1, how='all')
        assert np.allclose(exp.mean(axis=0), 0)
    del out

    # z-scoring: mean = 0, std = 1
    out = correct.normalize_expression(micro, norm='zscore')
    for exp, idx in zip(out, inds):
        assert np.all(np.isnan(exp.iloc[idx]))
        exp = exp.dropna(axis=1, how='all')
        assert np.allclose(exp.mean(axis=0), 0)
        assert np.allclose(exp.std(axis=0, ddof=1), 1)
    del out

    # # batch correct: force means identical
    # out = correct.normalize_expression(micro, norm='batch')
    # assert np.allclose(*[e.mean(axis=0, skipna=True) for e in out])
    # # the NaN values should still be there, though
    # for exp, idx in zip(out, inds):
    #     assert np.all(np.isnan(exp.iloc[idx]))

    # invalid norm parameter
    with pytest.raises(ValueError):
        correct.normalize_expression(micro, norm='notanorm')

    # # can't do batch correction with only one donor
    # with pytest.raises(ValueError):
    #     correct.normalize_expression(micro[0], norm='batch')


def test_remove_distance(donor_expression, atlas):
    expr = pd.concat(donor_expression).groupby('label').aggregate(np.mean)
    expr = expr.dropna(axis=1, how='any')
    coexpr = np.corrcoef(expr)
    for atlas_info in [None, atlas['info']]:
        out = correct.remove_distance(coexpr, atlas['image'], atlas_info)
        assert np.allclose(out, out.T)
        assert isinstance(out, np.ndarray)

    # subset expression data + and atlas_info
    coexpr = np.corrcoef(expr.iloc[:-1])
    removed_label = pd.read_csv(atlas_info).iloc[:-1]
    out = correct.remove_distance(coexpr, atlas['image'], removed_label,
                                  labels=removed_label.id)
    assert np.allclose(out, out.T)
    assert isinstance(out, np.ndarray)
    assert len(out) == len(removed_label)

    with pytest.raises(ValueError):
        correct.remove_distance(np.corrcoef(expr), atlas['image'],
                                removed_label, labels=removed_label.id)

    with pytest.raises(ValueError):
        correct.remove_distance(expr, atlas['image'], atlas['info'])


def test_resid_dist():
    dv = np.array([1, 2, 3, 4, 5])
    # residualizing against self should yield 0
    assert np.allclose(correct._resid_dist(dv, iv=dv), 0)
    # residualizing against perfectly anticorrelated should also yield 0
    assert np.allclose(correct._resid_dist(dv, iv=dv[::-1]), 0)
    # residualizing against scaled self should also yield 0 (intercept incl)
    assert np.allclose(correct._resid_dist(dv, iv=(dv + 10)), 0)
    # residualizing against constant should yield de-meaned input
    assert np.allclose(correct._resid_dist(dv, iv=np.ones_like(dv)),
                       dv - dv.mean())


def test_keep_stable_genes(donor_expression):
    for thr, per, rank in itertools.product(np.arange(0, 1, 0.1),
                                            [True, False],
                                            [True, False]):
        out = correct.keep_stable_genes(donor_expression, threshold=thr,
                                        percentile=per, rank=rank)
        assert all([isinstance(f, pd.DataFrame) for f in out])
        for df1, df2 in itertools.combinations(out, 2):
            assert df1.shape == df2.shape
