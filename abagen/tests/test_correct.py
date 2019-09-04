# -*- coding: utf-8 -*-
"""
Tests for abagen.correct module
"""

import itertools

import numpy as np
import pandas as pd
import pytest

from abagen import allen, correct, io


@pytest.fixture(scope='module')
def donor_expression(testdir, testfiles, atlas):
    return allen.get_expression_data(atlas['image'], atlas['info'],
                                     exact=False, return_donors=True,
                                     data_dir=testdir,
                                     donors=['12876', '15496'])


@pytest.mark.xfail
def test__batch():
    assert False


@pytest.mark.xfail
def test__srs():
    assert False


def test_normalize_expression_real(testfiles):
    # load in data and add some NaN values for "realness"
    micro = [io.read_microarray(f).T for f in testfiles['microarray']]
    inds = [[5, 15, 25], [0, 10, 20]]
    for n, idx in enumerate(inds):
        micro[n].iloc[idx] = np.nan

    # min-max scaling (with some extra pizzazz)
    srs = correct.normalize_expression(micro, norm='srs')
    for exp, idx in zip(srs, inds):
        assert np.all(np.isnan(exp.iloc[idx]))
        exp = exp.dropna(axis=1, how='all')
        assert np.allclose(exp.max(axis=0), 1)
        assert np.allclose(exp.min(axis=0), 0)

    # z-scoring: mean = 0, std = 1
    zscore = correct.normalize_expression(micro, norm='zscore')
    for exp, idx in zip(zscore, inds):
        assert np.all(np.isnan(exp.iloc[idx]))
        exp = exp.dropna(axis=1, how='all')
        assert np.allclose(exp.mean(axis=0), 0)
        assert np.allclose(exp.std(axis=0, ddof=1), 1)

    # this is effectively centering so mean = 0
    batch = correct.normalize_expression(micro, norm='batch')
    for exp, idx in zip(batch, inds):
        assert np.all(np.isnan(exp.iloc[idx]))
        exp = exp.dropna(axis=1, how='all')
        assert np.allclose(exp.mean(axis=0), 0)

    # invalid norm parameter
    with pytest.raises(ValueError):
        correct.normalize_expression(micro, norm='notanorm')

    # can't do batch correction with only one donor
    with pytest.raises(ValueError):
        correct.normalize_expression(micro[0], norm='batch')


def test_remove_distance(donor_expression, atlas):
    expr = pd.concat(donor_expression).groupby('label').aggregate(np.mean)
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
