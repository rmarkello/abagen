import itertools
import numpy as np
import pandas as pd
import pytest
from abagen import allen, correct
from abagen.datasets import fetch_desikan_killiany

ATLAS = fetch_desikan_killiany()


@pytest.fixture(scope='module')
def donor_expression(testdir, testfiles):
    return allen.get_expression_data(ATLAS['image'], ATLAS['info'],
                                     exact=False, return_donors=True,
                                     data_dir=testdir,
                                     donors=['12876', '15496'])


def test_remove_distance(donor_expression):
    expr = pd.concat(donor_expression).groupby('label').aggregate(np.mean)
    coexpr = np.corrcoef(expr)
    for atlas_info in [None, ATLAS['info']]:
        out = correct.remove_distance(coexpr, ATLAS['image'], atlas_info)
        assert np.allclose(out, out.T)
        assert isinstance(out, np.ndarray)

    # subset expression data + and atlas_info
    coexpr = np.corrcoef(expr.iloc[:-1])
    removed_label = pd.read_csv(atlas_info).iloc[:-1]
    out = correct.remove_distance(coexpr, ATLAS['image'], removed_label,
                                  labels=removed_label.id)
    assert np.allclose(out, out.T)
    assert isinstance(out, np.ndarray)
    assert len(out) == len(removed_label)

    with pytest.raises(ValueError):
        correct.remove_distance(np.corrcoef(expr), ATLAS['image'],
                                removed_label, labels=removed_label.id)

    with pytest.raises(ValueError):
        correct.remove_distance(expr, ATLAS['image'], ATLAS['info'])


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
