# -*- coding: utf-8 -*-
"""
Tests for abagen.process module
"""

import numpy as np
import pandas as pd

from abagen import io, process


def test_normalize_expression_real(testfiles):
    mi = io.read_microarray(testfiles['microarray'][0]).T

    srs = process.normalize_expression(mi, norm='srs')
    assert np.allclose(srs.max(axis=0), 1)
    assert np.allclose(srs.min(axis=0), 0)

    zscore = process.normalize_expression(mi, norm='zscore')
    assert np.allclose(zscore.mean(axis=0), 0)
    assert np.allclose(zscore.std(axis=0), 1)


def test_aggregate_donors_real(testfiles):
    microarray = [io.read_microarray(m) for m in testfiles['microarray']]

    # we need our expression arrays to have the same # of regions (rows)
    microarray = [m.T.iloc[:10] for m in microarray]
    for m in microarray:
        m.index.name = 'label'

    # only two donors so these should be identical (and should both work)
    mean = process.aggregate_donors(microarray, metric='mean')
    median = process.aggregate_donors(microarray, metric='median')
    pd.testing.assert_frame_equal(mean, median)
