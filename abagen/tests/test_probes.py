# -*- coding: utf-8 -*-
"""
Tests for abagen.probes module
"""

import numpy as np
import pandas as pd
import pytest

import abagen
from abagen import probes


def test_reannotate_probes(testfiles):
    # set up a few useful variables
    probe_file = testfiles['probes'][0]

    # should work with either a filename _or_ a dataframe
    reannot = probes.reannotate_probes(probe_file)
    probe_df = abagen.io.read_probes(probe_file)
    pd.testing.assert_frame_equal(reannot, probes.reannotate_probes(probe_df))

    # expected output
    cols = ['probe_name', 'gene_symbol', 'entrez_id']
    assert np.all(reannot.columns == cols)
    assert reannot.index.name == 'probe_id'
    assert reannot.shape == (45892, 3)


def test_filter_probes(testfiles):
    # set up a few useful variables
    pacall = testfiles['pacall']
    probe_file = testfiles['probes'][0]
    probe_df = abagen.io.read_probes(probe_file)

    # should work with either a filename _or_ a dataframe
    filtered = probes.filter_probes(pacall, probe_file)
    pd.testing.assert_frame_equal(
        filtered,
        probes.filter_probes(pacall, probe_df)
    )

    # expected output with default threshold
    cols = [
        'probe_name', 'gene_id', 'gene_symbol',
        'gene_name', 'entrez_id', 'chromosome'
    ]
    assert np.all(filtered.columns == cols)
    assert filtered.index.name == 'probe_id'
    assert filtered.shape == (38176, 6)

    # threshold of zero should just return the full probe dataframe
    no_filter = probes.filter_probes(pacall, probe_file, threshold=0.0)
    assert len(no_filter) == len(probe_df)

    # threshold of one should NOT return an empty dataframe (that's useless)
    # this will return all the probes that are greater than background noise
    # across ALL samples from ALL provided donors
    max_filter = probes.filter_probes(pacall, probe_file, threshold=1.0)
    assert len(max_filter) == 11878

    # threshold is clipped to the (0, 1) range, so ensure equivalence w/above
    below = probes.filter_probes(pacall, probe_file, threshold=-1.0)
    pd.testing.assert_frame_equal(no_filter, below)
    above = probes.filter_probes(pacall, probe_file, threshold=2.0)
    pd.testing.assert_frame_equal(max_filter, above)


@pytest.mark.xfail(run=False)
def test_groupby_and_apply():
    assert False


def test_max_idx():
    df = pd.DataFrame(dict(a=[1, 2, 3], b=[6, 5, 4]),
                      index=['one', 'two', 'three'])

    # should return max index of column specified
    assert probes._max_idx(df, column='a') == 'three'
    assert probes._max_idx(df, column='b') == 'one'

    # not specifying column should use first numerical column (i.e., 'a')
    assert probes._max_idx(df) == 'three'

    # what if we only have a single column dataframe?
    dfa = df[['a']]
    assert probes._max_idx(dfa, column='a') == 'three'
    assert probes._max_idx(dfa) == 'three'

    # string column anywhere in dataframe is bad
    with pytest.raises(TypeError):
        df['b'] = df['b'].astype(str)
        probes._max_idx(df)


def test_max_loading():
    def _get_ellipse_coords(ang, a=2, b=0.5, origin=(0, 0), rotate=0):
        # convert to radians
        ang = ang * (np.pi / 180)
        rotate = rotate * (np.pi / 180)

        # calculate new points
        x = (a * np.cos(ang) * np.cos(rotate)) \
            - (b * np.sin(ang) * np.sin(rotate)) \
            + origin[0]
        y = (a * np.cos(ang) * np.sin(rotate)) \
            + (b * np.sin(ang) * np.cos(rotate)) \
            + origin[1]

        return (x, y)

    # grab a few points from a rotated ellipse and make a dataframe
    angles = np.arange(0, 360, 45)
    ellipse = np.row_stack([_get_ellipse_coords(f, rotate=45) for f in angles])
    df = pd.DataFrame(ellipse)

    assert probes._max_loading(df) == 0


def test_correlate():
    df = pd.DataFrame(dict(a=[1.50, 1, 3],   # maximum correlation
                           b=[2.00, 2, 3],   # maximum variance
                           c=[2.75, 3, 4]),  # maximum intensity
                      index=['one', 'two', 'three'])

    # when there are >2 rows the provided method is ignored
    assert probes._correlate(df, method='variance') == 'one'
    assert probes._correlate(df, method='intensity') == 'one'

    # when there are exactly two rows we use the specified method
    df = df.iloc[1:]
    assert probes._correlate(df, method='variance') == 'two'
    assert probes._correlate(df, method='intensity') == 'three'

    # and one row just returns the index of that row
    df = df.iloc[1:]
    assert probes._correlate(df, method='variance') == 'three'
    assert probes._correlate(df, method='intensity') == 'three'


@pytest.mark.xfail(run=False)
def test_diff_stability(testfiles):
    assert False


@pytest.mark.xfail(run=False)
def test_average():
    assert False


@pytest.mark.parametrize('method', [
    'average',
    'max_intensity',
    'max_variance',
    'pc_loading',
    'corr_intensity',
    'corr_variance',
    'diff_stability'
])
@pytest.mark.xfail(run=False)
def test_collapse_probes(testfiles, method):
    assert False
