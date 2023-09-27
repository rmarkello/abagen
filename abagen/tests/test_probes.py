# -*- coding: utf-8 -*-
"""
Tests for abagen.probes module
"""

import numpy as np
import pandas as pd
import pytest

import abagen
from abagen import probes_
from abagen.utils import first_entry, flatten_dict
from functools import partial


def test_reannotate_probes(testfiles):
    # set up a few useful variables
    probe_file = first_entry(testfiles, 'probes')

    # should work with either a filename _or_ a dataframe
    reannot = probes_.reannotate_probes(probe_file)
    probe_df = abagen.io.read_probes(probe_file)
    pd.testing.assert_frame_equal(reannot, probes_.reannotate_probes(probe_df))

    # expected output
    cols = ['probe_name', 'gene_symbol', 'entrez_id']
    assert np.all(reannot.columns == cols)
    assert reannot.index.name == 'probe_id'
    assert reannot.shape == (45821, 3)


@pytest.mark.parametrize('threshold, expected_length', [
    # threshold of zero should just return the full probe dataframe
    (0.0, 58692),
    # default threshold
    (0.5, 38176),
    # threshold of one should NOT return an empty dataframe (that's useless)
    # this will return all the probes that are greater than background noise
    # across ALL samples from ALL provided donors
    (1.0, 11878),
    # threshold is clipped to [0, 1] so these should be identical to [0, 1]
    (-1.0, 58692),
    (2.0, 11878),
])
def test_filter_probes(testfiles, threshold, expected_length):
    # set up a few useful variables
    pacall = flatten_dict(testfiles, 'pacall')
    probe_file = first_entry(testfiles, 'probes')
    samples = flatten_dict(testfiles, 'annotation')
    probe_df = abagen.io.read_probes(probe_file)

    # should work with either a filename _or_ a dataframe
    filtered = probes_.filter_probes(pacall, samples, probe_file,
                                     threshold=threshold)
    pd.testing.assert_frame_equal(
        filtered,
        probes_.filter_probes(pacall, samples, probe_df, threshold=threshold)
    )

    # provided threshold returns expected output
    cols = [
        'probe_name', 'gene_id', 'gene_symbol', 'gene_name', 'entrez_id',
        'chromosome'
    ]
    assert np.all(filtered.columns == cols)
    assert filtered.index.name == 'probe_id'
    assert len(filtered) == expected_length


@pytest.mark.parametrize('func, expected', [
    (partial(probes_._max_idx, column="a"), [3000, 5000, 6000]),
    (lambda x: x.index[0], [1000, 4000, 6000])
])
def test_groupby_and_apply(func, expected):
    # lots of set up for this little function...
    index = pd.Series([1000, 2000, 3000, 4000, 5000, 6000], name='probe_id')
    gene_symbol = ['a', 'a', 'a', 'b', 'b', 'c']
    info = pd.DataFrame(dict(gene_symbol=gene_symbol,
                             a=range(6),
                             b=range(12, 6, -1)),
                        index=index)
    prb = pd.DataFrame(dict(gene_symbol=gene_symbol), index=index)
    exp = [info[['a', 'b']].copy()]

    mi = list(probes_._groupby_and_apply(exp, prb, info, func).values())[0]
    pd.testing.assert_frame_equal(mi, info.loc[expected, ['a', 'b']].T)


def test_max_idx():
    df = pd.DataFrame(dict(a=[1, 2, 3],
                           b=[6, 5, 4]),
                      index=['one', 'two', 'three'])

    # should return max index of column specified
    assert probes_._max_idx(df, column='a') == 'three'
    assert probes_._max_idx(df, column='b') == 'one'

    # not specifying column should raise error
    with pytest.raises(TypeError):
        probes_._max_idx(df)

    # what if we only have a single column dataframe?
    dfa = df[['a']]
    assert probes_._max_idx(dfa, column='a') == 'three'
    # assert probes_._max_idx(dfa) == 'three'

    # string column anywhere in dataframe is bad news bears
    # with pytest.raises(TypeError):
    #     df['b'] = df['b'].astype(str)
    #     probes_._max_idx(df)


@pytest.mark.parametrize('angles, index', [
    (np.arange(0, 360, 15), 0),
    (np.arange(180, -180, -15), 12),
])
def test_max_loading(angles, index):
    def _get_ellipse_coords(ang, a=2, b=0.5, origin=(0, 0), rotate=0):
        # convert to radians
        rad = np.pi / 180
        ang, rotate = ang * rad, rotate * rad
        cang, sang = np.cos(ang), np.sin(ang)
        crot, srot = np.cos(rotate), np.sin(rotate)

        # calculate new points
        x = (a * cang * crot) - (b * sang * srot) + origin[0]
        y = (a * cang * srot) + (b * sang * crot) + origin[1]

        return x, y

    # grab a few points from a rotated ellipse and make a dataframe of it
    # the max loading should be wherever angle = 0
    df = pd.DataFrame([_get_ellipse_coords(f, rotate=45) for f in angles])
    assert probes_._max_loading(df) == index


@pytest.mark.parametrize('from_row, method, expected', [
    # when there are >2 rows the provided method is ignored
    (0, 'variance', 'one'),
    (0, 'intensity', 'one'),
    # when there are exactly two rows we use the specified method
    (1, 'variance', 'two'),
    (1, 'intensity', 'three'),
    # and one row just returns the index of that row
    (2, 'variance', 'three'),
    (2, 'variance', 'three'),
])
def test_correlate(from_row, method, expected):
    # a = max correlation, b = max variance, c = max intensity
    df = pd.DataFrame([[1.50, 2.00, 2.75],
                       [1.00, 2.00, 3.00],
                       [3.00, 3.00, 4.00]],
                      index=['one', 'two', 'three'], columns=['a', 'b', 'c'])

    assert probes_._correlate(df.iloc[from_row:], method=method) == expected


def test_diff_stability():
    # so much set up...
    index = pd.Series([1000, 2000, 3000, 4000, 5000, 6000], name='probe_id')
    prb = pd.DataFrame(dict(gene_symbol=['a', 'a', 'a', 'b', 'b', 'c']),
                       index=index)
    expression = [
        pd.DataFrame([[0, 2, 4, 9],
                      [1, 3, 3, 2],
                      [2, 4, 4, 1],
                      [3, 5, 1, 0],
                      [4, 6, 6, 7],
                      [5, 7, 3, 4]],
                     index=index, columns=['one', 'two', 'three', 'four']),
        pd.DataFrame([[1, 2, 3],
                      [3, 2, 1],
                      [1, 3, 2],
                      [1, 3, 2],
                      [1, 2, 3],
                      [1, 2, 3]],
                     index=index, columns=['one', 'two', 'three'])
    ]

    annotation = [
        pd.DataFrame(dict(structure_id=[0, 0, 1, 2]),
                     index=expression[0].columns),
        pd.DataFrame(dict(structure_id=[0, 1, 2]),
                     index=expression[1].columns)
    ]

    expected = [
        pd.DataFrame([[0, 2, 4, 9],
                      [4, 6, 6, 7],
                      [5, 7, 3, 4]],
                     index=pd.Series([1000, 5000, 6000], name='probe_id'),
                     columns=expression[0].columns).T,
        pd.DataFrame([[1, 2, 3],
                      [1, 2, 3],
                      [1, 2, 3]],
                     index=pd.Series([1000, 5000, 6000], name='probe_id'),
                     columns=expression[1].columns).T,
    ]

    out = list(probes_._diff_stability(expression, prb, annotation).values())
    assert len(out) == len(expected)
    for df, exp in zip(out, expected):
        pd.testing.assert_frame_equal(df, exp)


def test_average():
    index = pd.Series([1000, 2000, 3000, 4000, 5000, 6000], name='probe_id')
    prb = pd.DataFrame(dict(gene_symbol=['a', 'a', 'a', 'b', 'b', 'c']),
                       index=index)
    exp = [pd.DataFrame([[0, 12], [1, 11], [2, 10], [3, 9], [4, 8], [5, 7]],
                        index=index, columns=['a', 'b'])]

    # output should be a list of dataframes (only one dataframe in this case)
    out = probes_._average(exp, prb)
    assert len(out) == 1
    out = out[0]

    # confirm output dataframe is as expected
    expected = pd.DataFrame(dict(a=[1.0, 3.5, 5.0],
                                 b=[11.0, 8.5, 7.0]),
                            index=pd.Series(['a', 'b', 'c'],
                                            name='gene_symbol')).T
    pd.testing.assert_frame_equal(out, expected)


@pytest.mark.parametrize('method', [
    'average',
    'mean',
    'max_intensity',
    'max_variance',
    'pc_loading',
    'corr_intensity',
    'corr_variance',
    'diff_stability',
])
def test_collapse_probes(testfiles, method):
    # we've aleady tested the underlying methods so here we just want to do
    # some smoke tests to make sure the function returns what we expected
    # regardless of the provided method
    out = probes_.collapse_probes(flatten_dict(testfiles, 'microarray'),
                                  flatten_dict(testfiles, 'annotation'),
                                  first_entry(testfiles, 'probes'),
                                  method=method)

    out = list(out.values())
    assert len(out) == 2  # number of donors
    assert np.all([len(exp) == n_samp for exp, n_samp in zip(out, [363, 470])])
    assert np.all([len(exp.columns) == 29131 for exp in out])
    assert out[0].index.name == 'sample_id'
    assert out[0].columns.name == 'gene_symbol'


@pytest.mark.parametrize('donor_probes', [
    'aggregate',
    'independent',
    'common'
])
def test_collapse_probes_donors(testfiles, donor_probes):
    out = probes_.collapse_probes(flatten_dict(testfiles, 'microarray'),
                                  flatten_dict(testfiles, 'annotation'),
                                  first_entry(testfiles, 'probes'),
                                  method='max_intensity',
                                  donor_probes=donor_probes)

    out = list(out.values())
    assert len(out) == 2  # number of donors
    assert np.all([len(exp) == n_samp for exp, n_samp in zip(out, [363, 470])])
    assert np.all([len(exp.columns) == 29131 for exp in out])


def test_collapse_probes_errors(testfiles):
    with pytest.raises(ValueError):
        probes_.collapse_probes(flatten_dict(testfiles, 'microarray'),
                                flatten_dict(testfiles, 'annotation'),
                                first_entry(testfiles, 'probes'),
                                method='notamethod')
