# -*- coding: utf-8 -*-
"""
Tests for abagen.utils module
"""

import numpy as np
import pandas as pd
import pytest

from abagen import utils


@pytest.mark.xfail
def test_check_img():
    assert False


def test_check_atlas_info(atlas):
    # general usage (providing two filenames) works as expected
    out = utils.check_atlas_info(atlas['image'], atlas['info'])
    assert all(out.columns == ['label', 'hemisphere', 'structure'])
    assert out.index.name == 'id'

    # can accept dataframe as input
    atlas_df = pd.read_csv(atlas['info'])
    out2 = utils.check_atlas_info(atlas['image'], atlas_df)
    pd.testing.assert_frame_equal(out, out2)

    # setting ID as index of dataframe is acceptable usage
    atlas_df = atlas_df.set_index('id')
    out3 = utils.check_atlas_info(atlas['image'], atlas_df)
    pd.testing.assert_frame_equal(out, out3)

    # check that coercion of different hemisphere designations works
    atlas_df.loc[atlas_df['hemisphere'] == "L", 'hemisphere'] = "lh"
    atlas_df.loc[atlas_df['hemisphere'] == "R", 'hemisphere'] = "r"
    out4 = utils.check_atlas_info(atlas['image'], atlas_df)
    pd.testing.assert_frame_equal(out, out4)

    # validate = True returns None
    none = utils.check_atlas_info(atlas['image'], atlas['info'], validate=True)
    assert none is None

    # providing labels allows for missing ids in atlas_info (i.e., does not
    # raise ValueError)
    drop_last_df = atlas_df.copy().iloc[:-1]
    out5 = utils.check_atlas_info(atlas['image'], drop_last_df,
                                  labels=range(1, 83))
    assert len(out5) == 82

    # not a filename or dataframe = failure
    with pytest.raises(TypeError):
        utils.check_atlas_info(atlas['image'], [1, 2, 3])

    # missing data = failure
    empty_df = pd.DataFrame(columns=['id', 'hemisphere', 'structure'])
    with pytest.raises(ValueError):
        utils.check_atlas_info(atlas['image'], empty_df)

    # invalid hemisphere designations
    bad_hemi_df = atlas_df.copy()
    bad_hemi_df.loc[1, 'hemisphere'] = 'notahemisphere'
    with pytest.raises(ValueError):
        utils.check_atlas_info(atlas['image'], bad_hemi_df)

    # invalid structural designation
    bad_struct_df = atlas_df.copy()
    bad_struct_df.loc[1, 'structure'] = 'notastructure'
    with pytest.raises(ValueError):
        utils.check_atlas_info(atlas['image'], bad_struct_df)


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


@pytest.mark.xfail
def test_get_unique_labels():
    assert False


@pytest.mark.xfail
def test_get_centroids():
    assert False


@pytest.mark.xfail
def test_closest_centroid():
    assert False


@pytest.mark.parametrize('ijk, xyz', [
    ([0, 0, 0],
     np.array([[-90, -150, -80]])),
    ([[10, 10, 10], [100, 50, 100]],
     np.array([[-80, -140, -70], [10, -100, 20]])),
    ([[54, 32, 20], [82, 205, 38], [32, 51, 82]],
     np.array([[-36, -118, -60], [-8, 55, -42], [-58, -99, 2]])),
])
def test_coords_transform(ijk, xyz):
    affine = np.array([[1, 0, 0, -90],
                       [0, 1, 0, -150],
                       [0, 0, 1, -80],
                       [0, 0, 0, 1]])

    assert np.all(utils.xyz_to_ijk(xyz, affine) == ijk)
    assert np.all(utils.ijk_to_xyz(ijk, affine) == xyz)

    with pytest.raises(ValueError):
        utils.xyz_to_ijk([[10, 10], [20, 30]], affine)
    with pytest.raises(ValueError):
        utils.ijk_to_xyz([[10, 10], [20, 30]], affine)


@pytest.mark.xfail
def test_expand_roi():
    assert False
