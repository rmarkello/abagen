# -*- coding: utf-8 -*-
"""
Tests for abagen.images module
"""

import numpy as np
import pandas as pd
import pytest

from abagen import images


def test_leftify_atlas(atlas):
    out = images.leftify_atlas(atlas['image'])
    assert len(np.unique(out.dataobj)) == 44
    assert np.all(np.asarray(out.dataobj)[98:] == 0)


def test_check_img(atlas):
    # some really basic, silly checks
    out = images.check_img(atlas['image'])
    assert out.header.get_data_dtype() == np.dtype('int32')
    assert len(out.shape) == 3


def test_check_atlas_info(atlas):
    # general usage (providing two filenames) works as expected
    out = images.check_atlas_info(atlas['image'], atlas['info'])
    assert all(out.columns == ['label', 'hemisphere', 'structure'])
    assert out.index.name == 'id'

    # can accept dataframe as input
    atlas_df = pd.read_csv(atlas['info'])
    out2 = images.check_atlas_info(atlas['image'], atlas_df)
    pd.testing.assert_frame_equal(out, out2)

    # setting ID as index of dataframe is acceptable usage
    atlas_df = atlas_df.set_index('id')
    out3 = images.check_atlas_info(atlas['image'], atlas_df)
    pd.testing.assert_frame_equal(out, out3)

    # check that coercion of different hemisphere designations works
    atlas_df.loc[atlas_df['hemisphere'] == "L", 'hemisphere'] = "lh"
    atlas_df.loc[atlas_df['hemisphere'] == "R", 'hemisphere'] = "r"
    out4 = images.check_atlas_info(atlas['image'], atlas_df)
    pd.testing.assert_frame_equal(out, out4)

    # validate = True returns None
    none = images.check_atlas_info(atlas['image'], atlas['info'],
                                   validate=True)
    assert none is None

    # providing labels allows for missing ids in atlas_info (i.e., does not
    # raise ValueError)
    drop_last_df = atlas_df.copy().iloc[:-1]
    out5 = images.check_atlas_info(atlas['image'], drop_last_df,
                                   labels=range(1, 83))
    assert len(out5) == 82

    # not a filename or dataframe = failure
    with pytest.raises(TypeError):
        images.check_atlas_info(atlas['image'], [1, 2, 3])

    # missing data = failure
    empty_df = pd.DataFrame(columns=['id', 'hemisphere', 'structure'])
    with pytest.raises(ValueError):
        images.check_atlas_info(atlas['image'], empty_df)

    # invalid hemisphere designations
    bad_hemi_df = atlas_df.copy()
    bad_hemi_df.loc[1, 'hemisphere'] = 'notahemisphere'
    with pytest.raises(ValueError):
        images.check_atlas_info(atlas['image'], bad_hemi_df)

    # invalid structural designation
    bad_struct_df = atlas_df.copy()
    bad_struct_df.loc[1, 'structure'] = 'notastructure'
    with pytest.raises(ValueError):
        images.check_atlas_info(atlas['image'], bad_struct_df)


@pytest.mark.xfail
def test_get_centroids():
    assert False


@pytest.mark.xfail
def test_closest_centroid():
    assert False
