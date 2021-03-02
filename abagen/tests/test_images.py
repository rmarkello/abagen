# -*- coding: utf-8 -*-
"""
Tests for abagen.images module
"""

import gzip

import nibabel as nib
import numpy as np
import pandas as pd
import pytest

from abagen import images
from abagen.matching import AtlasTree


def test_leftify_atlas(atlas):
    out = images.leftify_atlas(atlas['image'])
    assert len(np.unique(out.dataobj)) == 44
    assert np.all(np.asarray(out.dataobj)[98:] == 0)


@pytest.mark.xfail
def test_relabel_gifti():
    assert False


@pytest.mark.xfail
def test_annot_to_gifti():
    assert False


def test_check_img(atlas):
    # some really basic, silly checks
    out = images.check_img(atlas['image'])
    assert out.header.get_data_dtype() == np.dtype('int32')
    assert len(out.shape) == 3

    with pytest.raises(TypeError):
        images.check_img('doesnotexist.nii.gz')

    with pytest.raises(ValueError):
        images.check_img(nib.Nifti1Image(np.zeros((5, 5, 5, 2)), np.eye(4)))


def test_check_surface(surface):
    surface = surface['image']
    # default; load images
    atlas, info = images.check_surface(surface)
    assert atlas.shape == (20484,)
    assert isinstance(info, pd.DataFrame)
    assert info.shape == (68, 3)
    assert all(info.columns == ['label', 'hemisphere', 'structure'])
    assert info.index.name == 'id'
    assert all(info['structure'] == 'cortex')

    # load pre-loaded images
    imgs = []
    for hemi in surface:
        with gzip.GzipFile(hemi) as gz:
            imgs.append(nib.GiftiImage.from_bytes(gz.read()))
    atlas2, info2 = images.check_surface(imgs)
    assert np.allclose(atlas, atlas2)
    pd.testing.assert_frame_equal(info, info2)

    # array is simply returned
    atlas3, info3 = images.check_surface(atlas)
    assert np.allclose(atlas, atlas3)
    assert info3 is None

    with pytest.raises(TypeError):
        images.check_surface(surface[0])

    with pytest.raises(TypeError):
        images.check_surface(('lh.nii.gz', 'rh.nii.gz'))


def test_check_atlas(atlas, surface):
    # check loading volumetric atlas
    tree = images.check_atlas(atlas['image'])
    assert isinstance(tree, AtlasTree)
    assert tree.atlas_info is None
    assert not tree.surface
    assert len(tree.coords) == 819621

    # check loading volumetric atlas with info
    tree = images.check_atlas(atlas['image'], atlas['info'])
    assert isinstance(tree, AtlasTree)
    assert isinstance(tree.atlas_info, pd.DataFrame)
    assert not tree.surface
    assert len(tree.coords) == 819621

    # check loading surface (info is intuited)
    tree = images.check_atlas(surface['image'])
    assert isinstance(tree, AtlasTree)
    assert isinstance(tree.atlas_info, pd.DataFrame)
    assert tree.surface
    assert len(tree.coords) == 18426


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
