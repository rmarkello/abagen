# -*- coding: utf-8 -*-
"""
Tests for abagen.images module
"""

import gzip
from pkg_resources import resource_filename

import nibabel as nib
import numpy as np
import pandas as pd
import pytest

from abagen import images
from abagen.matching import AtlasTree


@pytest.fixture(scope='module')
def annotation(tmp_path_factory):
    labels = np.tile([0, 1, 2, 3, 4], 5)
    ctab = np.asarray([
        [25, 5, 25, 0, 1639705],
        [25, 100, 40, 0, 2647065],
        [125, 100, 160, 0, 10511485],
        [100, 25, 0, 0, 6500],
        [120, 70, 50, 0, 3294840]
    ])
    names = [b'background', b'label1', b'label2', b'label3', b'label4']
    fname = tmp_path_factory.mktemp('annot') / 'test.annot'
    nib.freesurfer.write_annot(fname, labels, ctab, names, False)

    return fname


@pytest.fixture(scope='module')
def fsgeometry():
    return (
        resource_filename('abagen', 'data/fsaverage5-pial-lh.surf.gii.gz'),
        resource_filename('abagen', 'data/fsaverage5-pial-rh.surf.gii.gz'),
    )


def test_leftify_atlas(atlas):
    out = images.leftify_atlas(atlas['image'])
    assert len(np.unique(out.dataobj)) == 44
    assert np.all(np.asarray(out.dataobj)[98:] == 0)


def test_relabel_gifti(surface):
    surface = surface['image']

    # basic usage (`surface` has gap between left + right hemi for subcortex)
    lh, rh = images.relabel_gifti(surface, background=None)
    data = np.hstack((lh.agg_data(), rh.agg_data()))
    assert np.allclose(np.unique(data), np.arange(69))

    # usage with unique "background"
    lh, rh = images.relabel_gifti(surface, background=['bankssts'])
    data = np.hstack((lh.agg_data(), rh.agg_data()))
    assert np.allclose(np.unique(data), np.arange(67))

    # usage with offset
    lh, rh = images.relabel_gifti(surface, offset=100)
    assert np.allclose(np.unique(lh.agg_data())[1:], np.arange(1, 35))
    assert np.allclose(np.unique(rh.agg_data())[1:], np.arange(100, 134))


def test_annot_to_gifti(annotation):
    labels = [
        'background', 'label1', 'label2', 'label3', 'label4'
    ]
    gii = images.annot_to_gifti(annotation)
    assert np.allclose(np.unique(gii.agg_data()), np.arange(5))
    lt = gii.labeltable.get_labels_as_dict()
    assert np.allclose(list(lt.keys()), np.arange(5))
    assert np.all(list(lt.values()) == labels)


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


def test_check_atlas(atlas, surface, fsgeometry):
    # check loading volumetric atlas
    tree = images.check_atlas(atlas['image'])
    assert isinstance(tree, AtlasTree)
    assert tree.atlas_info is None
    assert tree.volumetric
    assert len(tree.coords) == 819621

    # check loading volumetric atlas with info
    tree = images.check_atlas(atlas['image'], atlas['info'])
    assert isinstance(tree, AtlasTree)
    assert isinstance(tree.atlas_info, pd.DataFrame)
    assert tree.volumetric
    assert len(tree.coords) == 819621

    # check loading surface (info is intuited)
    tree = images.check_atlas(surface['image'])
    assert isinstance(tree, AtlasTree)
    assert isinstance(tree.atlas_info, pd.DataFrame)
    assert not tree.volumetric
    assert len(tree.coords) == 18426

    tree = images.check_atlas(surface['image'], geometry=fsgeometry,
                              space='fsaverage')
    assert isinstance(tree, AtlasTree)
    assert isinstance(tree.atlas_info, pd.DataFrame)
    assert not tree.volumetric
    assert len(tree.coords) == 18426

    with pytest.raises(ValueError):
        images.check_atlas(surface['image'], geometry=fsgeometry)

    # check loading donor-specific surface file
    fp = 'data/native_dk/12876/atlas-desikankilliany-{}.label.gii.gz'
    surf = [
        resource_filename('abagen', fp.format(hemi)) for hemi in ('lh', 'rh')
    ]
    tree = images.check_atlas(surf, donor='12876')
    assert isinstance(tree, AtlasTree)
    assert isinstance(tree.atlas_info, pd.DataFrame)
    assert not tree.volumetric
    assert len(tree.coords) == 386566


def test_check_geometry(fsgeometry):
    coords, triangles = images.check_geometry(fsgeometry, 'fsaverage')
    assert len(coords) == 20484
    assert len(triangles) == 40960

    with pytest.raises(ValueError):
        images.check_geometry(fsgeometry, 'notaspace')
    with pytest.raises(ValueError):
        images.check_geometry(fsgeometry, 'fsnative', donor=None)
    with pytest.raises(TypeError):
        images.check_geometry(fsgeometry[0], 'fsaverage')


def test_check_atlas_info(atlas):
    labels = np.trim_zeros(np.unique(nib.load(atlas['image']).dataobj))

    # general usage (providing two filenames) works as expected
    out = images.check_atlas_info(atlas['info'], labels)
    assert all(out.columns == ['label', 'hemisphere', 'structure'])
    assert out.index.name == 'id'

    # can accept dataframe as input
    atlas_df = pd.read_csv(atlas['info'])
    out2 = images.check_atlas_info(atlas_df, labels)
    pd.testing.assert_frame_equal(out, out2)

    # setting ID as index of dataframe is acceptable usage
    atlas_df = atlas_df.set_index('id')
    out3 = images.check_atlas_info(atlas_df, labels)
    pd.testing.assert_frame_equal(out, out3)

    # check that coercion of different hemisphere designations works
    atlas_df.loc[atlas_df['hemisphere'] == "L", 'hemisphere'] = "lh"
    atlas_df.loc[atlas_df['hemisphere'] == "R", 'hemisphere'] = "r"
    out4 = images.check_atlas_info(atlas_df, labels)
    pd.testing.assert_frame_equal(out, out4)

    # providing labels allows for missing ids in atlas_info (i.e., does not
    # raise ValueError)
    drop_last_df = atlas_df.copy().iloc[:-1]
    out5 = images.check_atlas_info(drop_last_df, labels=range(1, 83))
    assert len(out5) == 82

    # not a filename or dataframe = failure
    with pytest.raises(TypeError):
        images.check_atlas_info([1, 2, 3], labels)

    # missing data = failure
    empty_df = pd.DataFrame(columns=['id', 'hemisphere', 'structure'])
    with pytest.raises(ValueError):
        images.check_atlas_info(empty_df, labels)

    # invalid hemisphere designations
    bad_hemi_df = atlas_df.copy()
    bad_hemi_df.loc[1, 'hemisphere'] = 'notahemisphere'
    with pytest.raises(ValueError):
        images.check_atlas_info(bad_hemi_df, labels)

    # invalid structural designation
    bad_struct_df = atlas_df.copy()
    bad_struct_df.loc[1, 'structure'] = 'notastructure'
    with pytest.raises(ValueError):
        images.check_atlas_info(bad_struct_df, labels)


def test_coerce_atlas_to_dict(testfiles, atlas):
    img, info = atlas['image'], atlas['info']
    donors = ['12876', '15496']

    # test providing single atlas file
    atl, same = images.coerce_atlas_to_dict(img, donors, info)
    assert same
    assert sorted(atl.keys()) == sorted(donors)
    imgs = list(atl.values())
    assert all(imgs[0] is a for a in imgs[1:])
    assert isinstance(imgs[0].atlas_info, pd.DataFrame)

    # test providing pre-loaded atlas file
    atl, same = images.coerce_atlas_to_dict(nib.load(img), donors, info)
    assert same
    assert sorted(atl.keys()) == sorted(donors)
    imgs = list(atl.values())
    assert all(imgs[0] is a for a in imgs[1:])
    assert isinstance(imgs[0].atlas_info, pd.DataFrame)

    # test providing dictionary for atlas
    atlas_dict = {d: img for d in donors}
    atl, same = images.coerce_atlas_to_dict(atlas_dict, donors)
    assert not same
    assert sorted(atl.keys()) == sorted(donors)
    imgs = list(atl.values())
    assert not any(imgs[0] is a for a in imgs[1:])
    assert imgs[0].atlas_info is None

    with pytest.raises(ValueError):
        images.coerce_atlas_to_dict(atlas_dict, donors + ['9861'])

    with pytest.raises(TypeError):
        images.coerce_atlas_to_dict('notanatlas', donors)
