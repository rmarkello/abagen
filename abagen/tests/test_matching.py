# -*- coding: utf-8 -*-
"""
Tests for abagen.matching module
"""

import nibabel as nib
import numpy as np
import pandas as pd
import pytest

from abagen import images, matching, transforms


def test_check_label():
    sample_info = pd.DataFrame(dict(
        hemisphere=['R', 'R', 'L', 'L'],
        structure=['cortex', 'cortex', 'cortex', 'subcortex/brainstem']
    ), index=pd.Series([1, 2, 3, 4], name='label'))
    atlas_info = pd.DataFrame(dict(
        hemisphere=['R', 'L', 'L', 'B'],
        structure=['cortex', 'cortex', 'cortex', 'subcortex/brainstem']
    ), index=pd.Series([1, 2, 3, 4], name='label'))

    labels, expected = [1, 2, 3, 4], [1, 0, 3, 4]
    out = matching._check_label(labels, sample_info, atlas_info)
    assert np.allclose(out, expected)

    col = ['hemisphere']
    out = matching._check_label(labels, sample_info, atlas_info[col])
    assert np.allclose(out, labels)
    out = matching._check_label(labels, sample_info[col], atlas_info)
    assert np.allclose(out, labels)


def test_get_centroids():
    # basic test data
    data = np.array([1, 1, 1, 2, 2, 2])
    coords = np.array([
        [0, 0, 0], [1, 1, 1], [2, 2, 2],
        [3, 3, 3], [4, 4, 4], [5, 5, 5]
    ])
    expected = {
        1: [1, 1, 1], 2: [4, 4, 4]
    }

    # check output is as expected
    centroids = matching.get_centroids(data, coords)
    for k, v in centroids.items():
        assert k in expected and np.allclose(expected[k], v)

    # providing labels only returns centroids with corresponding label
    centroids = matching.get_centroids(data, coords, labels=[1])
    assert len(centroids) == 1 and np.allclose(centroids[1], expected[1])


def test_closest_centroids():
    centroids = np.array([[1, 1, 1], [4, 4, 4]])
    out, dist = matching.closest_centroid(centroids[0], centroids, True)
    assert len(out) == len(dist) == 1
    assert np.allclose(out, 0) and np.allclose(dist, 0)

    out = matching.closest_centroid([0, 0, 0], centroids)
    assert np.allclose(out, 0)


def test_AtlasTree(atlas, surface):
    # basic test data
    data = np.array([1, 1, 1, 2, 2, 2])
    coords = np.array([
        [0, 0, 0], [1, 1, 1], [2, 2, 2],
        [3, 3, 3], [4, 4, 4], [5, 5, 5]
    ])

    atlas_info = pd.DataFrame(dict(
        hemisphere=['R', 'L'],
        structure=['cortex', 'cortex']
    ), index=pd.Series([1, 2], name='id'))

    # check basic properties of tree
    tree = matching.AtlasTree(data, coords)
    assert str(tree) == 'AtlasTree[n_rois=2, n_vertex=6]'
    assert not tree.volumetric
    assert tree.atlas_info is None
    assert np.all(tree.atlas == data)
    assert np.all(tree.coords == coords)
    assert np.all(tree.labels == [1, 2])
    assert len(tree.centroids) == 2 and list(tree.centroids.keys()) == [1, 2]
    tree.atlas_info = atlas_info
    pd.testing.assert_frame_equal(tree.atlas_info, atlas_info)

    # check sample matching
    labels = tree.label_samples([[0, 0, 0], [3, 3, 3]])
    assert np.all(labels['label'] == [1, 2])

    # check coordinate assignment AND outlier removal in same go
    tree.coords = coords[::-1]
    labels = tree.label_samples(np.row_stack((coords, [1000, 1000, 1000])))
    assert np.all(labels['label'] == [2, 2, 2, 1, 1, 1, 0])

    with pytest.raises(ValueError):
        tree.coords = coords[:-1]

    # check centroid matching
    lab, dist = tree.match_closest_centroids([[-1, -1, -1]], return_dist=True)
    assert np.all(lab == 2)
    assert np.allclose(dist, np.sqrt(3) * 2)
    centinfo = pd.DataFrame(dict(mni_x=[-1], mni_y=[-1], mni_z=[-1],
                                 hemisphere='L', structure='cortex'))
    lab, dist = tree.match_closest_centroids([[-1, -1, -1]], return_dist=True)
    assert np.all(lab == 2)
    assert np.allclose(dist, np.sqrt(3) * 2)
    centinfo['structure'] = 'cerebellum'
    lab, dist = tree.match_closest_centroids(centinfo, return_dist=True)
    assert np.all(lab == -1)
    assert np.all(np.isinf(dist))

    # check label filling
    lab1 = tree.fill_label([[0, 0, 0], [3, 3, 3]], label=2)
    lab2, dist = tree.fill_label([[0, 0, 0], [3, 3, 3]], label=2,
                                 return_dist=True)
    assert np.allclose(lab1, lab2)
    assert np.allclose(lab2, [1, 0, 0])
    assert np.allclose(dist, [np.sqrt(3), np.sqrt(3), 0])

    # check providing niimg-like atlas
    tree = matching.AtlasTree(nib.load(atlas['image']))
    assert str(tree) == 'AtlasTree[n_rois=83, n_voxel=819621]'
    assert tree.volumetric
    assert tree.atlas_info is None
    assert np.all(tree.labels == np.arange(1, 84))
    assert len(tree.centroids) == 83
    tree.atlas_info = atlas['info']
    assert isinstance(tree.atlas_info, pd.DataFrame)
    labels = tree.label_samples([tree.centroids[1], tree.centroids[2]])
    assert np.all(labels['label'] == [1, 2])

    # coordinates supplied with volumetric image
    with pytest.warns(UserWarning):
        matching.AtlasTree(nib.load(atlas['image']), coords=coords)

    # check surface AtlasTree
    tree = images.check_atlas(surface['image'])
    assert str(tree) == 'AtlasTree[n_rois=68, n_vertex=18426]'
    assert not tree.volumetric
    assert tree.triangles is not None
    assert tree.graph is not None
    assert isinstance(tree.atlas_info, pd.DataFrame)
    assert len(tree.centroids) == 68
    labels = tree.label_samples([tree.centroids[1], tree.centroids[2]])
    assert np.all(labels['label'] == [1, 2])

    with pytest.raises(ValueError):
        tree.triangles = [[512423, 512312, 4213215]]

    # check negative surface tolerance
    labels = tree.label_samples([-72, -25, -13], tolerance=-4)
    assert np.all(labels['label'] == 14)
    labels = tree.label_samples([-72, -25, -13], tolerance=-3)
    assert np.all(labels['label'] == 0)

    # no coordinates
    with pytest.raises(ValueError):
        matching.AtlasTree(np.random.choice(10, size=(100,)))

    # different length coordinates
    with pytest.raises(ValueError):
        matching.AtlasTree(np.random.choice(10, size=(100,)),
                           coords=np.random.rand(99, 3))


def test_nonint_voxels(atlas):
    coord = [[-56.8, -50.6, 8.8]]
    affine = np.zeros((4, 4))
    affine[:-1, :-1] = np.eye(3) * 1.5
    affine[:, -1] = nib.load(atlas['image']).affine[:, -1]
    dataobj = np.zeros((100, 100, 100))
    i, j, k = transforms.xyz_to_ijk(coord, affine).squeeze()
    dataobj[i, j, k] = 1
    newatl = nib.Nifti1Image(dataobj, affine)

    tree = matching.AtlasTree(newatl, group_atlas=True)
    assert tree.volumetric
    assert tree.label_samples(coord, tolerance=0).loc[0, 'label'] == 1
