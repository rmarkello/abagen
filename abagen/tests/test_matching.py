# -*- coding: utf-8 -*-
"""
Tests for abagen.matching module
"""

import numpy as np
import pandas as pd
import pytest

from abagen import matching


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


def test_AtlasTree(atlas, surface, testfiles):
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
