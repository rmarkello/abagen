# -*- coding: utf-8 -*-
"""
Tests for abagen.surfaces module
"""

import numpy as np
import pytest
from scipy import sparse

from abagen import datasets, surfaces


@pytest.fixture(scope='module')
def surf():
    data = datasets.fetch_fsaverage5()
    coords = np.row_stack([hemi.vertices for hemi in data])
    triangles, offset = [], 0
    for hemi in data:
        triangles.append(hemi.faces + offset)
        offset += hemi.vertices.shape[0]
    triangles = np.row_stack(triangles)

    return coords, triangles


def test_make_surf_graph(surf):
    coords = np.array([[0, 0, 0], [0, 0, 1], [0, 1, 1], [1, 1, 0], [1, 1, 1]])
    triangles = np.array([[0, 1, 2], [0, 1, 3], [2, 3, 4]])
    graph = surfaces.make_surf_graph(coords, triangles)
    assert isinstance(graph, sparse.csr_matrix)
    assert graph.shape[0] == coords.shape[0]
    assert graph.nnz == 8
    assert np.allclose(graph[[0, 1], [-1, -1]], 0)

    # real graph
    coords, triangles = surf
    graph = surfaces.make_surf_graph(coords, triangles)
    assert isinstance(graph, sparse.csr_matrix)
    assert graph.shape[0] == coords.shape[0]


def test_get_graph_distance(surf):
    coords = np.array([[0, 0, 0], [0, 0, 1], [0, 1, 1], [1, 1, 0], [1, 1, 1]])
    triangles = np.array([[0, 1, 2], [0, 1, 3], [2, 3, 4]])
    graph = surfaces.make_surf_graph(coords, triangles)
    dist = surfaces.get_graph_distance(graph)
    assert dist.shape == (5, 5)
    assert dist[0, 4] == dist[0, 2] + dist[2, 4]

    dist = surfaces.get_graph_distance(graph, nodes=[0, 1])
    assert dist.shape == (2, 5)

    dist = surfaces.get_graph_distance(graph, nodes=[0, 1],
                                       labels=[1, 1, 1, 2, 2])
    assert dist.shape == (2, 2)

    # real graph
    graph = surfaces.make_surf_graph(*surf)
    dist = surfaces.get_graph_distance(graph, nodes=[0, 1, 2, 3])
    assert dist.shape == (4, graph.shape[1])
