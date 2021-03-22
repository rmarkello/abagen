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
    coords, triangles = surf
    graph = surfaces.make_surf_graph(coords, triangles)
    assert isinstance(graph, sparse.csr_matrix)
    assert graph.shape[0] == coords.shape[0]


def test_get_graph_distance(surf):
    graph = surfaces.make_surf_graph(*surf)
    dist = surfaces.get_graph_distance(graph, nodes=[0, 1, 2, 3])
    assert dist.shape == (4, graph.shape[1])
