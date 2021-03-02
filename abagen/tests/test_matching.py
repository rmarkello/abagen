# -*- coding: utf-8 -*-
"""
Tests for abagen.matching module
"""

import pytest

from abagen import matching


@pytest.mark.xfail
def test_check_label():
    assert False


@pytest.mark.xfail
def test_centroids():
    assert False
    centroids = matching.get_centroids()
    assert matching.closest_centroid(centroids[0], centroids) == 0


@pytest.mark.xfail
def test_AtlasTree():
    assert False
