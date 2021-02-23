# -*- coding: utf-8 -*-
"""
Tests for abagen.images module
"""

import numpy as np
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


@pytest.mark.xfail
def test_get_unique_labels():
    assert False


@pytest.mark.xfail
def test_get_centroids():
    assert False


@pytest.mark.xfail
def test_closest_centroid():
    assert False


@pytest.mark.xfail
def test_expand_roi():
    assert False
