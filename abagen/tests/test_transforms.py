# -*- coding: utf-8 -*-
"""
Tests for abagen.transforms module
"""

import numpy as np
import pytest

from abagen import transforms


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

    assert np.all(transforms.xyz_to_ijk(xyz, affine) == ijk)
    assert np.all(transforms.ijk_to_xyz(ijk, affine) == xyz)

    with pytest.raises(ValueError):
        transforms.xyz_to_ijk([[10, 10], [20, 30]], affine)
    with pytest.raises(ValueError):
        transforms.ijk_to_xyz([[10, 10], [20, 30]], affine)
