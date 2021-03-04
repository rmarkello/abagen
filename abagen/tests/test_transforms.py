# -*- coding: utf-8 -*-
"""
Tests for abagen.transforms module
"""

import numpy as np
import pytest

from abagen import transforms


@pytest.mark.parametrize('xyz, fsnative, donor', [
    ([0, 0, 0], [0, 1, 1], '12876'),
    ([0, 0, 0], [1, 18, -18], '15496')
])
def test_ijk_to_fsnative(xyz, fsnative, donor):
    assert np.allclose(transforms.xyz_to_fsnative(xyz, donor), fsnative)
    assert np.allclose(transforms.fsnative_to_xyz(fsnative, donor), xyz)


# honestly no clue how to actually test this beyond just trusting the affine
# matrices provided by FreeSurfer are accurate :man_shrugging:
@pytest.mark.parametrize('mni, fsavg', [
    ([0, 0, 0], [0.0528, -1.5519, -1.2012])
])
def test_mni152_fsaverage(mni, fsavg):
    assert np.allclose(transforms.mni152_to_fsaverage(mni), fsavg)
    assert np.allclose(transforms.fsaverage_to_mni152(fsavg), mni)

    with pytest.raises(ValueError):
        transforms.mni152_to_fsaverage([[0, 0]])
    with pytest.raises(ValueError):
        transforms.fsaverage_to_mni152([[0, 0]])


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
