import numpy as np
import pytest
from abagen import utils

EXAMPLE_COORDS = dict(
    affine=np.array([[1, 0, 0, -90],
                     [0, 1, 0, -150],
                     [0, 0, 1, -80],
                     [0, 0, 0, 1]]),
    coords=[
        dict(
            ijk=[0, 0, 0],
            xyz=np.array([[-90, -150, -80]])
        ),
        dict(
            ijk=[[10, 10, 10], [100, 50, 100]],
            xyz=np.array([[-80, -140, -70], [10, -100, 20]])
        ),
        dict(
            ijk=[[54, 32, 20], [82, 205, 38], [32, 51, 82]],
            xyz=np.array([[-36, -118, -60], [-8, 55, -42], [-58, -99, 2]])
        )
    ]
)


def test_coords_transform():
    aff = EXAMPLE_COORDS['affine']
    for coords in EXAMPLE_COORDS['coords']:
        ijk, xyz = coords['ijk'], coords['xyz']
        assert np.all(utils.xyz_to_ijk(xyz, aff) == ijk)
        assert np.all(utils.ijk_to_xyz(ijk, aff) == xyz)
    with pytest.raises(ValueError):
        utils.xyz_to_ijk([[10, 10], [20, 30]], aff)
    with pytest.raises(ValueError):
        utils.ijk_to_xyz([[10, 10], [20, 30]], aff)
