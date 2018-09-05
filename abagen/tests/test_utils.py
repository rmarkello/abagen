import numpy as np
import pandas as pd
import pytest
from abagen import utils
from abagen.tests.utils import get_resource

ATLAS = get_resource('atlas-desikankilliany.nii.gz')
ATLAS_INFO = get_resource('atlas-desikankilliany.csv')
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
EXAMPLE_METRICS = [
    dict(
        metric='mean',
        check=np.array([[1, 2], [3, 4], [5, 6]]),
        confirm=np.array([1.5, 3.5, 5.5]),
        kwargs=dict(axis=1)
    ),
    dict(
        metric='median',
        check=np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]),
        confirm=np.array([2, 5, 8]),
        kwargs=dict(axis=1)
    )
]


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


def test_check_atlas_info():
    # check appropriate usage
    out = utils.check_atlas_info(ATLAS, ATLAS_INFO)
    assert isinstance(out, pd.DataFrame)
    atlas_df = pd.read_csv(ATLAS_INFO)
    out = utils.check_atlas_info(ATLAS, atlas_df)
    atlas_df = atlas_df.set_index('id')
    out = utils.check_atlas_info(ATLAS, atlas_df)

    # check bad usage
    with pytest.raises(ValueError):
        utils.check_atlas_info(ATLAS, [1, 2, 3])

    bad_df = pd.DataFrame(columns=['id', 'hemisphere', 'structure'])
    with pytest.raises(ValueError):
        utils.check_atlas_info(ATLAS, bad_df)


def test_check_metric():
    for config in EXAMPLE_METRICS:
        metric = utils.check_metric(config['metric'])
        check, kwargs = config.get('check'), config.get('kwargs', {})
        assert np.allclose(metric(check, **kwargs), config['confirm'])

    # not an appropriate str
    with pytest.raises(ValueError):
        utils.check_metric('blargh')

    # doesn't accept axis kwarg
    with pytest.raises(TypeError):
        utils.check_metric(lambda x: np.mean(x))
