import numpy as np
import pandas as pd
import pytest
from abagen import utils
from abagen.datasets import fetch_desikan_killiany

ATLAS = fetch_desikan_killiany()
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


def test_check_atlas_info(atlas):
    # general usage (providing two filenames) works as expected
    out = utils.check_atlas_info(atlas['image'], atlas['info'])
    assert all(out.columns == ['label', 'hemisphere', 'structure'])
    assert out.index.name == 'id'

    # can accept dataframe as input
    atlas_df = pd.read_csv(atlas['info'])
    out2 = utils.check_atlas_info(atlas['image'], atlas_df)
    pd.testing.assert_frame_equal(out, out2)

    # setting ID as index of dataframe is acceptable usage
    atlas_df = atlas_df.set_index('id')
    out3 = utils.check_atlas_info(atlas['image'], atlas_df)
    pd.testing.assert_frame_equal(out, out3)

    # check that coercion of different hemisphere designations works
    atlas_df.loc[atlas_df['hemisphere'] == "L", 'hemisphere'] = "lh"
    atlas_df.loc[atlas_df['hemisphere'] == "R", 'hemisphere'] = "r"
    out4 = utils.check_atlas_info(atlas['image'], atlas_df)
    pd.testing.assert_frame_equal(out, out4)

    # validate = True returns None
    none = utils.check_atlas_info(atlas['image'], atlas['info'], validate=True)
    assert none is None

    # providing labels allows for missing ids in atlas_info (i.e., does not
    # raise ValueError)
    drop_last_df = atlas_df.copy().iloc[:-1]
    out5 = utils.check_atlas_info(atlas['image'], drop_last_df,
                                  labels=range(1, 83))
    assert len(out5) == 82

    # not a filename or dataframe = failure
    with pytest.raises(TypeError):
        utils.check_atlas_info(atlas['image'], [1, 2, 3])

    # missing data = failure
    empty_df = pd.DataFrame(columns=['id', 'hemisphere', 'structure'])
    with pytest.raises(ValueError):
        utils.check_atlas_info(atlas['image'], empty_df)

    # invalid hemisphere designations
    bad_hemi_df = atlas_df.copy()
    bad_hemi_df.loc[1, 'hemisphere'] = 'notahemisphere'
    with pytest.raises(ValueError):
        utils.check_atlas_info(atlas['image'], bad_hemi_df)

    # invalid structural designation
    bad_struct_df = atlas_df.copy()
    bad_struct_df.loc[1, 'structure'] = 'notastructure'
    with pytest.raises(ValueError):
        utils.check_atlas_info(atlas['image'], bad_struct_df)


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
