# -*- coding: utf-8 -*-
"""
Tests for abagen.allen module
"""

import nibabel as nib
import numpy as np
import pandas as pd
import pytest

from abagen import allen
from abagen.utils import check_img


def test_vanilla_get_expression_data(testfiles, atlas):
    out = allen.get_expression_data(atlas['image'], donors=['12876', '15496'])
    assert np.allclose(out.index, range(1, 84))
    assert out.index.name == 'label'
    assert out.columns.name == 'gene_symbol'


@pytest.mark.parametrize('opts', [
    ({'atlas_info': True}),
    ({'exact': False}),
    ({'reannotated': False}),
    ({'atlas_info': True, 'exact': False}),
])
def test_extra_get_expression_data(testfiles, atlas, opts):
    if opts.get('atlas_info'):
        opts['atlas_info'] = atlas['info']

    out = allen.get_expression_data(atlas['image'], donors=['12876', '15496'],
                                    **opts)
    assert out.index.name == 'label'
    assert out.columns.name == 'gene_symbol'


def test_individualized_get_expression_data(testfiles):
    atlas = allen.datasets.fetch_desikan_killiany(native=True)
    out = allen.get_expression_data(atlas['image'], donors=['12876', '15496'])
    assert out.index.name == 'label'
    assert out.columns.name == 'gene_symbol'


def test_get_expression_data_errors(testfiles, atlas):
    # invalid probe_selection method
    with pytest.raises(ValueError):
        allen.get_expression_data(atlas['image'], donors=['12876', '15496'],
                                  probe_selection='nonsense')

    # cannot use diff_stability with only one donor
    with pytest.raises(ValueError):
        allen.get_expression_data(atlas['image'], donors=['12876'],
                                  probe_selection='diff_stability')


def test_missing_labels(testfiles, atlas):
    # remove some labels from atlas image so numbers are non-sequential
    remove = [10, 20, 60]

    # subset atlas image
    img = check_img(atlas['image'])
    img_data = np.asarray(img.dataobj)
    for i in remove:
        img_data[img_data == i] = 0
    img = img.__class__(img_data, img.affine)

    # subset atlas info
    info = pd.read_csv(atlas['info'])
    info = info[~info['id'].isin(remove)]

    # test get expression
    out, counts = allen.get_expression_data(img, info,
                                            exact=False, return_counts=True,
                                            donors=['12876', '15496'])
    assert out.index.name == 'label'
    assert out.columns.name == 'gene_symbol'
    assert len(out) == len(info)

    assert isinstance(counts, pd.DataFrame)
    assert counts.shape == (len(info), len(testfiles))


def test_get_samples_in_mask(testfiles, atlas):
    allexp, allcoords = allen.get_samples_in_mask(donors=['12876', '15496'])
    cortexp, cortcoords = allen.get_samples_in_mask(mask=atlas['image'],
                                                    donors=['12876', '15496'])

    # exp + coords shape as expected?
    assert len(allexp) == len(allcoords) and len(cortexp) == len(cortcoords)
    assert allcoords.shape[-1] == 3 and cortcoords.shape[-1] == 3
    for df in [allexp, cortexp]:
        assert df.index.name == 'well_id'
        assert df.columns.name == 'gene_symbol'

    # providing the cortical mask (atlas) should reduce # of samples returned
    assert len(allexp) > len(cortexp)


def test_coerce_atlas_to_dict(testfiles, atlas):
    img, info = atlas['image'], atlas['info']
    donors = ['12876', '15496']

    # test providing single atlas file
    atl, inf, same = allen.coerce_atlas_to_dict(img, donors, info)
    assert same
    assert sorted(atl.keys()) == sorted(donors)
    assert isinstance(inf, pd.DataFrame)
    imgs = list(atl.values())
    assert all(imgs[0] is a for a in imgs[1:])

    # test providing pre-loaded atlas file
    atl, inf, same = allen.coerce_atlas_to_dict(nib.load(img), donors, info)
    assert same
    assert sorted(atl.keys()) == sorted(donors)
    assert isinstance(inf, pd.DataFrame)
    imgs = list(atl.values())
    assert all(imgs[0] is a for a in imgs[1:])

    # test providing dictionary for atlas
    atlas_dict = {d: img for d in donors}
    atl, inf, same = allen.coerce_atlas_to_dict(atlas_dict, donors)
    assert not same
    assert sorted(atl.keys()) == sorted(donors)
    assert inf is None
    imgs = list(atl.values())
    assert not any(imgs[0] is a for a in imgs[1:])

    with pytest.raises(ValueError):
        allen.coerce_atlas_to_dict(atlas_dict, donors + ['9861'])

    with pytest.raises(TypeError):
        allen.coerce_atlas_to_dict('notanatlas', donors)
