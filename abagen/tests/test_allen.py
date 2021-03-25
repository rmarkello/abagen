# -*- coding: utf-8 -*-
"""
Tests for abagen.allen module
"""

import nibabel as nib
import numpy as np
import pandas as pd
import pytest

from abagen import allen
from abagen.images import check_img, check_surface


def test_vanilla_get_expression_data(testfiles, atlas):
    out = allen.get_expression_data(atlas['image'], donors=testfiles.keys())
    assert np.allclose(out.index, range(1, 84))
    assert out.index.name == 'label'
    assert out.columns.name == 'gene_symbol'


def test_vanilla_surface_expression_data(testfiles, surface):
    out = allen.get_expression_data(surface['image'], donors=testfiles.keys())
    assert out.index.name == 'label'
    assert out.columns.name == 'gene_symbol'


@pytest.mark.parametrize('opts', [
    ({'atlas_info': True}),
    ({'missing': 'centroids'}),
    ({'atlas_info': True, 'missing': 'centroids'}),
])
def test_extra_get_expression_data(testfiles, atlas, opts):
    if opts.get('atlas_info', False):
        opts['atlas_info'] = atlas['info']

    out = allen.get_expression_data(atlas['image'], donors=testfiles.keys(),
                                    **opts)
    assert out.index.name == 'label'
    assert out.columns.name == 'gene_symbol'
    if 'missing' in opts and 'atlas_info' not in opts:
        assert not np.any(out.isna())
    elif 'missing' in opts and 'atlas_info' in opts:
        assert np.any(out.isna())


def test_individualized_get_expression_data(testfiles):
    atlas = allen.datasets.fetch_desikan_killiany(native=True)
    out = allen.get_expression_data(atlas['image'], donors=testfiles.keys())
    assert out.index.name == 'label'
    assert out.columns.name == 'gene_symbol'


def test_individualized_surface_expression_data(testfiles):
    atlas = allen.datasets.fetch_desikan_killiany(native=True, surface=True)
    out = allen.get_expression_data(atlas['image'], donors=testfiles.keys())
    assert out.index.name == 'label'
    assert out.columns.name == 'gene_symbol'


def test_get_expression_data_errors(testfiles, atlas):
    # invalid probe_selection method
    with pytest.raises(ValueError):
        allen.get_expression_data(atlas['image'], donors=testfiles.keys(),
                                  probe_selection='nonsense')

    # cannot use diff_stability with only one donor
    with pytest.raises(ValueError):
        donor = list(testfiles.keys())[0]
        allen.get_expression_data(atlas['image'], donors=[donor],
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
                                            missing='interpolate',
                                            return_counts=True,
                                            donors=testfiles.keys())
    assert out.index.name == 'label'
    assert out.columns.name == 'gene_symbol'
    assert len(out) == len(info)

    assert isinstance(counts, pd.DataFrame)
    assert counts.shape == (len(info), len(testfiles))


def test_get_samples_in_mask(testfiles, atlas):
    allexp, allcoords = allen.get_samples_in_mask(donors=testfiles.keys())
    cortexp, cortcoords = allen.get_samples_in_mask(mask=atlas['image'],
                                                    donors=testfiles.keys())

    # exp + coords shape as expected?
    assert len(allexp) == len(allcoords) and len(cortexp) == len(cortcoords)
    assert allcoords.shape[-1] == 3 and cortcoords.shape[-1] == 3
    for df in [allexp, cortexp]:
        assert df.index.name == 'well_id'
        assert df.columns.name == 'gene_symbol'

    # providing the cortical mask (atlas) should reduce # of samples returned
    assert len(allexp) > len(cortexp)


def test_get_interpolated_map(testfiles, atlas, surface):
    img = np.asarray(nib.load(atlas['image']).dataobj)
    out = allen.get_interpolated_map('A1BG', mask=atlas['image'],
                                     donors=testfiles.keys())
    assert isinstance(out, dict)
    assert len(out) == 1
    assert 'A1BG' in out
    assert out['A1BG'].shape == img.shape
    assert np.allclose(out['A1BG'].nonzero(), img.nonzero())

    img = check_surface(surface['image'])[0]
    out = allen.get_interpolated_map(('A1BG', 'ZZZ3'), mask=surface['image'],
                                     donors=testfiles.keys())
    assert isinstance(out, dict)
    assert len(out) == 2
    for f in ('A1BG', 'ZZZ3'):
        assert f in out
        assert out[f].shape == img.shape
        assert np.allclose(out[f].nonzero(), img.nonzero())
