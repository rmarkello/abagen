# -*- coding: utf-8 -*-
"""
Tests for abagen.allen module
"""

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
    assert counts.shape == (len(info), len(testfiles['probes']))
