import pandas as pd
import pytest
from abagen import allen
from abagen.datasets import fetch_desikan_killiany

ATLAS = fetch_desikan_killiany()


def test_label_samples(testfiles):
    out = allen.label_samples(testfiles.annotation[0], ATLAS.image)
    assert isinstance(out, pd.DataFrame)
    assert out.index.name == 'sample_id'
    assert out.columns == ['label']


def test_vanilla_get_expression_data(testfiles):
    out = allen.get_expression_data(testfiles, ATLAS.image)
    assert isinstance(out, pd.DataFrame)
    assert out.index.name == 'label'
    assert out.columns.name == 'gene_symbol'

    with pytest.raises(KeyError):
        allen.get_expression_data({'microarray': [1, 2, 3]}, ATLAS.image)


def test_extra_get_expression_data(testfiles):
    for opts in [{'atlas_info': ATLAS.info},
                 {'exact': False},
                 {'atlas_info': ATLAS.info, 'exact': False}]:
        out = allen.get_expression_data(testfiles, ATLAS.image, **opts)
        assert isinstance(out, pd.DataFrame)
        assert out.index.name == 'label'
        assert out.columns.name == 'gene_symbol'
