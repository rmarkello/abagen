import pandas as pd
import pytest
from abagen import allen
from abagen.datasets import fetch_desikan_killiany

ATLAS = fetch_desikan_killiany()


def test_get_expression_data(testfiles):
    kwargs = {}
    for opts in [{}, {'atlas_info': ATLAS.info}, {'exact': False}]:
        kwargs.update(opts)
        out = allen.get_expression_data(testfiles, ATLAS.image, **kwargs)
        assert isinstance(out, pd.DataFrame)
        assert out.index.name == 'label'
        assert out.columns.name == 'gene_symbol'

    with pytest.raises(KeyError):
        allen.get_expression_data({'microarray': [1, 2, 3]}, ATLAS.image)
