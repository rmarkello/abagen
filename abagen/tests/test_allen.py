import pandas as pd
import pytest
from abagen import allen
from abagen.tests.utils import get_resource

ATLAS = get_resource('atlas-desikankilliany.nii.gz')
ATLAS_INFO = get_resource('atlas-desikankilliany.csv')


def test_get_expression_data(testfiles):
    out = allen.get_expression_data(testfiles, ATLAS, ATLAS_INFO, exact=False)
    assert isinstance(out, pd.DataFrame)

    with pytest.raises(KeyError):
        allen.get_expression_data({'microarray': [1, 2, 3]}, ATLAS, ATLAS_INFO)
