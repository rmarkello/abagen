import numpy as np
import pandas as pd
import pytest
from abagen import datasets

KEYS = [
    'microarray', 'annotation', 'pacall', 'probes', 'ontology'
]


def test_fetch_datasets(testdir):
    # check downloading for a subset of donors
    files = datasets.fetch_microarray(data_dir=str(testdir),
                                      donors=['12876'])
    assert isinstance(files, dict)
    for k in KEYS:
        assert len(files.get(k)) == 1

    # check downloading incorrect donor
    with pytest.raises(ValueError):
        datasets.fetch_microarray(donors=['notadonor'])

    files = datasets.fetch_microarray(data_dir=str(testdir),
                                      donors=None)


def test_fetch_alleninf_coords():
    coords = datasets._fetch_alleninf_coords()
    assert isinstance(coords, pd.DataFrame)
    assert coords.index.name == 'well_id'
    assert np.all(coords.columns == ['mni_x', 'mni_y', 'mni_z'])
    assert coords.shape == (3702, 3)


def test_fetch_mri():
    with pytest.raises(NotImplementedError):
        datasets.fetch_mri()
