import os.path as op
import numpy as np
import pandas as pd
import pytest
from abagen import datasets, io

KEYS = [
    'microarray', 'annotation', 'pacall', 'probes', 'ontology'
]


def test_fetch_datasets(tmpdir):
    # try downloading donor w/smallest data
    files = datasets.fetch_microarray(data_dir=str(tmpdir), donors=['12876'])
    assert isinstance(files, dict)
    for k in KEYS:
        assert hasattr(files, k)
        assert len(files.get(k)) == 1
    # check if parquet was created
    if io.use_parq:
        assert op.exists(files.microarray[0].rpartition('.csv')[0] + '.parq')
    with pytest.raises(ValueError):
        datasets.fetch_microarray(donors=['notadonor'])


def test_fetch_alleninf_coords():
    coords = datasets._fetch_alleninf_coords()
    assert isinstance(coords, pd.DataFrame)
    assert coords.index.name == 'well_id'
    assert np.all(coords.columns == ['mni_x', 'mni_y', 'mni_z'])
    assert coords.shape == (3702, 3)


def test_fetch_mri():
    with pytest.raises(NotImplementedError):
        datasets.fetch_mri()
