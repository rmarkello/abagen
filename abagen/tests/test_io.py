# -*- coding: utf-8 -*-
"""
Tests for abagen.io module
"""

import os.path as op

import numpy as np
import pandas as pd
import pytest

from abagen import io


annot_cols = [
    'structure_id', 'slab_num', 'well_id', 'slab_type', 'structure_acronym',
    'structure_name', 'polygon_id', 'mri_voxel_x', 'mri_voxel_y',
    'mri_voxel_z', 'mni_x', 'mni_y', 'mni_z'
]

ont_cols = [
    'id', 'acronym', 'name', 'parent_structure_id', 'hemisphere',
    'graph_order', 'structure_id_path', 'color_hex_triplet'
]

probe_cols = [
    'probe_name', 'gene_id', 'gene_symbol', 'gene_name', 'entrez_id',
    'chromosome'
]


@pytest.mark.parametrize('key, has_parq, columns', [
    ('microarray', True, None),
    ('pacall', True, None),
    ('annotation', False, annot_cols),
    ('ontology', False, ont_cols),
    ('probes', False, probe_cols),
])
def test_readfiles(testfiles, key, has_parq, columns):
    for fn in testfiles.get(key):
        func = getattr(io, 'read_{}'.format(key))

        # check file (CSV + parquet) exist
        assert op.exists(fn)
        if has_parq and io.use_parq:
            assert op.exists(fn.rpartition('.csv')[0] + '.parq')

        # check loading from filepath
        data = func(fn, parquet=True) if has_parq else func(fn)
        assert isinstance(data, pd.DataFrame)

        # check loading from dataframe (should return same object)
        data2 = func(data)
        assert id(data) == id(data2)

        # check that copy parameter works as expected
        data3 = func(data, copy=True)
        assert isinstance(data3, pd.DataFrame)
        assert id(data) != id(data3)

        # confirm columns are as expected
        if columns is not None:
            assert np.all(columns == data.columns)

        # confirm errors
        with pytest.raises(TypeError):
            func(1)

        with pytest.raises(TypeError):
            func([1, 2, 3])

        with pytest.raises(FileNotFoundError):
            func('notafile')
