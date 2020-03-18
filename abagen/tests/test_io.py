# -*- coding: utf-8 -*-
"""
Tests for abagen.io module
"""

import os.path as op

import numpy as np
import pandas as pd
import pytest

from abagen import io
from abagen.utils import flatten_dict


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
gene_cols = [
    'gene_id', 'entrez_id', 'chromosome', 'strand', 'number_of_transcripts',
    'median_transcriptome_length', 'median_genome_length',
    'median_number_of_exons', 'median_gene_start', 'median_gene_end'
]
rna_annot_cols = [
    'RNAseq_sample_name', 'replicate_sample', 'sample_name', 'well_id',
    'microarray_run_id', 'ontology_color', 'main_structure',
    'sub_structure', 'structure_id', 'structure_acronym', 'hemisphere',
    'brain', 'million_clusters', 'clip_percentage', 'RIN_RNA_quality',
    'rnaseq_run_id', 'A.Pct', 'C.Pct', 'G.Pct', 'T.Pct', 'N.Pct'
]


@pytest.mark.parametrize('key, has_parq, columns', [
    ('microarray', True, None),
    ('pacall', True, None),
    ('annotation', False, annot_cols),
    ('ontology', False, ont_cols),
    ('probes', False, probe_cols),
])
def test_readfiles(testfiles, key, has_parq, columns):
    for d, fn in flatten_dict(testfiles, key).items():
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


@pytest.mark.parametrize('key, columns', [
    ('genes', gene_cols),
    ('ontology', ont_cols),
    ('counts', None),
    ('tpm', None),
    ('annotation', rna_annot_cols)
])
def test_readrnaseq(rnafiles, key, columns):
    for d, fn in flatten_dict(rnafiles, key).items():
        func = getattr(io, 'read_{}'.format(key))

        # check file exists
        assert op.exists(fn)

        # check loading from filepath
        data = func(fn)
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
