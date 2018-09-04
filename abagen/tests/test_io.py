import os.path as op
import numpy as np
import pandas as pd
import pytest
from abagen import io


IO_FUNC = [
    dict(
        key='microarray',
        has_parq=True,
        columns=None,
        func=io.read_microarray
    ),
    dict(
        key='pacall',
        has_parq=True,
        columns=None,
        func=io.read_pacall
    ),
    dict(
        key='annotation',
        has_parq=False,
        columns=[
            'structure_id', 'slab_num', 'well_id', 'slab_type',
            'structure_acronym', 'structure_name', 'polygon_id', 'mri_voxel_x',
            'mri_voxel_y', 'mri_voxel_z', 'mni_x', 'mni_y', 'mni_z'
        ],
        func=io.read_annotation
    ),
    dict(
        key='ontology',
        has_parq=False,
        columns=[
            'id', 'acronym', 'name', 'parent_structure_id', 'hemisphere',
            'graph_order', 'structure_id_path', 'color_hex_triplet'
        ],
        func=io.read_ontology
    ),
    dict(
        key='probes',
        has_parq=False,
        columns=[
            'probe_name', 'gene_id', 'gene_symbol', 'gene_name', 'entrez_id',
            'chromosome'
        ],
        func=io.read_probes
    )
]


def test_readfiles(testfiles):
    for config in IO_FUNC:
        for fn in testfiles.get(config['key']):
            assert op.exists(fn)
            if config['has_parq'] and io.use_parq:
                assert op.exists(fn.rpartition('.csv')[0] + '.parq')

            # check loading from filepath
            data = config['func'](fn, parquet=True)
            assert isinstance(data, pd.DataFrame)

            # check loading from dataframe
            data = config['func'](data)
            assert isinstance(data, pd.DataFrame)

            if config['columns'] is not None:
                assert np.all(config['columns'] == data.columns)

        with pytest.raises(TypeError):
            config['func']([1, 2, 3])
