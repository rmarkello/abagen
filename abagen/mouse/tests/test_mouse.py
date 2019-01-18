from ..mouse import (
    get_experiment_id_from_gene,
    _get_single_unionization_attribute,
    get_unionization_from_experiment,
    get_unionization_from_gene
)

import pytest
import numpy as np
import requests
from xml.etree import ElementTree as ET

RANDOM_ID = -1
RANDOM_STRING = 'random_string'
TEST_GENE_ACRONYM = {
    'Snca': {
        'sagittal': [988, 990, 79904550],
        'coronal': [986, 989, 79908848]
    },  # mutiple experiment IDs found
    'Gba': {
        'sagittal': [1612], 'coronal': [1611]
    },  # one experiment ID found
    'Elf4': {
        'sagittal': [73834415, 77464840], 'coronal': []
    }  # no experiment ID found (while the gene is valid)
}
TEST_GENE_ID = {
    84193: {
        'sagittal': [70238925, 71147924, 71213117],
        'coronal': [74047443]
    },  # one experiment ID found
    18608: {
        'sagittal': [69289721], 'coronal': []
    }  # no experiment ID found (while the gene is valid)
}
# test single experiment ID
TEST_STRUCTURE = [182305713, 182305709,182305705, 15564]
TEST_ATTRIBUTE = ['expression-energy', 'expression-density', 'sum-pixels']
TEST_EXPERIMENT = {
    986: {
        'expression-energy': np.array([7.73432, 7.28206, 3.82741, 2.52219]),
        'expression-density':np.array([0.0603072, 0.0553335, 0.0298628, 0.0188711]),
        'sum-pixels':np.array([419628.0, 2238000.0, 3629050.0, 1534150000.0])
    },
    69782969: {
        'expression-energy': np.array([1.00879, 1.40198, 2.34988, 2.729]),
        'expression-density':np.array([0.00650952, 0.00896709, 0.015409, 0.0178506]),
        'sum-pixels':np.array([195863.0, 1051230.0, 1440430.0, 998731000.0])
    }
}


def test_get_experiment_id_from_gene():
    # exception: gene is invalid
    with pytest.raises(ValueError):
        get_experiment_id_from_gene(
            gene_acronym=RANDOM_STRING,
            slicing_direction='coronal'
        )
    with pytest.raises(ValueError):
        get_experiment_id_from_gene(
            gene_id=RANDOM_ID
        )
    # test acronym entry type
    for gene in TEST_GENE_ACRONYM:
        # test sagittal
        experiment_id_list = get_experiment_id_from_gene(
            gene_acronym=gene
        )
        assert len(experiment_id_list) == len(
            TEST_GENE_ACRONYM[gene]['sagittal']
        )
        for experiment_id in experiment_id_list:
            assert experiment_id in TEST_GENE_ACRONYM[gene]['sagittal']
        # exception: slicing direction is invalid
        with pytest.raises(ValueError):
            get_experiment_id_from_gene(
                gene_acronym=gene,
                slicing_direction=RANDOM_STRING
            )

        # test coronal
        experiment_id_list = get_experiment_id_from_gene(
            gene_acronym=gene,
            slicing_direction='coronal'
        )
        assert len(experiment_id_list) == len(
            TEST_GENE_ACRONYM[gene]['coronal']
        )
        for experiment_id in experiment_id_list:
            assert experiment_id in TEST_GENE_ACRONYM[gene]['coronal']

    # test id entry type
    for gene in TEST_GENE_ID:
        # test sagittal
        experiment_id_list = get_experiment_id_from_gene(
            gene_id=gene, slicing_direction='sagittal'
        )
        assert len(experiment_id_list) == len(
            TEST_GENE_ID[gene]['sagittal']
        )
        for experiment_id in experiment_id_list:
            assert experiment_id in TEST_GENE_ID[gene]['sagittal']
        # test coronal
        experiment_id_list = get_experiment_id_from_gene(
            gene_id=gene,
            slicing_direction='coronal'
        )
        assert len(experiment_id_list) == len(
            TEST_GENE_ID[gene]['coronal']
        )

        # exception: slicing direction is invalid
        with pytest.raises(ValueError):
            get_experiment_id_from_gene(
                gene_id=gene,
                slicing_direction=RANDOM_STRING
            )


def test_get_single_unionization_attribute():
    for experiment_id in TEST_EXPERIMENT:
        url = "http://api.brain-map.org/api/v2/data/SectionDataSet/query.xml?id={}" \
            "&include=structure_unionizes%28structure%29"\
            .format(experiment_id)
        r = requests.get(url)
        root = ET.fromstring(r.content)
        with pytest.raises(AttributeError):
            _get_single_unionization_attribute(
                root=root,
                attr=RANDOM_STRING,
                structure_list=TEST_STRUCTURE
            )
        vals = _get_single_unionization_attribute(
            root=root, attr='expression-energy',
            structure_list=TEST_STRUCTURE
        )
        assert np.allclose(vals, TEST_EXPERIMENT[experiment_id]['expression-energy'])
        vals = _get_single_unionization_attribute(
            root=root, attr='expression-density',
            structure_list=TEST_STRUCTURE
        )
        assert np.allclose(vals, TEST_EXPERIMENT[experiment_id]['expression-density'])
        vals = _get_single_unionization_attribute(
            root=root, attr='sum-pixels',
            structure_list=TEST_STRUCTURE
        )
        assert np.allclose(vals, TEST_EXPERIMENT[experiment_id]['sum-pixels'])


def test_get_unionization_from_experiment():
    with pytest.raises(ValueError):
        get_unionization_from_experiment(
            RANDOM_ID, structure_list=TEST_STRUCTURE
        )
    for experiment_id in TEST_EXPERIMENT:
        unionization = get_unionization_from_experiment(
            experiment_id=experiment_id,
            structure_list=TEST_STRUCTURE, attributes='expression-energy'
        )
        assert np.allclose(
            unionization,
            TEST_EXPERIMENT[experiment_id]['expression-energy']
        )
        unionization = get_unionization_from_experiment(
            experiment_id=experiment_id,
            structure_list=TEST_STRUCTURE
        )
        assert np.allclose(
            unionization['expression-density'],
            TEST_EXPERIMENT[experiment_id]['expression-density']
        )
        assert np.allclose(
            unionization['sum-pixels'],
            TEST_EXPERIMENT[experiment_id]['sum-pixels']
        )


def test_get_unionization_from_gene():
    pass
