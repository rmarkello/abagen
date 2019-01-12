from ..mouse import get_experiment_id_from_gene
import pytest
import numpy as np

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
TEST_ATTRIBUTE = ['expression-energy', 'expression-density', 'sum-pixels']
TEST_EXPERIMENT = {
    69782969: {'expression-energy': np.array([]), },

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


def test_get_unionization_from_experiment():
    pass


def test_get_unionization_from_gene():
    pass
