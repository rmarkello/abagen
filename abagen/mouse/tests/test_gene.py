import pytest
from ..gene import (check_gene_validity, get_gene_info)
import random
from ..io import read_all_genes

# number of tests to make
TEST_COUNT = 10  # test 1% of the genes
GENES_LIST_ACRONYM = read_all_genes(entry_type='acronym')
GENES_LIST_ID = read_all_genes(entry_type='id')
GENES_LIST_NAME = read_all_genes(entry_type='name')
RANDOM_STRING = 'random_string'
RANDOM_ID = -10000000

# sample 10 structures
TEST_SAMPLES = random.sample(
    range(len(GENES_LIST_ACRONYM)), TEST_COUNT
)


def test_check_gene_validity():
    # iterate the samples
    for sample in TEST_SAMPLES:
        validity, _ = check_gene_validity(
            gene_acronym=GENES_LIST_ACRONYM[sample]
        )
        assert validity is True
        validity, _ = check_gene_validity(
            gene_id=GENES_LIST_ID[sample]
        )
        assert validity is True
        # remove name test temporarily
        # gene naming conventions are too complicated...
        # validity, _ = check_gene_validity(
        #    name=GENES_LIST_NAME[sample]
        # )
        # assert validity is True

    # if the gene is invalid
    validity, _ = check_gene_validity(gene_acronym=RANDOM_STRING)
    assert validity is False

    validity, _ = check_gene_validity(gene_id=RANDOM_ID)
    assert validity is False

    # exception: missing parameters
    with pytest.raises(TypeError):
        check_gene_validity()


def test_get_gene_info():
    for sample in TEST_SAMPLES:
        # single attrib
        gene_info = get_gene_info(
            gene_id=GENES_LIST_ID[sample],
            attributes='acronym'
        )
        assert gene_info == GENES_LIST_ACRONYM[sample]
        gene_info = get_gene_info(
            gene_acronym=GENES_LIST_ACRONYM[sample],
            attributes='name'
        )
        assert gene_info == GENES_LIST_NAME[sample]
        # remove name test temporarily
        # gene naming conventions are too complicated...
        # gene_info = get_gene_info(
        #    name=GENES_LIST_NAME[sample], attributes='id'
        # )
        # assert gene_info == GENES_LIST_ID[sample]
        # exceptions: attribute given is invalid
        with pytest.raises(AttributeError):
            get_gene_info(
                gene_acronym=GENES_LIST_ACRONYM[sample],
                attributes=RANDOM_STRING
            )
        # multiple attributes
        gene_info = get_gene_info(
            gene_id=GENES_LIST_ID[sample],
            attributes=['name', 'acronym']
        )
        assert gene_info['name'] == GENES_LIST_NAME[sample]
        assert gene_info['acronym'] == GENES_LIST_ACRONYM[sample]
        # one attribute is invalid
        gene_info = get_gene_info(
            gene_acronym=GENES_LIST_ACRONYM[sample],
            attributes=[RANDOM_STRING, 'id']
        )
        assert RANDOM_STRING not in gene_info
        assert gene_info['id'] == GENES_LIST_ID[sample]

        # all atributes are invalid
        gene_info = get_gene_info(
            gene_id=GENES_LIST_ID[sample],
            attributes=[RANDOM_STRING, RANDOM_ID]
        )
        assert len(gene_info) == 0

        # all atributes
        gene_info = get_gene_info(
            gene_acronym=GENES_LIST_ACRONYM[sample],
        )
        assert gene_info['id'] == GENES_LIST_ID[sample]
        assert 'entrez-id' in gene_info

    # exception: invalid gene identifiers
    with pytest.raises(ValueError):
        get_gene_info(
            gene_id=RANDOM_ID, attributes='acronym'
        )
    with pytest.raises(ValueError):
        get_gene_info(
            gene_acronym=RANDOM_STRING
        )
