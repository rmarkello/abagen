from ..io import read_all_genes, read_all_structures
import numpy as np
import pytest
import random


def test_read_all_genes():
    with pytest.raises(KeyError):
        read_all_genes(entry_type='wrong_entry_type')

    all_genes = read_all_genes(entry_type='id')
    assert isinstance(all_genes, np.ndarray)

    all_genes = read_all_genes(
        entry_type=['id', 'acronym', 'name']
    )
    assert all_genes.shape == (19991, 3)


def test_real_all_structures():
    with pytest.raises(KeyError):
        read_all_structures(entry_type='wrong_entry_type')
    default_structures = \
        read_all_structures(entry_type='acronym')
    assert isinstance(default_structures, np.ndarray)
    assert isinstance(
        random.choice(default_structures), str
    )
    default_structures = \
        read_all_structures(
            entry_type=['id', 'acronym', 'name']
        )
    assert default_structures.shape == (56, 3)


