# -*- coding: utf-8 -*-

import pandas as pd
import pytest

from abagen.mouse import io

entries = ['id', 'acronym', 'name']


def test_fetch_allenref_genes():
    # invalid entry_type fails
    with pytest.raises(ValueError):
        io.fetch_allenref_genes(entry_type='wrong_entry_type')

    # specifying entry_type returns list
    genes = io.fetch_allenref_genes(entry_type='acronym')
    assert isinstance(genes, list)
    assert all(isinstance(f, str) for f in genes)

    # no specification returns data frame
    genes = io.fetch_allenref_genes()
    assert isinstance(genes, pd.DataFrame)
    assert genes.shape == (19991, 3)
    assert list(genes.columns) == entries


def test_fetch_allenref_structures():
    with pytest.raises(ValueError):
        io.fetch_allenref_structures(entry_type='wrong_entry_type')

    structures = io.fetch_allenref_structures(entry_type='acronym')
    assert isinstance(structures, list)
    assert all(isinstance(f, str) for f in structures)

    structures = io.fetch_allenref_structures()
    assert isinstance(structures, pd.DataFrame)
    assert structures.shape == (1327, 3)
    assert list(structures.columns) == entries


def test_fetch_rubinov2015_structures():
    with pytest.raises(ValueError):
        io.fetch_rubinov2015_structures(entry_type='wrong_entry_type')

    structures = io.fetch_rubinov2015_structures(entry_type='acronym')
    assert isinstance(structures, list)
    assert all(isinstance(f, str) for f in structures)

    structures = io.fetch_rubinov2015_structures()
    assert isinstance(structures, pd.DataFrame)
    assert structures.shape == (56, 3)
