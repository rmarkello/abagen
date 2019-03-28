# -*- coding: utf-8 -*-

import random
import pytest
from abagen.mouse import io, gene

# number of tests to make
RANDOM = random.Random(1234)
GENES = io.fetch_allenref_genes()


def test_get_gene_info():
    samples = GENES.loc[RANDOM.sample(range(len(GENES)), 10)]

    acronyms = gene.get_gene_info(id=samples['id'], attributes='acronym')
    assert sorted(acronyms['acronym']) == sorted(samples['acronym'])

    names = gene.get_gene_info(acronym=samples['acronym'], attributes='name')
    assert sorted(names['name']) == sorted(samples['name'])

    with pytest.raises(ValueError):
        gene.get_gene_info(acronym=samples['acronym'], attributes='invalid')

    # multiple attributes
    info = gene.get_gene_info(id=samples['id'], attributes=['name', 'acronym'])
    assert sorted(acronyms['acronym']) == sorted(samples['acronym'])
    assert sorted(info['name']) == sorted(samples['name'])

    # all attributes
    info = gene.get_gene_info(acronym=samples['acronym'])
    assert len(info.columns) == 15

    # exception: invalid gene identifiers
    with pytest.raises(ValueError):
        gene.get_gene_info(id=-100000)
    with pytest.raises(ValueError):
        gene.get_gene_info(acronym='notagene')
