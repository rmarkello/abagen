# -*- coding: utf-8 -*-
"""
Tests for abagen.datasets.fetchers module
"""

import pytest

from abagen.datasets import fetchers


def test_fetch_microarray():
    # test different `donors` inputs
    f1 = fetchers.fetch_microarray(donors=['12876'])
    f2 = fetchers.fetch_microarray(donors='12876')
    f3 = fetchers.fetch_microarray(donors='H0351.1009')
    f4 = fetchers.fetch_microarray(donors=None)

    # test n_proc
    fetchers.fetch_microarray(donors=['12876', '15496'], n_proc=2)

    # don't test this -- it will take a wicked long time
    # f5 = datasets.fetch_microarray(donors='all')

    assert f1 == f2 == f3 == f4
    assert len(f1) == 1  # only one donor
    for k in ['microarray', 'annotation', 'pacall', 'probes', 'ontology']:
        assert k in f1['12876']

    # check downloading incorrect donor
    with pytest.raises(ValueError):
        fetchers.fetch_microarray(donors='notadonor')

    with pytest.raises(ValueError):
        fetchers.fetch_microarray(donors=['notadonor'])


def test_fetch_raw_mri():
    # test different `donors` inputs
    f1 = fetchers.fetch_raw_mri(donors=['12876'])
    f2 = fetchers.fetch_raw_mri(donors='12876')
    f3 = fetchers.fetch_raw_mri(donors='H0351.1009')
    f4 = fetchers.fetch_raw_mri(donors=None)

    assert f1 == f2 == f3 == f4
    assert len(f1) == 1
    for k in ['t1w', 't2w']:
        assert k in f1['12876']

    # check downloading incorrect donor
    with pytest.raises(ValueError):
        fetchers.fetch_raw_mri(donors='notadonor')

    with pytest.raises(ValueError):
        fetchers.fetch_raw_mri(donors=['notadonor'])


def test_fetch_rnaseq():
    f1 = fetchers.fetch_rnaseq(donors=['9861'])
    f2 = fetchers.fetch_rnaseq(donors='9861')
    f3 = fetchers.fetch_rnaseq(donors='H0351.2001')
    f4 = fetchers.fetch_rnaseq(donors=None)

    assert f1 == f2 == f3 == f4
    assert len(f1) == 1
    for k in ['genes', 'ontology', 'counts', 'tpm', 'annotation']:
        assert k in f1['9861']

    with pytest.raises(ValueError):
        fetchers.fetch_rnaseq(donors='notadonor')

    with pytest.raises(ValueError):
        fetchers.fetch_rnaseq(donors=['9861', 'notadonor'])


@pytest.mark.parametrize('group, expected', [
    ('brain', 2413),
    ('neuron', 2530),
    ('oligodendrocyte', 1769),
    ('synaptome', 1886),
    ('layers', 46)
])
def test_get_gene_group(group, expected):
    assert len(fetchers.fetch_gene_group(group)) == expected

    with pytest.raises(ValueError):
        fetchers.fetch_gene_group('notagroup')
