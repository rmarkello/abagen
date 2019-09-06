import os
import shutil

import pytest

from abagen import datasets


def test_get_dataset_dir(testdir):
    os.environ.pop('ABAGEN_DATA', None)

    # check that data dir defaults to $HOME/abagen-data assuming no env var
    expected_base = os.path.expanduser('~/abagen-data')
    data_dir = datasets._get_dataset_dir('test', verbose=0)
    assert data_dir == os.path.join(expected_base, 'test')
    assert os.path.isdir(data_dir) and os.path.exists(data_dir)
    shutil.rmtree(data_dir)

    # if env var is set, we should default to that
    expected_base = os.path.join(testdir, 'test-abagen-data')
    os.environ['ABAGEN_DATA'] = expected_base
    data_dir = datasets._get_dataset_dir('test', verbose=0)
    assert data_dir == os.path.join(expected_base, 'test')
    assert os.path.isdir(data_dir) and os.path.exists(data_dir)
    shutil.rmtree(data_dir)

    # test explicitly setting the data_dir
    expected_dir = testdir
    data_dir = datasets._get_dataset_dir('test', data_dir=expected_dir,
                                         verbose=0)
    assert data_dir == os.path.join(expected_dir, 'test')
    assert os.path.isdir(data_dir) and os.path.exists(data_dir)
    # test that providing the returned data_dir gets us the same thing
    data_dir2 = datasets._get_dataset_dir('test', data_dir=data_dir,
                                          verbose=0)
    assert data_dir == data_dir2
    shutil.rmtree(data_dir)


def test_fetch_datasets(testdir):
    # test different `donors` inputs
    f1 = datasets.fetch_microarray(data_dir=str(testdir), donors=['12876'])
    f2 = datasets.fetch_microarray(data_dir=str(testdir), donors='12876')
    f3 = datasets.fetch_microarray(data_dir=str(testdir), donors='H0351.1009')
    f4 = datasets.fetch_microarray(data_dir=str(testdir), donors=None)

    # don't test this -- it will take a wicked long time
    # f5 = datasets.fetch_microarray(data_dir=str(testdir), donors='all')

    assert f1 == f2 == f3 == f4
    for k in ['microarray', 'annotation', 'pacall', 'probes', 'ontology']:
        assert len(f1.get(k)) == 1

    # check downloading incorrect donor
    with pytest.raises(ValueError):
        datasets.fetch_microarray(donors='notadonor')

    with pytest.raises(ValueError):
        datasets.fetch_microarray(donors=['notadonor'])


def test_fetch_raw_mri(testdir):
    # test different `donors` inputs
    f1 = datasets.fetch_raw_mri(data_dir=str(testdir), donors=['12876'])
    f2 = datasets.fetch_raw_mri(data_dir=str(testdir), donors='12876')
    f3 = datasets.fetch_raw_mri(data_dir=str(testdir), donors='H0351.1009')
    f4 = datasets.fetch_raw_mri(data_dir=str(testdir), donors=None)
    f5 = datasets.fetch_raw_mri(data_dir=str(testdir), donors='all')

    assert f1 == f2 == f3 == f4
    for k in ['t1w', 't2w']:
        assert len(f1.get(k)) == 1
        assert len(f5.get(k)) == 6

    # check downloading incorrect donor
    with pytest.raises(ValueError):
        datasets.fetch_raw_mri(donors='notadonor')

    with pytest.raises(ValueError):
        datasets.fetch_raw_mri(donors=['notadonor'])


@pytest.mark.parametrize('group, expected', [
    ('brain', 2413),
    ('neuron', 2530),
    ('oligodendrocyte', 1769),
    ('synaptome', 1886),
    ('layers', 46)
])
def test_get_gene_group(group, expected):
    assert len(datasets.fetch_gene_group(group)) == expected

    with pytest.raises(ValueError):
        datasets.fetch_gene_group('notagroup')
