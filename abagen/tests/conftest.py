import pytest
from abagen.datasets import fetch_microarray


@pytest.fixture(scope='session')
def testdir(tmpdir_factory):
    data_dir = tmpdir_factory.mktemp('data')
    return data_dir


@pytest.fixture(scope='session')
def testfiles(testdir):
    files = fetch_microarray(data_dir=str(testdir),
                             donors=['12876', '15496'],
                             convert=True)
    return files
