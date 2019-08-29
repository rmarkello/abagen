import pytest
from abagen.datasets import fetch_desikan_killiany, fetch_microarray


@pytest.fixture(scope='session')
def testdir(tmpdir_factory):
    data_dir = tmpdir_factory.mktemp('data')
    return str(data_dir)


@pytest.fixture(scope='session')
def testfiles(testdir):
    files = fetch_microarray(data_dir=testdir,
                             donors=['12876', '15496'],
                             convert=True)
    return files


@pytest.fixture(scope='session')
def atlas():
    return fetch_desikan_killiany()
