# -*- coding: utf-8 -*-
"""
Fixtures for all abagen tests
"""

import os
import pytest
from pathlib import Path
from abagen.datasets import (fetch_desikan_killiany,
                             fetch_microarray,
                             fetch_raw_mri,
                             fetch_rnaseq)


@pytest.fixture(scope='session')
def datadir(tmp_path_factory):
    dd = os.environ.get('ABAGEN_DATA')
    if dd is None:
        dd = str(tmp_path_factory.mktemp('abagen-data'))
    return str(Path(dd).expanduser())


@pytest.fixture(scope='session')
def testfiles(datadir):
    return fetch_microarray(data_dir=datadir, donors=['12876', '15496'],
                            n_proc=2)


@pytest.fixture(scope='session')
def rnafiles(datadir):
    return fetch_rnaseq(data_dir=datadir, donors=['9861'])


@pytest.fixture(scope='session')
def rawmri(datadir):
    return fetch_raw_mri(data_dir=datadir, donors=['12876', '15496'])


@pytest.fixture(scope='session')
def atlas():
    return fetch_desikan_killiany(native=False, surface=False)


@pytest.fixture(scope='session')
def surface():
    return fetch_desikan_killiany(native=False, surface=True)
