# -*- coding: utf-8 -*-
"""
Fixtures for all abagen tests
"""

import os
import pytest

from abagen.datasets import fetch_desikan_killiany, fetch_microarray


@pytest.fixture(scope='session')
def datadir(tmp_path_factory):
    dd = os.environ.get('ABAGEN_DATA')
    if dd is None:
        dd = str(tmp_path_factory.mktemp('abagen-data'))
    return dd


@pytest.fixture(scope='session')
def testfiles(datadir):
    return fetch_microarray(data_dir=datadir, donors=['12876', '15496'],
                            n_proc=2)


@pytest.fixture(scope='session')
def atlas():
    return fetch_desikan_killiany()
