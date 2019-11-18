# -*- coding: utf-8 -*-
"""
Fixtures for all abagen tests
"""

import os
import pytest

from abagen.datasets import fetch_desikan_killiany, fetch_microarray


@pytest.fixture(scope='session')
def datadir():
    dd = os.environ.get('ABAGEN_DATA')
    if dd is None:
        return os.path.join(os.environ['HOME'], 'abagen-data')
    return dd


@pytest.fixture(scope='session')
def testfiles():
    return fetch_microarray(donors=['12876', '15496'])


@pytest.fixture(scope='session')
def atlas():
    return fetch_desikan_killiany()
