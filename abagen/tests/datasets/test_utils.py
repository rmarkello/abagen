# -*- coding: utf-8 -*-
"""
Tests for abagen.datasets.utils module
"""

import os
import shutil

from abagen.datasets import utils


def test_get_dataset_dir(datadir):
    old = os.environ.pop('ABAGEN_DATA', None)

    # check that data dir defaults to $HOME/abagen-data assuming no env var
    expected_base = os.path.join(os.path.expanduser('~'), 'abagen-data')
    data_dir = utils._get_dataset_dir('test', verbose=0)
    assert data_dir == os.path.join(expected_base, 'test')
    assert os.path.isdir(data_dir) and os.path.exists(data_dir)
    shutil.rmtree(data_dir)

    # if env var is set, we should default to that
    expected_base = os.path.join(datadir, 'test-abagen-data')
    os.environ['ABAGEN_DATA'] = expected_base
    data_dir = utils._get_dataset_dir('test', verbose=0)
    assert data_dir == os.path.join(expected_base, 'test')
    assert os.path.isdir(data_dir) and os.path.exists(data_dir)
    shutil.rmtree(expected_base)

    # test explicitly setting the data_dir
    expected_dir = datadir
    data_dir = utils._get_dataset_dir('test', data_dir=expected_dir,
                                      verbose=0)
    assert data_dir == os.path.join(expected_dir, 'test')
    assert os.path.isdir(data_dir) and os.path.exists(data_dir)
    # test that providing the returned data_dir gets us the same thing
    data_dir2 = utils._get_dataset_dir('test', data_dir=data_dir,
                                       verbose=0)
    assert data_dir == data_dir2
    shutil.rmtree(data_dir)

    if old is not None:
        os.environ['ABAGEN_DATA'] = old
