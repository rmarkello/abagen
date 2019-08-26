# -*- coding: utf-8 -*-
"""
Tests for abagen.cli module
"""

import os
import pytest

from abagen.cli import run


def test_run_get_parser(atlas, testdir):
    parser = run.get_parser()

    # since we have a positional argument this fails, hard
    with pytest.raises(SystemExit):
        parser.parse_args([])

    # providing a valid arg succeeds though!
    args = parser.parse_args([atlas['image']])
    assert args.atlas == atlas['image']

    # providing an invalid arg fails again
    with pytest.raises(SystemExit):
        parser.parse_args(['notanatlas.nii.gz'])


def test_run_main(atlas, testdir):
    run.main([
        '--data-dir', testdir,
        '--donors', '12876', '15496',
        '--output-file', os.path.join(str(testdir), 'abagen_expression.csv'),
        atlas['image']
    ])

    assert os.path.exists(os.path.join(str(testdir), 'abagen_expression.csv'))
