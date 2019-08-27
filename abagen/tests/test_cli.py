# -*- coding: utf-8 -*-
"""
Tests for abagen.cli module
"""

import os
from pkg_resources import resource_filename
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


def test_run_main(capsys, atlas, testdir):
    outputfile = os.path.join(str(testdir), 'abagen_expression.csv')

    # check basic usage
    run.main([
        '--data-dir', testdir,
        '--donors', '12876', '15496',
        '--output-file', outputfile,
        atlas['image']
    ])
    assert os.path.exists(outputfile)

    # check that save donors/counts outputs desired files
    run.main([
        '--data-dir', testdir,
        '--donors', '12876', '15496',
        '--output-file', outputfile,
        '--save-donors', '--save-counts',
        atlas['image']
    ])
    for rep in ['_counts.csv', '_12876.csv', '_15496.csv']:
        assert os.path.exists(outputfile.replace('.csv', rep))

    # check stdout (BLARGH)
    run.main([
        '--data-dir', testdir,
        '--donors', '12876', '15496',
        '--output-file', outputfile,
        '--stdout',
        atlas['image']
    ])
    stdout = capsys.readouterr().out
    with open(outputfile, 'r') as src:
        assert stdout == src.read()


def test_exec_run_fail():
    executable = resource_filename('abagen', 'cli/run.py')

    # need to set this otherwise it won't fail
    __name__ = '__main__'  # noqa
    with pytest.raises(RuntimeError):
        with open(executable, 'r') as src:
            exec(src.read())
