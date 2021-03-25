# -*- coding: utf-8 -*-
"""
Tests for abagen.cli module
"""

import os
from pkg_resources import resource_filename
import pytest

from abagen import __version__ as version
from abagen.cli import run


def test_run_get_parser(capsys, atlas, datadir):
    parser = run.get_parser()

    # since we have a positional argument this fails, hard
    with pytest.raises(SystemExit):
        parser.parse_args([])
    assert "following arguments are required: atlas" in capsys.readouterr().err

    # providing the positional succeeds!
    args = parser.parse_args([atlas['image']])
    assert os.path.normcase(args.atlas) == os.path.normcase(atlas['image'])

    # some data directories/files need to exist!
    with pytest.raises(SystemExit):
        parser.parse_args(['notanatlas.nii.gz'])
    assert "does not exist" in capsys.readouterr().err
    with pytest.raises(SystemExit):
        parser.parse_args(['--data-dir', 'notadir', atlas['image']])
    assert "does not exist" in capsys.readouterr().err
    with pytest.raises(SystemExit):
        parser.parse_args(['--atlas-info', 'notafile', atlas['image']])
    assert "does not exist" in capsys.readouterr().err

    # does version print correctly?
    with pytest.raises(SystemExit):
        parser.parse_args(['--version'])
    assert 'abagen {}'.format(version) == capsys.readouterr().out.strip()

    # arguments with invalid choices (probe_selection) raise errors
    for arg in ('--probe-selection', '--gene-norm',
                '--sample-norm', '--lr-mirror'):
        with pytest.raises(SystemExit):
            parser.parse_args([arg, 'notamethod', atlas['image']])
        assert "invalid choice: 'notamethod'" in capsys.readouterr().err

    # just test every option
    args = parser.parse_args([
        '-v', '-v', '-v',
        '--debug',
        '--atlas-info', atlas['info'],
        '--donors', '12876', '15496',
        '--data-dir', datadir,
        '--missing', 'centroids',
        '--tol', '5',
        '--ibf-threshold', '0.6',
        '--region-agg', 'donors',
        '--agg-metric', 'median',
        '--probe-selection', 'average',
        '--norm-all',
        '--norm-structures',
        '--lr_mirror', 'bidirectional',
        '--gene-norm', 'srs',
        '--sample-norm', 'srs',
        '--no-reannotated', '--no-corrected-mni',
        '--output-file', 'test.csv',
        '--save-counts', '--save-donors',
        atlas['image']
    ])


def test_run_main(capsys, atlas, datadir):
    outputfile = os.path.join(datadir, 'abagen_expression.csv')

    # check basic usage
    run.main([
        '--data-dir', datadir,
        '--donors', '12876', '15496',
        '--output-file', outputfile,
        atlas['image']
    ])
    assert os.path.exists(outputfile)

    # check that save donors/counts outputs desired files
    run.main([
        '--data-dir', datadir,
        '--donors', '12876', '15496',
        '--output-file', outputfile,
        '--save-donors', '--save-counts',
        atlas['image']
    ])
    for rep in ['_counts.csv', '_12876.csv', '_15496.csv']:
        assert os.path.exists(outputfile.replace('.csv', rep))

    # check stdout (BLARGH)
    run.main([
        '--data-dir', datadir,
        '--donors', '12876', '15496',
        '--output-file', outputfile,
        '--stdout',
        atlas['image']
    ])
    stdout = capsys.readouterr().out.replace('\r\n', '\n').replace('\r', '\n')
    with open(outputfile, 'r') as src:
        data = src.read().replace('\r\n', '\n').replace('\r', '\n')
        assert stdout == data


def test_exec_run_fail():
    executable = resource_filename('abagen', 'cli/run.py')

    # need to set this otherwise it won't fail
    __name__ = '__main__'  # noqa
    with pytest.raises(RuntimeError):
        with open(executable, 'r') as src:
            exec(src.read())
