# -*- coding: utf-8 -*-

from argparse import ArgumentParser
from pathlib import Path
import sys

from ..allen import get_expression_data


def get_parser():
    from ..info import __version__

    verstr = 'abagen v{}'.format(__version__)
    parser = ArgumentParser(description='')

    parser.add_argument('atlas', action='store', type=Path,
                        help='An image in MNI space, where each parcel in the'
                             'image is identified by a unique integer ID.')

    parser.add_argument('--version', action='version', version=verstr)

    g_data = parser.add_argument_group('Options to specify what AHBA data to '
                                       'use')
    g_data.add_argument('--donors', action='store', nargs='+', default='all',
                        help='List of donors to use as sources of expression '
                             'data. Specified IDs can be either donor numbers '
                             '(i.e., 9861, 10021) or UIDs (i.e., H0351.2001). '
                             'If not specified all available donors will be '
                             'used.')
    g_data.add_argument('--data-dir', '--data_dir', action='store', type=Path,
                        help='Directory where expression data should be '
                             'downloaded to (if it does not already exist) / '
                             'loaded from. If not specified this will check '
                             'the environmental variable ABAGEN_DATA, the '
                             '$HOME/abagen-data directory, and the current '
                             'working directory. If data does not already '
                             'exist at one of those locations then it will be '
                             'downloaded to the first of these location that '
                             'exists and for which write access is enabled.')

    return parser


def main():
    opts = get_parser().parse_args()

    expression = get_expression_data(**opts)

    if opts.stdout:  # WHY?!
        expression.to_csv(sys.stdout)
        return
