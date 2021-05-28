# -*- coding: utf-8 -*-

import argparse
import logging
import os
from pathlib import Path
from typing import Iterable
import sys

LGR = logging.getLogger('abagen')


def isiterable(val):
    """ Helper function to check whether value is iterable (but not string)
    """
    return isinstance(val, Iterable) and not isinstance(val, str)


def _resolve_path(path):
    """ Helper function for get_parser() to resolve paths
    """

    if path is not None:
        if isiterable(path):
            return [_resolve_path(p) for p in path]
        try:
            return str(Path(path).expanduser().resolve())
        except FileNotFoundError:
            return os.path.abspath(os.path.expanduser(path))


def _resolve_none(inp):
    """ Helper function to allow 'None' as input from argparse
    """

    if inp == "None":
        return
    return inp


class CheckExists(argparse.Action):
    """ Helper class to check that provided paths exist
    """
    def __call__(self, parser, namespace, values, option_string=None):
        values = self.type(values)
        missing = False
        if isiterable(values):
            missing = any(not os.path.exists(val) for val in values)
            if len(values) == 1:
                values = values[0]
        else:
            missing = not os.path.exists(values)
        if missing:
            parser.error('Provided value for {} does not exist: {}'
                         .format(option_string, values))
        setattr(namespace, self.dest, values)


def get_parser():
    """ Gets command-line arguments for primary get_expression_data workflow
    """

    from .. import __version__
    from ..correct import NORMALIZATION_METHODS
    from ..probes_ import SELECTION_METHODS

    verstr = 'abagen {}'.format(__version__)
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
Assigns microarray expression data to ROIs defined in the specified `atlas`

This command aims to provide a workflow for generating pre-processed microarray
expression data from the Allen Human Brain Atlas for arbitrary atlas
designations. First, some basic filtering of genetic probes is performed,
including:

    1. Intensity-based filtering of microarray probes to remove probes that do
       not exceed a certain level of background noise (specified via the
       `--ibf_threshold` parameter),
    2. Selection of a single, representative probe (or collapsing across
       probes) for each gene, specified via the `--probe_selection` parameter
       (and influenced by the `--donor_probes` parameter), and
    3. Optional mirroring of the tissue samples across the left/right
       hemisphere boundary, as specified via the `--lr_mirror` parameter
       (turned off by default).

Tissue samples are then matched to parcels in the defined `atlas` for each
donor. If `--atlas_info` is provided then this matching is constrained by both
hemisphere and tissue class designation (e.g., cortical samples from the left
hemisphere are only matched to ROIs in the left cortex, subcortical samples
from the right hemisphere are only matched to ROIs in the left subcortex); see
the `atlas_info` parameter description for more information.

Matching of microarray samples to parcels in `atlas` is done via a multi-step
process:

    1. Determine if the sample falls directly within a parcel,
    2. Check to see if there are nearby parcels by slowly expanding the search
       space to include nearby voxels, up to a specified distance (specified
       via the `--tolerance` parameter),
    3. If there are multiple nearby parcels, the sample is assigned to the
       closest parcel, as determined by the parcel centroid.

If at any step a sample can be assigned to a parcel the matching process is
terminated. When the provided atlas is not volumetric (i.e., is surface-based)
the samples are simply matched to the nearest vertex, and `--tolerance` is used
as a standard deviation threshold. More control over the sample matching can be
obtained by setting the `--missing` parameter.

Once all samples have been matched to parcels for all supplied donors, the
microarray expression data are optionally normalized via the provided
`--sample_norm` and `--gene_norm` functions (which are influenced by the
`--norm_matched` and `--norm_structures` parameters) before being aggregated
across donors via the supplied `--region_agg` and `--agg_metric` parameters.
"""
    )

    parser.add_argument('atlas', action=CheckExists, type=_resolve_path,
                        nargs='+',
                        help='A NIFTI image in MNI152 space or two GIFTI '
                             'images in fsaverage5 space, where each parcel '
                             'is identified by a unique integer ID.')

    # because I like consistency in capitalization and punctuation...
    for act in parser._actions:
        if isinstance(act, argparse._HelpAction):
            act.help = act.help.capitalize() + '.'
            break

    parser.add_argument('--version', action='version', version=verstr,
                        help='Show program version and exit.')
    parser.add_argument('-v', '--verbose', action='count', default=0,
                        help='Increase verbosity of status messages to '
                             'display during workflow.')
    parser.add_argument('--debug', action='store_true', help=argparse.SUPPRESS)

    a_data = parser.add_argument_group('Options to specify information about '
                                       'the atlas used')
    a_data.add_argument('--atlas_info', '--atlas-info', action=CheckExists,
                        type=_resolve_path, default=None, metavar='PATH',
                        help='Filepath to CSV file containing information '
                             'about `atlas`. The CSV file must have at least '
                             'columns ["id", "hemisphere", "structure"] which '
                             'contain information mapping the atlas IDs to '
                             'hemispheres (i.e, "L", "R", or "B") and broad '
                             'structural groups (i.e., "cortex", "subcortex/'
                             'brainstem", "cerebellum"). If provided, this '
                             'will constrain matching of tissue samples to '
                             'regions in `atlas`. If the supplied `atlas` is '
                             'a pair of GIFTI files with valid label tables '
                             'this information will be intuited.')

    g_data = parser.add_argument_group('Options to specify which AHBA data to '
                                       'use during processing')
    g_data.add_argument('--donors', action='store', nargs='+',
                        default='all', metavar='DONOR_ID',
                        help='List of donors to use as sources of expression '
                             'data. Specified IDs can be either donor numbers '
                             '(i.e., 9861, 10021) or UIDs (i.e., H0351.2001). '
                             'Can specify "all" to use all available donors. '
                             'Default: "all"')
    g_data.add_argument('--data_dir', '--data-dir', action=CheckExists,
                        type=_resolve_path, metavar='PATH',
                        help='Directory where expression data should be '
                             'downloaded to (if it does not already exist) / '
                             'loaded from. If not specified this will check '
                             'the environmental variable $ABAGEN_DATA, the '
                             '$HOME/abagen-data directory, and the current '
                             'working directory. If data does not already '
                             'exist at one of those locations then it will be '
                             'downloaded to the first of these location that '
                             'exists and for which write access is enabled.')
    g_data.add_argument('--n_proc', '--n-proc', action='store', type=int,
                        default=1,
                        help='Number of processors to use to download AHBA '
                             'data. Can paralellize up to six times if all '
                             'donors are requested. Default: 1')

    w_data = parser.add_argument_group('Options to specify processing options')
    w_data.add_argument('--ibf_threshold', '--ibf-threshold', action='store',
                        default=0.5, metavar='THRESHOLD',
                        help='Threshold for intensity-based filtering of '
                             'probes. This number should specify the ratio of '
                             'samples, across all supplied donors, for which '
                             'a probe must have signal above background noise '
                             'in order to be retained. Default: 0.5')
    w_data.add_argument('--probe_selection', '--probe-selection',
                        action='store', default='diff_stability',
                        metavar='METHOD', choices=sorted(SELECTION_METHODS),
                        help='Selection method for subsetting (or collapsing '
                             'across) probes that index the same gene. Must '
                             'be one of {"average", "mean", "max_intensity", '
                             '"max_variance", "pc_loading", "corr_variance", '
                             '"corr_intensity", "diff_stability", "rnaseq"}. '
                             'Default: "diff_stability"')
    w_data.add_argument('--lr_mirror', '--lr-mirror', metavar='METHOD',
                        type=_resolve_none, default=None, choices=(
                            None, 'bidirectional', 'leftright', 'rightleft'),
                        help='Whether to mirror microarray expression samples '
                             'across hemispheres to increase spatial coverage.'
                             ' Using "bidirectional" will mirror samples '
                             'across both hemispheres, "leftright" will '
                             'mirror samples in the left hemisphere to the '
                             'right, and "rightleft" will mirror the right to '
                             'the left. Default: None')
    w_data.add_argument('--sim_threshold', '--sim-threshold',
                        type=_resolve_none, default=None, metavar='THRESHOLD',
                        help='Threshold for inter-areal similarity filtering. '
                             'Samples are correlated across probes and those '
                             'samples with a total correlation less than the '
                             'the provided threshold s.d. below the mean '
                             'across samples are excluded from futher'
                             'analysis. If not specified no filtering is '
                             'performed. Default: None')
    w_data.add_argument('--missing', dest='missing', metavar='METHOD',
                        type=_resolve_none, default=None, choices=(
                            None, 'centroids', 'interpolate'),
                        help='How to handle regions in `atlas` that are not '
                             'assigned any tissue samples. If "centroids", '
                             'any empty regions will be assigned the '
                             'expression value of the nearest tissue sample '
                             '(defined as the sample with the closest '
                             'Euclidean distance to the parcel centroid). If '
                             '"interpolate", expression values will be '
                             'interpolated in the empty regions by assigning '
                             'every node in the region the expression of the '
                             'nearest sample and taking a weighted (inverse '
                             'distance) average. If not specified empty '
                             'regions will be returned with expression values '
                             'of NaN. Default: None')
    w_data.add_argument('--tol', '--tolerance', dest='tolerance',
                        action='store', type=float, default=2,
                        help='Distance (in mm) that a sample can be from a '
                             'parcel for it to be matched to that parcel. If '
                             '`atlas` is GIFTI files then this measure is a '
                             'standard deviation threshold (i.e., samples '
                             'greater than `tolerance` SDs away from the mean '
                             'matched distance are ignored). Default: 2')
    w_data.add_argument('--sample_norm', '--sample-norm', action='store',
                        default='srs', metavar='METHOD', type=_resolve_none,
                        choices=sorted(NORMALIZATION_METHODS) + ['None', None],
                        help='Method by which to normalize microarray '
                             'expression values for each sample prior to '
                             'collapsing into regions in `atlas`. Expression '
                             'values are normalized separately for each '
                             'sample and donor across genes. If None is '
                             'specified then no normalization is performed. '
                             'Default: "srs"')
    w_data.add_argument('--gene_norm', '--gene-norm', action='store',
                        default='srs', metavar='METHOD', type=_resolve_none,
                        choices=sorted(NORMALIZATION_METHODS) + ['None', None],
                        help='Method by which to normalize microarray '
                             'expression values for each donor prior to '
                             'collapsing across donors. Expression values are '
                             'normalized separately for each gene for each '
                             'donor across all expression samples. If None is '
                             'specified then no normalization is performed. '
                             'Default: "srs"')
    w_data.add_argument('--norm_all', '--norm-all', dest='norm_matched',
                        action='store_false', default=True,
                        help='Whether to perform gene normalization '
                             '(`gene_norm`) across all available samples '
                             'instead of only across samples that were '
                             'matched to regions in `atlas`. If `atlas` is '
                             'very small (i.e., only a few regions of '
                             'interest) using `--norm_all` is suggested.')
    w_data.add_argument('--norm_structures', '--norm-structures',
                        action='store_true', default=False,
                        help='Whether to perform gene normalization '
                             '(`gene_norm`) within structural classes (i.e., '
                             '"cortex", "subcortex/brainstem", "cerebellum") '
                             'instead of across all available samples.')
    w_data.add_argument('--region_agg', '--region-agg', action='store',
                        default='donors', metavar='METHOD',
                        choices=['donors', 'samples'],
                        help='When multiple samples are identified as '
                             'belonging to a region in `atlas` this '
                             'determines how they are aggegated. If '
                             '\'samples\', expression data from all samples '
                             'for all donors assigned to a given region are '
                             'combined. If \'donors\', expression values for '
                             'all samples assigned to a given region are '
                             'combined independently for each donor before '
                             'being combined across donors. See `agg_metric` '
                             'for mechanism by which samples are combined. '
                             'Default: \'donors\'')
    w_data.add_argument('--agg_metric', '--agg-metric', action='store',
                        default='mean', metavar='METHOD',
                        choices=['mean', 'median'],
                        help='Mechanism by which to (1) reduce expression '
                             'data of multiple samples in the same `atlas` '
                             'region, and (2) reduce donor-level expression '
                             'data into a single "group" expression '
                             'dataframe. Must be one of {"mean", "median"}. '
                             'Default: "mean"')

    p_data = parser.add_argument_group('Options to modify the AHBA data used')
    p_data.add_argument('--no-reannotated', '--no_reannotated',
                        dest='reannotated', action='store_false', default=True,
                        help='Whether to use the original probe information '
                             'from the AHBA dataset instead of the '
                             'reannotated probe information from '
                             'Arnatkevic̆iūtė et al., 2019. Using reannotated '
                             'probe information discards probes that could '
                             'not be reliably matched to genes. Default: '
                             'False (i.e., use reannotations)')
    p_data.add_argument('--no-corrected-mni', '--no_corrected_mni',
                        dest='corrected_mni', action='store_false',
                        default=True,
                        help='Whether to use the original MNI coordinates '
                             'provided with the AHBA data instead of the '
                             '"corrected" MNI coordinates shipped with the '
                             '`alleninf` package when matching tissue samples '
                             'to anatomical regions. Default: False (i.e., '
                             'use corrected coordinates)')

    o_data = parser.add_argument_group('Options to modify how data are output')
    o_data.add_argument('--stdout', action='store_true',
                        help='Generated region x gene dataframes will be '
                             'printed to stdout for piping to other things. '
                             'You should REALLY consider just using --output-'
                             'file instead and working with the generated '
                             'CSV file(s). Incompatible with `--save-counts` '
                             'and `--save-donors` (i.e., this will override '
                             'those options). Default: False')
    o_data.add_argument('--output-file', '--output_file', action='store',
                        type=_resolve_path, metavar='PATH',
                        default='abagen_expression.csv',
                        help='Path to desired output file. The generated '
                             'region x gene dataframe will be saved here. '
                             'Default: $PWD/abagen_expression.csv')
    o_data.add_argument('--save-counts', '--save_counts', action='store_true',
                        help='Whether to save dataframe containing number of '
                             'samples from each donor that were assigned '
                             'to each region in `atlas`. If specified, will '
                             'be saved to the path specified by '
                             '`output-file`, appending "counts" to the end of '
                             'the filename. Default: False')
    o_data.add_argument('--save-donors', '--save_donors', action='store_true',
                        help='Whether to save donor-level expression '
                             'dataframes instead of aggregating expression '
                             'across donors with provided `agg_metric`. If '
                             'specified, dataframes will be saved to path '
                             'specified by `output-file`, appending donor IDs '
                             'to the end of the filename. Default: False')

    return parser


def main(args=None):
    """ Runs primary get_expression_data workflow
    """

    from ..allen import get_expression_data

    opts = get_parser().parse_args(args)

    # debugging is fun
    if opts.debug:
        print(opts)
        return

    # run the workflow
    expression = get_expression_data(atlas=opts.atlas,
                                     atlas_info=opts.atlas_info,
                                     ibf_threshold=opts.ibf_threshold,
                                     probe_selection=opts.probe_selection,
                                     sim_threshold=opts.sim_threshold,
                                     lr_mirror=opts.lr_mirror,
                                     missing=opts.missing,
                                     tolerance=opts.tolerance,
                                     sample_norm=opts.sample_norm,
                                     gene_norm=opts.gene_norm,
                                     norm_matched=opts.norm_matched,
                                     norm_structures=opts.norm_structures,
                                     region_agg=opts.region_agg,
                                     agg_metric=opts.agg_metric,
                                     corrected_mni=opts.corrected_mni,
                                     reannotated=opts.reannotated,
                                     return_counts=opts.save_counts,
                                     return_donors=opts.save_donors,
                                     donors=opts.donors,
                                     data_dir=opts.data_dir,
                                     verbose=opts.verbose,
                                     n_proc=opts.n_proc)

    output_path = os.path.dirname(opts.output_file)
    fname_pref = os.path.splitext(os.path.basename(opts.output_file))[0]

    # WHY?!?
    if opts.stdout and not (opts.save_counts or opts.save_donors):
        expression.to_csv(sys.stdout)
        return

    # expand the tuple, if needed
    if opts.save_counts:
        expression, counts = expression
        counts_fname = os.path.join(output_path, fname_pref + '_counts.csv')
        LGR.info('Saving samples counts to {}'.format(counts_fname))
        counts.to_csv(counts_fname)

    # determine how best to save expression output files
    if opts.save_donors:
        # save each donor dataframe as a separate file
        for donor, exp in expression.items():
            exp_fname = os.path.join(output_path,
                                     fname_pref + '_{}.csv'.format(donor))
            LGR.info('Saving donor {} info to {}'.format(donor, exp_fname))
            exp.to_csv(exp_fname)
    else:
        expression.to_csv(opts.output_file)


if __name__ == '__main__':
    raise RuntimeError('abagen/cli/run.py should not be run directly.\nPlease '
                       '`pip install` abagen and use the `abagen` command.')
