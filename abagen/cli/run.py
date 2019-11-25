# -*- coding: utf-8 -*-

import argparse
import logging
import os
from pathlib import Path
import sys

lgr = logging.getLogger('abagen')


def _resolve_path(path):
    """ Helper function for get_parser() to resolve paths
    """

    if path is not None:
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
        if not os.path.exists(_resolve_path(values)):
            parser.error('Provided value for {} does not exist: {}'
                         .format(option_string, values))
        setattr(namespace, self.dest, values)


def get_parser():
    """ Gets command-line arguments for primary get_expression_data workflow
    """

    from .. import __version__
    from ..correct import NORMALIZATION_METHODS
    from ..probes import SELECTION_METHODS

    verstr = 'abagen {}'.format(__version__)
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
Assigns microarray expression data to ROIs defined in the specific `atlas`

This command aims to provide a workflow for generating pre-processed microarray
expression data from the Allen Human Brain Atlas for arbitrary atlas
designations. First, some basic filtering of genetic probes is performed,
including:

    1. Intensity-based filtering of microarray probes to remove probes that do
       not exceed a certain level of background noise (specified via the
       `ibf_threshold` parameter), and
    2. Selection of a single, representative probe (or collapsing across
       probes) for each gene, specified via the `probe_selection` parameter.

Tissue samples are then matched to parcels in the defined `atlas` for each
donor. If `atlas_info` is provided then this matching is constrained by both
hemisphere and tissue class designation (e.g., cortical samples from the left
hemisphere are only matched to ROIs in the left cortex, subcortical samples
from the right hemisphere are only matched to ROIs in the left subcortex); see
the `atlas_info` parameter description for more information.

Matching of microarray samples to parcels in `atlas` is done via a multi-step
process:

    1. Determine if the sample falls directly within a parcel,
    2. Check to see if there are nearby parcels by slowly expanding the search
       space to include nearby voxels, up to a specified distance (specified
       via the `tolerance` parameter),
    3. If there are multiple nearby parcels, the sample is assigned to the
       closest parcel, as determined by the parcel centroid.

If at any step a sample can be assigned to a parcel the matching process is
terminated. If multiple sample are assigned to the same parcel they are
aggregated with the metric specified via the `agg_metric` parameter. More
control over the sample matching can be obtained by setting the `inexact`
parameter; see the parameter description for more information.

Once all samples have been matched to parcels for all supplied donors, the
microarray expression data are optionally normalized via the provided
`sample_norm` and `gene_norm` functions before being combined within parcels
and across donors via the supplied `agg_metric`.
"""
    )

    parser.add_argument('atlas', action=CheckExists, type=_resolve_path,
                        help='An image in MNI space, where each parcel in the '
                             'image is identified by a unique integer ID.')

    # because I like consistency in capitalization and punctuation...
    for act in parser._actions:
        if isinstance(act, argparse._HelpAction):
            act.help = act.help.capitalize() + '.'
            break

    parser.add_argument('--version', action='version', version=verstr,
                        help='Show program\'s version number and exit.')
    parser.add_argument('-v', '--verbose', action='count', default=1,
                        help='Increase verbosity of status messages to '
                             'display during workflow.')
    parser.add_argument('-q', '--quiet', action='store_true', default=False,
                        help='Suppress all status messages during workflow.')
    parser.add_argument('--debug', action='store_true', help=argparse.SUPPRESS)

    a_data = parser.add_argument_group('Options to specify information about '
                                       'the atlas used')
    a_data.add_argument('--atlas_info', '--atlas-info', action=CheckExists,
                        type=_resolve_path, default=None, metavar='PATH',
                        help='Filepath to CSV files containing information '
                             'about `atlas`. The CSV file must have at least '
                             'columns ["id", "hemisphere", "structure"] which '
                             'contain information mapping the atlas IDs to '
                             'hemispheres (i.e, "L", "R") and broad '
                             'structural groups (i.e., "cortex", "subcortex", '
                             '"cerebellum", "brainstem").')

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

    w_data = parser.add_argument_group('Options to specify processing options')
    w_data.add_argument('--inexact', dest='exact', action='store_false',
                        default=True,
                        help='Whether to use inexact matching of donor tissue '
                             'samples to parcels in `atlas`. By default, the '
                             'workflow will match tissue samples to parcels '
                             'within `tolerance` mm of the sample; any '
                             'samples that are beyond `tolerance` mm of a '
                             'parcel will be discarded, which may result in '
                             'some parcels having no assigned sample / '
                             'expression data. If --inexact, the matching '
                             'procedure will be performed and followed by a '
                             'check for parcels with no assigned samples; any '
                             'such parcels will be matched to the nearest '
                             'sample (defined as the sample with the closest '
                             'Euclidean distance to the parcel centroid). '
                             'Default: False')
    w_data.add_argument('--tol', '--tolerance', dest='tolerance',
                        action='store', type=float, default=2,
                        help='Distance (in mm) that a sample can be from a '
                             'parcel for it to be matched to that parcel. '
                             'This is only considered if the sample is not '
                             'directly within a parcel. Default: 2')
    w_data.add_argument('--ibf_threshold', '--ibf-threshold', action='store',
                        default=0.5, metavar='THRESHOLD',
                        help='Threshold for intensity-based filtering of '
                             'probes. This number should specify the ratio of '
                             'samples, across all supplied donors, for which '
                             'a probe must have signal above background noise '
                             'in order to be retained. Default: 0.5')
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
    w_data.add_argument('--probe_selection', '--probe-selection',
                        action='store', default='diff_stability',
                        metavar='METHOD', choices=sorted(SELECTION_METHODS),
                        help='Selection method for subsetting (or collapsing '
                             'across) probes that index the same gene. Must '
                             'be one of {"average", "mean", "max_intensity", '
                             '"max_variance", "pc_loading", "corr_variance", '
                             '"corr_intensity", "diff_stability"}. Default: '
                             '"diff_stability"')
    w_data.add_argument('--lr_mirror', '--lr-mirror', action='store_true',
                        help='Whether to mirror microarray expression samples '
                             'across hemispheres to increase spatial '
                             'coverage. This will duplicate samples across '
                             'both hemispheres (i.e., L->R and R->L), '
                             'approximately doubling the number of available '
                             'samples. Default: False (i.e., no mirroring)')
    w_data.add_argument('--gene_norm', '--gene-norm', action='store',
                        default='srs', metavar='METHOD', type=_resolve_none,
                        choices=sorted(NORMALIZATION_METHODS) + ['None'],
                        help='Method by which to normalize microarray '
                             'expression values for each donor prior to '
                             'collapsing across donors. Expression values are '
                             'normalized separately for each gene for each '
                             'donor across all expression samples. If None is '
                             'specified then no normalization is performed. '
                             'Default: "srs"')
    w_data.add_argument('--sample_norm', '--sample-norm', action='store',
                        default='srs', metavar='METHOD', type=_resolve_none,
                        choices=sorted(NORMALIZATION_METHODS) + ['None'],
                        help='Method by which to normalize microarray '
                             'expression values for each sample prior to '
                             'collapsing into regions in `atlas`. Expression '
                             'values are normalized separately for each '
                             'sample and donor across genes. If None is '
                             'specified then no normalization is performed. '
                             'Default: "srs"')

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
    from ..datasets import WELL_KNOWN_IDS as donors

    opts = get_parser().parse_args(args)

    # quiet overrides any verbosity setting
    if opts.quiet:
        opts.verbose = 0

    # debugging is fun
    if opts.debug:
        print(opts)
        return

    # run the workflow
    expression = get_expression_data(atlas=opts.atlas,
                                     atlas_info=opts.atlas_info,
                                     exact=opts.exact,
                                     tolerance=opts.tolerance,
                                     region_agg=opts.region_agg,
                                     agg_metric=opts.agg_metric,
                                     ibf_threshold=opts.ibf_threshold,
                                     probe_selection=opts.probe_selection,
                                     lr_mirror=opts.lr_mirror,
                                     gene_norm=opts.gene_norm,
                                     sample_norm=opts.sample_norm,
                                     corrected_mni=opts.corrected_mni,
                                     reannotated=opts.reannotated,
                                     return_counts=opts.save_counts,
                                     return_donors=opts.save_donors,
                                     donors=opts.donors,
                                     data_dir=opts.data_dir,
                                     verbose=opts.verbose)

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
        lgr.info('Saving samples counts to {}'.format(counts_fname))
        counts.to_csv(counts_fname)

    # determine how best to save expression output files
    if opts.save_donors:
        if opts.donors == 'all':
            donors = list(donors.value_set('subj'))
        else:
            donors = [donors[f] for f in opts.donors]

        # save each donor dataframe as a separate file
        for donor, exp in zip(donors, expression):
            exp_fname = os.path.join(output_path,
                                     fname_pref + '_{}.csv'.format(donor))
            lgr.info('Saving donor {} info to {}'.format(donor, exp_fname))
            exp.to_csv(exp_fname)
    else:
        expression.to_csv(opts.output_file)


if __name__ == '__main__':
    raise RuntimeError('abagen/cli/run.py should not be run directly.\nPlease '
                       '`pip install` abagen and use the `abagen` command.')
