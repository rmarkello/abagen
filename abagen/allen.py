# -*- coding: utf-8 -*-
"""
Functions for mapping AHBA microarray dataset to atlases and and parcellations
"""

from functools import reduce
import logging
import warnings

import nibabel as nib
import numpy as np
import pandas as pd

from . import (correct, datasets, images, io, matching, probes_, reporting,
               samples_, utils)
from .utils import first_entry, flatten_dict

logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.WARNING)
LGR = logging.getLogger('abagen')


def get_expression_data(atlas,
                        atlas_info=None,
                        *,
                        ibf_threshold=0.5,
                        probe_selection='diff_stability',
                        donor_probes='aggregate',
                        sim_threshold=None,
                        lr_mirror=None,
                        exact=None, missing=None,
                        tolerance=2,
                        sample_norm='srs',
                        gene_norm='srs',
                        norm_matched=True,
                        norm_structures=False,
                        region_agg='donors',
                        agg_metric='mean',
                        corrected_mni=True,
                        reannotated=True,
                        return_counts=False,
                        return_donors=False,
                        return_report=False,
                        donors='all',
                        data_dir=None,
                        verbose=0,
                        n_proc=1):
    """
    Assigns microarray expression data to ROIs defined in `atlas`

    This function aims to provide a workflow for generating pre-processed,
    microarray expression data from the Allen Human Brain Atlas ([A2]_) for
    abitrary `atlas` designations. First, some basic filtering of genetic
    probes is performed, including:

        1. Intensity-based filtering of microarray probes to remove probes that
           do not exceed a certain level of background noise (specified via the
           `ibf_threshold` parameter),
        2. Selection of a single, representative probe (or collapsing across
           probes) for each gene, specified via the `probe_selection`
           parameter (and influenced by the `donor_probes` parameter), and
        3. Optional mirroring of the tissue samples across the left/right
           hemisphere boundary, as specified via the `lr_mirror` parameter
           (turned off by default).

    Tissue samples are then matched to parcels in the defined `atlas` for each
    donor. If `atlas_info` is provided then this matching is constrained by
    both hemisphere and tissue class designation (e.g., cortical samples from
    the left hemisphere are only matched to ROIs in the left cortex,
    subcortical samples from the right hemisphere are only matched to ROIs in
    the left subcortex); see the `atlas_info` parameter description for more
    information.

    Matching of microarray samples to parcels in `atlas` is done via a multi-
    step process:

        1. Determine if the sample falls directly within a parcel,
        2. Check to see if there are nearby parcels by slowly expanding the
           search space to include nearby voxels, up to a specified distance
           (specified via the `tolerance` parameter),
        3. If there are multiple nearby parcels, the sample is assigned to the
           closest parcel, as determined by the parcel centroid.

    If at any step a sample can be assigned to a parcel the matching process is
    terminated. When the provided atlas is not volumetric (i.e., surface-based)
    the samples are simply matched to the nearest vertex, and `tolerance` is
    used as a standard deviation threshold. More control over the sample
    matching can be obtained by setting the `exact` parameter; see the
    parameter description for more information.

    Once all samples have been matched to parcels for all supplied donors, the
    microarray expression data are optionally normalized via the provided
    `sample_norm` and `gene_norm` functions (which are influenced by the
    `norm_matched` and `norm_structures` parameters) before being aggregated
    across donors via the supplied `region_agg` and `agg_metric` parameters.

    Parameters
    ----------
    atlas : niimg-like object or dict
        A parcellation image in MNI space or a tuple of GIFTI images in
        fsaverage5 space, where each parcel is identified by a unique integer
        ID. Alternatively, a dictionary where keys are donor IDs and values are
        parcellation images (or surfaces) in the native space of each donor.
    atlas_info : os.PathLike or pandas.DataFrame, optional
        Filepath to or pre-loaded dataframe containing information about
        `atlas`. Must have at least columns 'id', 'hemisphere', and 'structure'
        containing information mapping atlas IDs to hemisphere (i.e, "L", "R",
        "B") and broad structural class (i.e., "cortex", "subcortex/brainstem",
        "cerebellum"). If provided, this will constrain matching of tissue
        samples to regions in `atlas`. If `atlas` is a tuple of GIFTI images
        with valid label tables this will be intuited from the data. Default:
        None
    ibf_threshold : [0, 1] float, optional
        Threshold for intensity-based filtering. This number specifies the
        ratio of samples, across all supplied donors, for which a probe must
        have signal significantly greater than background noise in order to be
        retained. Default: 0.5
    probe_selection : str, optional
        Selection method for subsetting (or collapsing across) probes that
        index the same gene. Must be one of 'average', 'max_intensity',
        'max_variance', 'pc_loading', 'corr_variance', 'corr_intensity', or
        'diff_stability', 'rnaseq'; see Notes for more information on different
        options. Default: 'diff_stability'
    donor_probes : {'aggregate', 'independent', 'common'}, optional
        Whether specified `probe_selection` method should be performed with
        microarray data from all donors ('aggregate'), independently for each
        donor ('independent'), or based on the most common selected probe
        across donors ('common'). Not all combinations of `probe_selection`
        and `donor_probes` methods are viable. Default: 'aggregate'
    sim_threshold : (0, inf) float, optional
        Threshold for inter-areal similarity filtering. Samples are correlated
        across probes and those samples with a total correlation less than
        `sim_threshold` standard deviations below the mean across samples are
        excluded from futher analysis. If not specified no filtering is
        performed. Default: None
    lr_mirror : {None, 'bidirectional', 'leftright', 'rightleft'}, optional
        Whether to mirror microarray expression samples across hemispheres to
        increase spatial coverage. Using 'bidirectional' will mirror samples
        across both hemispheres, 'leftright' will mirror samples in the left
        hemisphere to the right, and 'rightleft' will mirror the right to the
        left. Default: None
    missing : {'centroids', 'interpolate', None}, optional
        How to handle regions in `atlas` that are not assigned any tissue
        samples. If 'centroids', any empty regions will be assigned the
        expression value of the nearest tissue sample (defined as the sample
        with the closest Euclidean distance to the parcel centroid). If
        'interpolate', expression values will be interpolated in the empty
        regions by assigning every node in the region the expression of the
        nearest sample and taking a weighted (inverse distance) average. If
        not specified empty regions will be returned with expression values of
        NaN. Default: None
    tolerance : int, optional
        Distance (in mm) that a sample must be from a parcel for it to be
        matched to that parcel. If `atlas` is a tuple of surface files then
        this measure is a standard deviation threshold (i.e., samples greater
        than `tolerance` SDs away from the mean matched distance are ignored).
        Default: 2
    sample_norm : {'rs', 'srs', 'minmax', 'center', 'zscore', None}, optional
        Method by which to normalize microarray expression values for each
        sample. Expression values are normalized separately for each sample and
        donor across all genes; see Notes for more information on different
        methods. If None is specified then no normalization is performed.
        Default: 'srs'
    gene_norm : {'rs', 'srs', 'minmax', 'center', 'zscore', None}, optional
        Method by which to normalize microarray expression values for each
        donor. Expression values are normalized separately for each gene and
        donor across all samples; see Notes for more information on different
        methods. If None is specified then no normalization is performed.
        Default: 'srs'
    norm_matched : bool, optional
        Whether to perform gene normalization (`gene_norm`) across only those
        samples matched to regions in `atlas` instead of all available samples.
        If `atlas` is very small (i.e., only a few regions of interest), using
        `norm_matched=False` is suggested. Default: True
    norm_structures : bool, optional
        Whether to perform gene normalization (`gene_norm`) within structural
        classes (i.e., 'cortex', 'subcortex/brainstem', 'cerebellum') instead
        of across all available samples. Default: False
    region_agg : {'samples', 'donors'}, optional
        When multiple samples are identified as belonging to a region in
        `atlas` this determines how they are aggegated. If 'samples',
        expression data from all samples for all donors assigned to a given
        region are combined. If 'donors', expression values for all samples
        assigned to a given region are combined independently for each donor
        before being combined across donors. See `agg_metric` for mechanism by
        which samples are combined. Default: 'donors'
    agg_metric : {'mean', 'median'} or callable, optional
        Mechanism by which to reduce sample-level expression data into region-
        level expression (see `region_agg`). If a callable, should be able to
        accept an `N`-dimensional input and the `axis` keyword argument and
        return an `N-1`-dimensional output. Default: 'mean'
    corrected_mni : bool, optional
        Whether to use the "corrected" MNI coordinates shipped with the
        `alleninf` package instead of the coordinates provided with the AHBA
        data when matching tissue samples to anatomical regions. Default: True
    reannotated : bool, optional
        Whether to use reannotated probe information provided by [A1]_ instead
        of the default probe information from the AHBA dataset. Using
        reannotated information will discard probes that could not be reliably
        matched to genes. Default: True
    return_counts : bool, optional
        Whether to return dataframe containing information on how many samples
        were assigned to each parcel in `atlas` for each donor. Default: False
    return_donors : bool, optional
        Whether to return donor-level expression arrays instead of aggregating
        expression across donors with provided `agg_metric`. Default: False
    return_report : bool, optional
        Whether to return a string containing longform text describing the
        processing procedures used to generate the `expression` DataFrames
        returned by this function. Default: False
    donors : list, optional
        List of donors to use as sources of expression data. Can be either
        donor numbers or UID. If not specified will use all available donors.
        Note that donors '9861' and '10021' have samples from both left + right
        hemispheres; all other donors have samples from the left hemisphere
        only. Default: 'all'
    data_dir : os.PathLike, optional
        Directory where expression data should be downloaded (if it does not
        already exist) / loaded. If not specified will use the current
        directory. Default: None
    verbose : int, optional
        Specifies verbosity of status messages to display during workflow.
        Higher numbers increase verbosity of messages while zero suppresses all
        messages. Default: 1
    n_proc : int, optional
        Number of processors to use to download AHBA data. Can parallelize up
        to six times. Default: 1

    Returns
    -------
    expression : (R, G) pandas.DataFrame
        Microarray expression for `R` regions in `atlas` for `G` genes,
        aggregated across donors, where the index corresponds to the unique
        integer IDs of `atlas` and the columns are gene names. If
        ``return_donors=True`` then this is a list of (R, G) dataframes, one
        for each donor.
    counts : (R, D) pandas.DataFrame
        Number of samples assigned to each of `R` regions in `atlas` for each
        of `D` donors (if multiple donors were specified); only returned if
        ``return_counts=True``.
    report : str
        Methods describing processing procedures implemented to generate
        `expression`, suitable to be used in a manuscript Methods section. Only
        returned if ``return_report=True``.

    Notes
    -----
    The following methods can be used for collapsing across probes when
    multiple probes are available for the same gene:

        1. ``probe_selection='average'``

        Takes the average of expression data across all probes indexing the
        same gene. Providing 'mean' as the input method will return the same
        thing. This method can only be used when `donor_probes='aggregate'`.

        2. ``probe_selection='max_intensity'``

        Selects the probe with the maximum average expression across samples
        from all donors.

        3. ``probe_selection='max_variance'``

        Selects the probe with the maximum variance in expression across
        samples from all donors.

        4. ``probe_selection='pc_loading'``

        Selects the probe with the maximum loading along the first principal
        component of a decomposition performed across samples from all donors.

        5. ``probe_selection='corr_intensity'``

        Selects the probe with the maximum correlation to other probes from the
        same gene when >2 probes exist; otherwise, uses the same procedure as
        `max_intensity`.

        6. ``probe_selection='corr_variance'``

        Selects the probe with the maximum correlation to other probes from the
        same gene when >2 probes exist; otherwise, uses the same procedure as
        `max_varance`.

        7. ``probe_selection='diff_stability'``

        Selects the probe with the most consistent pattern of regional
        variation across donors (i.e., the highest average correlation across
        brain regions between all pairs of donors). This method can only be
        used when `donor_probes='aggregate'`.

        8. ``method='rnaseq'``

        Selects probes with most consistent pattern of regional variation to
        RNAseq data (across the two donors with RNAseq data). This method can
        only be used when `donor_probes='aggregate'`.

    Note that for incompatible combinations of `probe_selection` and
    `donor_probes` (as detailed above), the `probe_selection choice will take
    precedence. For example, providing ``probe_selection='diff_stability'`` and
    ``donor_probes='independent'`` will cause `donor_probes` to be reset to
    `'aggregate'`.

    The following methods can be used for normalizing microarray expression
    values prior to aggregating:

        1. ``{sample,gene}_norm=='rs'``

        Uses a robust sigmoid function as in [A3]_ to normalize values

        2. ``{sample,gene}_norm='srs'``

        Same as 'rs' but scales output to the unit normal (i.e., range 0-1)

        3. ``{sample,gene}_norm='minmax'``

        Scales data to the unit normal (i.e., range 0-1)

        4. ``{sample,gene}_norm='center'``

        Removes the mean of expression values

        5. ``{sample,gene}_norm='zscore'``

        Applies a basic z-score (subtract mean, divide by standard deviation);
        uses degrees of freedom equal to one for standard deviation

    References
    ----------
    .. [A1] Arnatkevic̆iūtė, A., Fulcher, B. D., & Fornito, A. (2019). A
       practical guide to linking brain-wide gene expression and neuroimaging
       data. NeuroImage, 189, 353-367.
    .. [A2] Hawrylycz, M.J. et al. (2012) An anatomically comprehensive atlas
       of the adult human transcriptome. Nature, 489, 391-399.
    .. [A3] Fulcher, B. D., & Fornito, A. (2016). A transcriptional signature
       of hub connectivity in the mouse connectome. Proceedings of the National
       Academy of Sciences, 113(5), 1435-1440.
    """

    # set logging verbosity level
    LGR.setLevel(dict(zip(range(3), [40, 20, 10])).get(int(verbose), 10))

    # get combination functions
    agg_metric = utils.check_metric(agg_metric)

    # check probe_selection input
    if probe_selection not in probes_.SELECTION_METHODS:
        raise ValueError('Provided probe_selection method is invalid, must be '
                         f'one of {list(probes_.SELECTION_METHODS)}. Received '
                         f'value: {probe_selection}')
    dprobe_opts = ('aggregate', 'independent', 'common')
    if donor_probes not in dprobe_opts:
        raise ValueError('Provided donor_probes method is invalid, must be one'
                         f' of {dprobe_opts}. Received value: {donor_probes}')

    # check lr_mirror input
    if isinstance(lr_mirror, bool):
        warnings.warn('Setting lr_mirror to a boolean value will be '
                      'deprecated in an upcoming release. Use either '
                      '`lr_mirror=None` or `lr_mirror="bidirectional"` '
                      'instead. See documentation for more details.',
                      FutureWarning)
        lr_mirror = 'bidirectional' if lr_mirror else None
    mirror_opts = (None, 'bidirectional', 'leftright', 'rightleft')
    if lr_mirror not in mirror_opts:
        raise ValueError('Provided lr_mirror method is invalid, must be one '
                         f'of {mirror_opts}. Received value: {lr_mirror}')

    # check exact/missing input
    if exact is not None:
        warnings.warn('The exact parameter will be deprecated in an upcoming '
                      'release. Use the `missing` parameter instead. See '
                      'documentation for more details.',
                      FutureWarning)
        missing = None if exact else 'centroids'
    missing_opts = (None, 'interpolate', 'centroids')
    if missing not in missing_opts:
        raise ValueError('Provided missing method is invalid, must be one '
                         f'of {missing_opts}. Received value: {missing}')

    # check that region_agg is viable with return_donors
    if return_donors and region_agg == 'samples':
        raise ValueError('Cannot return donor-level expresison data when '
                         'region_agg parameter is set to "samples".')

    # load atlas and atlas_info, if provided, and coerce to dict
    atlas, group_atlas = images.coerce_atlas_to_dict(atlas, donors, atlas_info)

    # fetch files (downloading if necessary) and unpack to variables
    files = datasets.fetch_microarray(data_dir=data_dir, donors=donors,
                                      verbose=verbose, n_proc=n_proc)

    if probe_selection == 'diff_stability' and len(files) == 1:
        raise ValueError('Cannot use diff_stability for probe_selection with '
                         'only one donor. Please specify a different probe_'
                         'selection method or use more donors.')
    elif probe_selection == 'rnaseq':  # fetch RNAseq if we're gonna need it
        datasets.fetch_rnaseq(data_dir=data_dir, donors=donors,
                              verbose=verbose)

    # get some info on labels in `atlas`
    all_labels = utils.first_entry(atlas).labels
    n_gb = (8 * len(all_labels) * 30000) / (1024 ** 3)
    if n_gb > 1:
        LGR.warning(f'Output matrix may require up to {n_gb:.2f} GB RAM')

    # get dataframe of probe information (reannotated or otherwise)
    # the Probes.csv files are the same for every donor so just grab the first
    probe_info = io.read_probes(first_entry(files, 'probes'))
    if reannotated:
        probe_info = probes_.reannotate_probes(probe_info)

    # drop probes with no/invalid Entrez ID
    probe_info = probe_info.dropna(subset=['entrez_id'])

    # update the annotation "files". this handles updating the MNI coordinates,
    # dropping mistmatched samples (where MNI coordinates don't match the
    # provided ontology), and mirroring samples across hemispheres, if desired
    for donor, data in files.items():
        annot, ontol = data['annotation'], data['ontology']
        t1w = None
        if not group_atlas:
            t1w = datasets.fetch_raw_mri(donors=donor,
                                         data_dir=data_dir,
                                         verbose=verbose)[donor]['t1w']
        annot = samples_.update_coords(annot, corrected_mni=corrected_mni,
                                       native_space=t1w, atlas=atlas[donor])
        if lr_mirror is not None:
            annot = samples_.mirror_samples(annot, ontol, swap=lr_mirror)
        annot = samples_.drop_mismatch_samples(annot, ontol)
        if sim_threshold is not None:
            annot = samples_.similarity_threshold(data['microarray'],
                                                  annot, probe_info,
                                                  threshold=sim_threshold)
        data['annotation'] = annot
    annotation = flatten_dict(files, 'annotation')

    # intensity-based filtering of probes
    probe_info = probes_.filter_probes(flatten_dict(files, 'pacall'),
                                       annotation, probe_info,
                                       threshold=ibf_threshold)

    # get probe-reduced microarray expression data for all donors based on
    # selection method; this will be a list of gene x sample dataframes (one
    # for each donor)
    microarray = probes_.collapse_probes(flatten_dict(files, 'microarray'),
                                         annotation, probe_info,
                                         method=probe_selection,
                                         donor_probes=donor_probes)
    centroids = []
    counts = pd.DataFrame(np.zeros((len(all_labels) + 1, len(microarray)),
                                   dtype=int),
                          index=np.append([0], all_labels),
                          columns=microarray.keys())
    for subj in microarray:
        if lr_mirror is not None:  # reset index
            idx = pd.Series(range(1, len(microarray[subj]) + 1),
                            name='sample_id')
            microarray[subj].index = idx.copy()
            annotation[subj].index = idx.copy()

        # assign samples to regions
        labels = atlas[subj].label_samples(annotation[subj], tolerance)

        # if we're doing exact matching and want to aggregate samples w/i
        # regions, remove the non-labelled samples prior to normalization.
        # otherwise, we'll remove the non-labelled samples after normalization
        nz = np.asarray(labels != 0).squeeze()
        if nz.sum() == 0:
            LGR.warning(f'No samples matched to atlas for donor {subj}')
            microarray[subj].index = labels['label']
            if missing == 'centroids':
                missing += [(pd.DataFrame(), {})]
            continue

        if norm_matched:
            microarray[subj] = microarray[subj].loc[nz]
            annotation[subj] = annotation[subj].loc[nz]
            labels = labels.loc[nz]

        # if normalizing by structural class get annotation dataframe
        annot = annotation[subj][['structure']] if norm_structures else None
        if sample_norm is not None:
            microarray[subj] = correct.normalize_expression(microarray[subj].T,
                                                            norm=sample_norm,
                                                            ignore_warn=True).T
        if gene_norm is not None:
            microarray[subj] = correct.normalize_expression(microarray[subj],
                                                            norm=gene_norm,
                                                            structures=annot,
                                                            ignore_warn=True)

        # get counts of samples collapsed into each ROI
        labs, num = np.unique(labels, return_counts=True)
        counts.loc[labs, subj] = num
        LGR.info(f'{counts.iloc[1:][subj].sum():>3} / {len(nz)} '
                 f'samples matched to regions for donor #{subj}')

        # if we don't want to do exact matching then cache which parcels are
        # missing data and the expression data for the closest sample to that
        # parcel; we'll use this once we've iterated through all donors
        if missing is not None:
            if missing == 'interpolate':
                exp, lab = _interpolate_missing(atlas[subj], labels,
                                                microarray[subj],
                                                annotation[subj])
                microarray[subj] = microarray[subj].append(exp)
                labels = labels.append(lab)
            elif missing == 'centroids':
                centroids += [_match_centroids(atlas[subj], labels,
                                               microarray[subj],
                                               annotation[subj])]

        microarray[subj] = pd.merge(microarray[subj], labels, on='sample_id') \
                             .set_index('label') \
                             .rename_axis('gene_symbol', axis=1)

    # if we don't want to aggregate over regions return sample-level results
    if region_agg is None:
        microarray = {
            donor: micro.set_index(annotation[donor]['well_id'], append=True)
                        .dropna(axis=1, how='any')
            for donor, micro in microarray.items()
        }
        if not return_donors:
            microarray = pd.concat(microarray.values())
        return microarray

    if missing == 'centroids':
        # labels that are missing across all donors
        empty = reduce(set.intersection, [set(f.index) for f, d in centroids])
        LGR.info(f'Matching {len(empty)} region(s) with no data to the '
                 'nearest tissue sample(s)')
        for roi in empty:
            # take weighted average of samples closest to parcel centroid
            # across subjects
            exp, dist = zip(*[(micro.loc[[roi]], dist.get(roi))
                              for micro, dist in centroids])
            weights = _get_weights(np.asarray(dist)[:, None])
            exp = pd.DataFrame(
                (pd.concat(exp) * weights).sum(axis=0) / weights.sum(),
                columns=pd.Series([roi], name='label')
            ).T
            # assign same value to every subject
            for subj in microarray:
                microarray[subj] = microarray[subj].append(exp)

    microarray = samples_.aggregate_samples(microarray,
                                            labels=all_labels,
                                            region_agg=region_agg,
                                            agg_metric=agg_metric,
                                            return_donors=return_donors)

    counts = counts.drop([0], axis=0)
    if return_report:  # generate report
        report = reporting.Report(atlas, atlas_info=atlas[subj].atlas_info,
                                  ibf_threshold=ibf_threshold,
                                  probe_selection=probe_selection,
                                  donor_probes=donor_probes,
                                  lr_mirror=lr_mirror, missing=missing,
                                  tolerance=tolerance, sample_norm=sample_norm,
                                  gene_norm=gene_norm,
                                  norm_matched=norm_matched,
                                  norm_structures=norm_structures,
                                  region_agg=region_agg, agg_metric=agg_metric,
                                  corrected_mni=corrected_mni,
                                  reannotated=reannotated, donors=donors,
                                  return_donors=return_donors,
                                  data_dir=data_dir, counts=counts,
                                  n_probes=len(probe_info),
                                  n_genes=(first_entry(microarray).shape[1]
                                           if return_donors
                                           else microarray.shape[1])).body

    # pack outputs
    out = (microarray,)
    if return_counts:
        out += (counts,)
    if return_report:
        out += (report,)
    if len(out) == 1:
        out = out[0]

    return out


def get_samples_in_mask(mask=None, **kwargs):
    """
    Returns preprocessed microarray expression data for samples in `mask`

    Uses the same processing workflow as :func:`abagen.get_expression_data` but
    instead of aggregating samples within regions simply returns sample-level
    expression data for all samples that fall within boundaries of `mask`.

    Parameters
    ----------
    mask : niimg-like object or dict, optional
        A mask image in MNI space or a tuple of GIFTI images in fsaverage5
        space (where 0 is the background). Alternatively, a dictionary where
        keys are donor IDs and values are mask images (or surfaces) in the
        native space of each donor. If not supplied, all available samples will
        be returned. Default: None
    kwargs : key-value pairs
        All key-value pairs from :func:`abagen.get_expression_data` except for:
        `atlas`, `atlas_info`, `region_agg`, and `agg_metric`, which will be
        ignored. If `atlas` is supplied instead of `mask` then `atlas` will be
        used instead as a modified binary image. If both `atlas` and `mask` are
        supplied then `mask` will be used

    Returns
    -------
    expression : (S, G) pandas.DataFrame
        Microarray expression for `S` samples for `G` genes, aggregated across
        donors, where the columns are gene names
    coords : (S,) numpy.ndarray
        MNI coordinates of samples in `expression`. Even if donor-specific
        masks are provided MNI coordinates will be returned to ensure
        comparability between subjects
    """

    # fetch files (downloading if necessary) to get coordinates
    files = datasets.fetch_microarray(data_dir=kwargs.get('data_dir', None),
                                      donors=kwargs.get('donors', 'all'),
                                      verbose=kwargs.get('verbose', 1),
                                      n_proc=kwargs.get('n_proc', 1))

    # get updated coordinates
    for donor, data in files.items():
        annot = data['annotation']
        if kwargs.get('corrected_mni', True):
            annot = samples_.update_mni_coords(annot)
        files[donor]['annotation'] = annot
    cols = ['well_id', 'mni_x', 'mni_y', 'mni_z']
    coords = np.asarray(pd.concat(flatten_dict(files, 'annotation'))[cols])
    well_id, coords = np.asarray(coords[:, 0], 'int'), coords[:, 1:]
    coords = pd.DataFrame(coords, columns=['x', 'y', 'z'], index=well_id)

    # in case people mix things up and use atlas instead of mask, use that
    if kwargs.get('atlas') is not None and mask is None:
        mask = kwargs['atlas']
    elif mask is None:
        # generate atlas image where each voxel with sample has value of 1
        affine = np.eye(4)
        affine[:-1, -1] = -1 * np.floor(np.max(np.abs(coords), axis=0))
        img = np.ones(np.asarray(-2 * affine[:-1, -1], dtype=int))
        mask = nib.Nifti1Image(img, affine=affine)

    # reset these parameters
    kwargs['atlas'] = mask
    kwargs['atlas_info'] = None
    kwargs['region_agg'] = None
    # soft reset this parameter
    kwargs.setdefault('norm_matched', False)

    # get expression data + drop sample coordinates that weren't in atlas
    exp = get_expression_data(**kwargs)
    if kwargs.get('return_donors'):
        exp = {
            donor: micro.drop(index=[0], level='label', errors='ignore')
                        .droplevel('label')
            for donor, micro in exp.items()
        }
        coords = {
            donor: coords.loc[micro.index]
            for donor, micro in exp.items()
        }
    else:
        exp = exp.drop(index=[0], level='label', errors='ignore') \
                 .droplevel('label')
        coords = coords.loc[exp.index]

    return exp, coords


def get_interpolated_map(genes, mask, n_neighbors=10, **kwargs):
    """
    Generates dense (i.e., interpolated) expression maps for `genes`

    Uses k-nearest neighbors regression to estimate values at every voxel or
    vertex in the supplied `mask`. Note that by default `lr_mirror` is set to
    'bidirectional' and `norm_matched` is set to False, unless explicitly
    specified.

    Parameters
    ----------
    genes : (G,) list-of-str
        List of gene acronyms for which dense maps are desired
    mask : niimg-like object, optional
        A mask image in MNI space or a tuple of GIFTI images in fsaverage5
        space (where 0 is the background)
    n_neighbors : int, optional
        Number of neighboring tissue samples to use when interpolating data in
        dense map. Default: 10
    kwargs : key-value pairs
        All key-value pairs from :func:`abagen.get_expression_data` except for:
        `atlas`, `atlas_info`, `region_agg`, and `agg_metric`, which will be
        ignored.

    Returns
    -------
    dense : (G,) dict
        Dictionary where keys are `genes` and values are dense maps in space of
        provided `mask`
    """

    genes = np.atleast_1d(genes)
    mask = images.check_atlas(mask)
    if 'atlas' in kwargs:
        del kwargs['atlas']

    # soft reset here and then get ALL the samples in the brain
    # (we don't just want samples _inside_ the mask since we're using NN)
    kwargs.setdefault('lr_mirror', 'bidirectional')
    expression, coords = get_samples_in_mask(mask=None, **kwargs)

    # this is a "naive" matching based only on Euclidean distance
    # TODO: use surface distances when working with surfaces
    # TODO: constrain by hemisphere (and structure?)
    exptree = matching.AtlasTree(np.asarray(expression.index),
                                 coords=np.asarray(coords))
    dist, idx = exptree.tree.query(mask.coords, k=n_neighbors)
    dist = _get_weights(dist)

    # get average of nearest neighbors
    dense = np.zeros(((genes.shape[0],) + mask._shape))
    denom = np.sum(dist, axis=1)
    for n, j in enumerate(genes):
        num = np.sum(np.asarray(expression[j])[idx] * dist, axis=1)
        dense[n][mask._nz] = num / denom

    return dict(zip(genes, dense))


def _get_weights(dist):
    """ Gets inverse of `dist`, handling potential infs

    Parameters
    ----------
    dist : array_like
        Distances to be converted to weights

    Returns
    -------
    weights : np.ndarray
        Inverse of `dist`
    """

    with np.errstate(divide='ignore'):
        dist = 1. / dist
    isinf = np.isinf(dist)
    infrow = np.any(isinf, axis=1)
    dist[infrow] = isinf[infrow]

    return dist


def _interpolate_missing(atlas, labels, microarray, annotation):
    """ Interpolate missing data in `atlas` with weighted nearest-neighbors
    """

    data, labs = [], []
    for label in np.setdiff1d(atlas.labels, labels):
        # it's more cost efficient to do this check here
        if atlas.atlas_info is not None:
            cols = ['structure', 'hemisphere']
            match = atlas.atlas_info.loc[label, cols] == annotation[cols]
            if not np.any(np.all(match, axis=1)):
                continue

        # get sample indices for every node of `lab` in `atlas`
        idx, dist = atlas.fill_label(annotation, label, return_dist=True)
        if np.all(np.isinf(dist)):
            continue

        # get expression of selected tissue samples and weights
        dist = _get_weights(dist[:, None])

        # get weighted average of sample expression and store
        exp = (microarray.loc[idx] * dist).sum(axis=0) / dist.sum()
        sampid = pd.Series([f'interp_{label}'], name='sample_id')
        data.append(pd.DataFrame(exp, columns=sampid).T)
        labs.append(pd.DataFrame(label, columns=['label'], index=sampid))

    return pd.concat(data), pd.concat(labs)


def _match_centroids(atlas, labels, microarray, annotation):
    """ Matches missing data in `atlas` based on centroids
    """

    cols = ['mni_x', 'mni_y', 'mni_z']
    tree = matching.AtlasTree(
        np.asarray(annotation.index),
        coords=np.asarray(annotation[cols]),
        atlas_info=annotation.rename_axis('id')
    )
    empty = np.setdiff1d(atlas.labels, labels)

    # get coordinates of centroids for each missing parcel
    centroids = np.r_[[atlas.centroids[lab] for lab in empty]]
    if atlas.atlas_info is not None:
        centroid_info = atlas.atlas_info.loc[empty]
        centroid_info[cols] = centroids
        centroids = centroid_info

    # get idx of tissue sample closest to each centroid
    idx, dist = tree.match_closest_centroids(centroids, return_dist=True)

    # get expression of selected tissue samples
    mask = idx != -1
    exp = microarray.loc[idx[mask]]
    exp.index = pd.Series(empty[mask], name='label')
    # if some labels weren't matched (due to e.g., hemispheric/structural
    # constraints) then append empty rows to the dataframe
    exp = exp.append(pd.DataFrame(index=empty[~mask])) \
             .sort_index() \
             .rename_axis('gene_symbol', axis=1)

    return (exp, dict(zip(empty, dist)))
