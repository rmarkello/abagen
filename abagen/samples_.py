# -*- coding: utf-8 -*-
"""
Functions for cleaning and processing the AHBA microarray dataset
"""

import logging
import os
from pkg_resources import resource_filename

import nibabel as nib
import numpy as np
import pandas as pd

from . import io, transforms, utils

LGR = logging.getLogger('abagen')

# AHBA structure IDs corresponding to different brain parts
ONTOLOGY = nib.volumeutils.Recoder(
    (('4008', 'cerebral cortex', 'cortex'),
     ('4249', 'hippocampal formation', 'subcortex/brainstem'),
     ('4275', 'cerebral nuclei', 'subcortex/brainstem'),
     ('4391', 'diencephalon', 'subcortex/brainstem'),
     ('4696', 'cerebellum', 'cerebellum'),
     ('9001', 'mesencephalon', 'subcortex/brainstem'),
     ('9131', 'pons', 'subcortex/brainstem'),
     ('9218', 'white matter', 'white matter'),
     ('9352', 'sulci & spaces', 'other'),
     ('9512', 'myelencephalon', 'subcortex/brainstem')),
    fields=('id', 'name', 'structure')
)


def update_mni_coords(annotation):
    """
    Replaces MNI coords in `annotation` with corrected coords from `alleninf`

    Parameters
    ----------
    annotation : str
        Annotation file from Allen Brain Institute. Optimally obtained by
        calling `abagen.fetch_microarray()` and accessing the `annotation`
        attribute on the resulting object

    Returns
    -------
    corrected : pandas.DataFrame
        Annotation data with corrected MNI coordinates

    References
    ----------
    Updated MNI coordinates taken from https://github.com/chrisfilo/alleninf,
    which is licensed under the BSD-3 (reproduced here):

    Copyright (c) 2018, Krzysztof Gorgolewski
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.
    * Neither the name of the copyright holder nor the names of its
      contributors may be used to endorse or promote products derived from
      this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
    LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.
    """

    coords = resource_filename(
        'abagen', os.path.join('data', 'corrected_mni_coordinates.csv.gz')
    )
    coords = pd.read_csv(coords).rename(dict(corrected_mni_x='mni_x',
                                             corrected_mni_y='mni_y',
                                             corrected_mni_z='mni_z'),
                                        axis=1)
    coords = coords.set_index('well_id')

    annotation = io.read_annotation(annotation, copy=True)

    # basic check that all well_ids in annotation are present in coords
    # a future pandas update may cause this to raise a KeyError but we want
    # this to raise a KeyError NOW
    diff = np.setdiff1d(annotation['well_id'], coords.index)
    if len(diff) > 0:
        raise KeyError('Provided annotation file has well IDs that do not '
                       'exist in updated MNI coordinate file from `alleninf`. '
                       'Please check input annotation file and try again. '
                       'Unknown well IDs: {}'.format(diff))

    mni_coords = coords.loc[annotation.well_id]
    annotation[['mni_x', 'mni_y', 'mni_z']] = np.asarray(mni_coords)

    return annotation


def update_coords(annotation, corrected_mni=True, native_space=None,
                  atlas=None):
    """
    Updates coordinates in `annotation`

    Parameters
    ----------
    annotation : str
        Annotation file from Allen Brain Institute. Optimally obtained by
        calling `abagen.fetch_microarray()` and accessing the `annotation`
        value for one of the returned donors
    corrected_mni : bool
        Whether to replace MNI coordinates in `annotation` with coordinates
        from `alleninf` (uses :func:`update_mni_coords`)
    native_space : niimg-like
        Whether to replace MNI coordinates in `annotation` with native-space
        coordinates; uses `affine` of provided image to convert coordinates

    Returns
    -------
    corrected : pandas.DataFrame
        Annotation data with updated coordinates
    """

    from .images import check_atlas

    annotation = io.read_annotation(annotation, copy=True)

    if corrected_mni:
        annotation = update_mni_coords(annotation)

    if atlas is not None:
        atlas = check_atlas(atlas)
        if atlas.volumetric:
            cols = ['mni_x', 'mni_y', 'mni_z']
            vox_size = 1 / atlas._volumetric
            annotation[cols] = np.floor(annotation[cols] * vox_size) / vox_size

    if native_space is not None:
        try:
            native_space = nib.load(native_space)
        except TypeError:
            native_space = native_space
        affine = native_space.affine
        ijk = annotation[['mri_voxel_x', 'mri_voxel_y', 'mri_voxel_z']]
        annotation[['mni_x', 'mni_y', 'mni_z']] = \
            transforms.ijk_to_xyz(ijk, affine)

    return annotation


def _get_struct(structure_path):
    """
    Gets overarching "structure" of region defined by `structure_path`

    Structure here is defined as being one of {'cortex', 'subcortex',
    'cerebellum', 'brainstem', 'white matter', 'other'}.

    Parameters
    ----------
    structure_path : str
        As obtained from an ontology file from Allen Brain Institute. Optimally
        obtained by calling `abagen.fetch_microarray()`, loading the
        `ontology` attribute on the resulting object, and querying the
        `structure_id_path` column from the resulting dataframe

    Returns
    -------
    structure : str
        Structure, or None if unable to identify a corresponding structure
    """

    structures = set(structure_path.strip('/').split('/'))
    ids = sorted(set(ONTOLOGY.value_set('id')) & structures,
                 key=lambda x: structure_path.index(x))

    try:
        return ONTOLOGY.structure[ids[-1]]
    except IndexError:
        return


def drop_mismatch_samples(annotation, ontology):
    """
    Removes samples from `annotation` whose coordinates do not match `ontology`

    Checks MNI coordinates specified in `annotation` and matches them to L/R
    hemisphere designation in `ontology`; samples who do not match (e.g., L
    hemisphere designation in `ontology` but X coordinate < 0) are removed.

    Parameters
    ----------
    annotation : str
        Annotation file from Allen Brain Institute. Optimally obtained by
        calling `abagen.fetch_microarray()` and accessing the `annotation`
        attribute on the resulting object
    ontology : str
        Ontology file from Allen Brain Institute. Optimally obtained by
        calling `abagen.fetch_microarray()` and accessing the `ontology`
        attribute on the resulting object

    Returns
    -------
    annot : pandas.DataFrame
        Annotation data
    """

    # read in data files
    annotation = io.read_annotation(annotation)
    ontology = io.read_ontology(ontology).set_index('id')
    sid = np.asarray(annotation['structure_id'])

    # get hemisphere and structure path
    hemisphere = np.asarray(ontology.loc[sid, 'hemisphere']
                                    .replace({np.nan: 'B'}))
    structure = np.asarray(ontology.loc[sid, 'structure_id_path']
                                   .apply(_get_struct))

    # add hemisphere + brain "structure" designation to annotation data and
    # only keep samples with consistent hemisphere + MNI coordinate designation
    annot = annotation.assign(hemisphere=hemisphere, structure=structure) \
                      .query('(hemisphere == "L" & mni_x < 0) '
                             '| (hemisphere == "R" & mni_x > 0) '
                             '| (hemisphere == "B" & mni_x == 0)',
                             engine='python')

    return annot


def _mirror_ontology(annotation, ontology):
    """
    Assumes all hemisphere assignments of structures in `annotation` are wrong

    Parameters
    ----------
    annotation,ontology : str
        Filepath to {annotation,ontology} file

    Returns
    -------
    annotation : pandas.DataFrame
        Loaded annotation with updated structure information
    """

    HEMI_SWAP = dict(L='R', R='L')

    annotation = io.read_annotation(annotation, copy=True)
    ontology = io.read_ontology(ontology)

    # structure IDs are specific to structure + hemisphere, so we can use this
    # to grab the hemisphere designation for each sample and flip that L<->R
    hemi = ontology.set_index('id') \
                   .loc[annotation['structure_id'], 'hemisphere'] \
                   .replace(HEMI_SWAP)

    # structure acronyms are distinct to structure but not to hemisphere, so we
    # can use this to grab all variations of a given structure and, with the
    # flipped hemisphere designations we just generated, find the ID/name of
    # the relevant structure in the new hemisphere
    acr = ontology.set_index(['acronym', 'hemisphere']) \
                  .loc[zip(annotation['structure_acronym'], hemi)] \
                  .reset_index()

    # update the original annotation with the new values
    annotation['structure_acronym'] = np.asarray(acr['acronym'])
    annotation['structure_id'] = np.asarray(acr['id'])
    annotation['structure_name'] = np.asarray(acr['name'])

    return annotation


def mirror_samples(annotation, ontology, swap='bidirectional'):
    """
    Mirrors all tissue samples across hemispheres

    Parameters
    ----------
    annotation : str or pandas.DataFrame
        Filepath to annotation file from Allen Brain Institute (i.e., as
        obtained by calling :func:`abagen.fetch_microarray` and accessing the
        `annotation` attribute on the resulting object).
    ontology : str or pandas.DataFrame
        Filepath to ontology file from Allen Brain Institute (i.e., as
        obtained by calling :func:`abagen.fetch_microarray` and accessing the
        `ontology` attribute on the resulting object).
    swap : {'bidirectional', 'leftright', 'rightleft'}, optional
        Which hemispheres to mirror, where 'bidirectional' will mirror both
        hemispheres, 'leftright' will mirror the left to the right, and
        'rightleft' will mirror the right to the left. Default: 'bidirectional'

    Returns
    -------
    annotation : pandas.DataFrame
        Loaded input `annotation` data with all samples duplicated across
        hemispheres (where structure IDs and MNI coordinates are updated,
        accordingly)

    References
    ----------
    .. [SA1] Gryglewski, G., Seiger, R., James, G. M., Godbersen, G. M.,
       Komorowski, A., Unterholzner, J., ... & Kasper, S. (2018). Spatial
       analysis and high resolution mapping of the human whole-brain
       transcriptome for integrative analysis in neuroimaging. NeuroImage, 176,
       259-267.

    .. [SA2] Romero-Garcia, R., Whitaker, K. J., Váša, F., Seidlitz, J., Shinn,
       M., Fonagy, P., ... & Vértes, P. E. (2018). Structural covariance
       networks are coupled to expression of genes enriched in supragranular
       layers of the human cortex. NeuroImage, 171, 256-267.
    """

    avail_swap = ('bidirectional', 'leftright', 'rightleft')
    if swap not in avail_swap:
        raise ValueError('Provided hemisphere `swap` value not recognized. '
                         f'Must be one of {avail_swap}. Received: {swap}')

    annotation = io.read_annotation(annotation)
    ontology = io.read_ontology(ontology)

    # take all lh and rh samples and flip x-coordinate
    # also update ontology information (structure_id/acronym/name) as this is
    # used when dropping mismatched samples later in the workflow
    lh = rh = pd.DataFrame()
    if swap in ('leftright', 'bidirectional'):
        lh = _mirror_ontology(annotation[annotation['mni_x'] < 0], ontology)
        lh['mni_x'] *= -1
    if swap in ('rightleft', 'bidirectional'):
        rh = _mirror_ontology(annotation[annotation['mni_x'] > 0], ontology)
        rh['mni_x'] *= -1

    # grow microarray and pacall based on duplicated samples
    annotation = pd.concat([annotation, lh, rh])

    return annotation


def similarity_threshold(microarray, annotation, probes, threshold=5):
    """
    Performs inter-areal similarity filtering of tissue samples

    Parameters
    ----------
    microarray : str or pandas.DataFrame
        Dataframe with `P` rows representing probes and `S` columns
        representing distinct samples, with values indicating microarray
        expression levels
    annotation : str or pandas.DataFrame
        DataFrame with information about `S` tissue samples in `microarray`
    probes : str or pandas.DataFrame
        Dataframe with information about `P` microarray probes in `microarray`
    threshold : (0, np.inf) float, optional
        Threshold for filtering samples. Specifies the standard deviation
        cutoff used to remove samples. Default: 5

    Returns
    -------
    samples : pandas.DataFrame
        Dataframe containing information on samples that should be retained
        according to inter-areal similarity filtering
    """

    annotation = io.read_annotation(annotation)
    probes = io.read_probes(probes)
    microarray = io.read_microarray(microarray, copy=False)

    corrs = np.corrcoef(microarray.loc[probes.index, annotation.index].T)
    sim = np.sum(corrs, axis=1)
    thresh = np.mean(sim) - (threshold * np.std(sim, ddof=1))

    return annotation.loc[sim >= thresh]


def groupby_index(microarray, labels=None, metric='mean'):
    """
    Averages expression data in `microarray` over samples with same label

    Parameters
    ----------
    microarray : (S, G) pandas.DataFrame
        Microarray expression data, where `S` is samples and `G` is genes
    labels : (L,) array_like, optional
        All possible labels for parcellation (to account for possibility that
        some parcels have NO expression data). Default: None
    metric : str or func, optional
        Mechanism by which to collapse across samples within a parcel. If a
        str, should be in ['mean', 'median']; if a function, should be able to
        accept an `N`-dimensional input and the `axis` keyword argument and
        return an `N-1`-dimensional output. Default: 'mean'

    Returns
    -------
    gene_by_label : (L, G) pandas.DataFrame
        Microarray expression data
    """

    # get combination function
    metric = utils.check_metric(metric)

    # get missing labels
    if labels is not None:
        missing = np.setdiff1d(labels, np.unique(microarray.index))
        labels = pd.DataFrame(columns=microarray.columns,
                              index=pd.Series(missing, name='label'))

    gene_by_label = (microarray.groupby('label')
                               .aggregate(metric)
                               .append(labels)
                               .sort_index()
                               .rename_axis('label'))

    # remove "zero" label (if it exists)
    if 0 in gene_by_label.index:
        gene_by_label = gene_by_label.drop([0], axis=0)

    return gene_by_label


def aggregate_samples(microarray, labels=None, region_agg='donors',
                      agg_metric='mean', return_donors=False):
    """
    Aggregates samples in `microarray` belonging to same regions

    Parameters
    ----------
    microarray : list of (S, G) pandas.DataFrame
        Microarray expression data for multiple donor, where `S` is samples and
        `G` is genes. Index of dataframes should identify to which region each
        sample was assigned
    labels : array_like
        Array containing IDs of all possible region labels. If no samples from
        any donor were assigned to a region NaNs are inserted in place of
        expression values
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
    return_donors : bool, optional
        Whether to return donor-level expression arrays instead of aggregating
        expression across donors with provided `agg_metric`. Default: False

    Returns
    -------
    microarray : (R, G) pandas.DataFrame or list of (R, G) pandas.DataFrame
        Microarray expression data aggregated across samples within each of `R`
        regions. If `return_donors=True` then a list of dataframes is returned
    """

    region_aggs = ['donors', 'samples']
    if region_agg not in region_aggs:
        raise ValueError('Provided `region_agg` {} invalid. Must be one of {}'
                         .format(region_agg, region_aggs))
    if region_agg == 'samples' and return_donors:
        raise ValueError('Cannot use region_agg=\'samples\' when '
                         '`return_donors=True`.')
    metric = utils.check_metric(agg_metric)

    LGR.info('Aggregating samples to regions with provided region_agg: {}'
             .format(region_agg))

    if region_agg == 'donors':
        microarray = [
            groupby_index(e, labels=labels, metric=metric) for e in microarray
        ]
    elif region_agg == 'samples':
        microarray = [
            groupby_index(pd.concat(microarray), labels=labels, metric=metric)
        ]

    if not return_donors:
        microarray = pd.concat(microarray).groupby('label').aggregate(metric)
        # some genes may have been poorly normalized; remove these
        drop = microarray.dropna(axis=0, how='all').isna().any(axis=0)
        LGR.info('Dropping {} gene from concatenated expression data due to '
                 'poor normalization'.format(drop.sum()))
        microarray = microarray.drop(drop[drop].index, axis=1)

    return microarray
