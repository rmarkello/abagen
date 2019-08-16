# -*- coding: utf-8 -*-
"""
Functions for cleaning and processing the AHBA microarray dataset
"""

from nibabel.volumeutils import Recoder
import numpy as np
import pandas as pd

from . import io, utils
from .datasets import _fetch_alleninf_coords

# AHBA structure IDs corresponding to different brain parts
ONTOLOGY = Recoder(
    (('4008', 'cerebral cortex', 'cortex'),
     ('4275', 'cerebral nuclei', 'subcortex'),
     ('4391', 'diencephalon', 'subcortex'),
     ('9001', 'mesencephalon', 'subcortex'),
     ('4696', 'cerebellum', 'cerebellum'),
     ('4833', 'metencephalon', 'cerebellum'),
     ('9512', 'myelencephalon', 'cerebellum')),
    fields=('id', 'name', 'structure')
)


def drop_mismatch_samples(annotation, ontology, corrected=True):
    """
    Removes samples from `annotation` whose coordinates do not match `ontology`

    Checks MNI coordinates specified in `annotation` and matches them to L/R
    hemisphere designation in `ontology`; samples who do not match (e.g., L
    hemisphere designation in `ontology` but X coordinate < 0) are removed.

    Optionally updates MNI coordinates in `annotation` (see `corrected`).

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
    corrected : bool, optional
        Whether to use the "corrected" MNI coordinates shipped with the
        `alleninf` package instead of the coordinates provided with the AHBA
        data. Default: True

    Returns
    -------
    keep : pandas.DataFrame
        Annotation data
    """

    # read in data files
    annot = io.read_annotation(annotation)
    ont = io.read_ontology(ontology)

    # add hemisphere + brain "structure" designation to annotation data
    hemi = dict(zip(ont['id'], ont['hemisphere']))
    stru = dict(zip(ont['id'], ont['structure_id_path'].apply(_get_struct)))
    annot = annot.assign(hemisphere=annot['structure_id'].replace(hemi),
                         structure=annot['structure_id'].replace(stru))

    # correct MNI coordinates, if requested
    if corrected:
        annot = _replace_mni_coords(annot)

    # only keep samples with consistent hemisphere + MNI coordinate designation
    keep = annot.query('(hemisphere == "L" & mni_x < 0)'
                       '| (hemisphere == "R" & mni_x > 0)')

    return keep


def _replace_mni_coords(annotation):
    """
    Replaces MNI coords in `annotation` with corrected coords from alleninf

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
    """

    annotation = io.read_annotation(annotation)
    mni_coords = _fetch_alleninf_coords().loc[annotation.well_id]
    annotation[['mni_x', 'mni_y', 'mni_z']] = np.asarray(mni_coords)

    return annotation


def _get_struct(structure_path):
    """
    Gets overarching "structure" of region defined by `structure_path`

    Structure here is defined as one of ['cortex', 'subcortex', 'cerebellum']

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
        One of ['cortex', 'subcortex', 'cerebellum'] or None
    """

    structure_path = set(structure_path.split('/'))
    ids = list(set(ONTOLOGY.value_set('id')).intersection(structure_path))

    try:
        return ONTOLOGY.structure[ids[0]]
    except IndexError:
        return


def normalize_expression(expression):
    """
    Performs scaled robust sigmoid (SRS) normalization on `expression` data

    Parameters
    ----------
    expression : (S, G) pandas.DataFrame
        Microarray expression data, where `S` is samples (or regions) and `G`
        is genes

    Returns
    -------
    normalized : (S, G) pandas.DataFrame
        Data from `expression` normalized separately for each gene
    """

    # get non-NaN values
    data = np.asarray(expression.dropna(axis=0, how='all'))

    # calculate sigmoid normalization
    norm = (data - np.median(data, axis=0)) / np.std(data, axis=0)
    srs = 1 / (1 + np.exp(-norm))

    # rescale normalized values to a unit interval
    srs_min = srs.min(axis=0, keepdims=True)
    srs_max = srs.max(axis=0, keepdims=True)
    scaled = (srs - srs_min) / (srs_max - srs_min)

    # recreate dataframe and fill non-NaN values
    normalized = pd.DataFrame(np.nan,
                              columns=expression.columns,
                              index=expression.index)
    inds = expression[expression.notna().all(axis=1)].index
    normalized.loc[inds] = scaled

    return normalized


def aggregate_donors(expression, metric='mean'):
    """
    Aggregates microarray `expression` across donors using `metric`

    Parameters
    ----------
    expression : list of (R, G) pandas.DataFrame
        Where each entry is the microarray expression of `R` regions across `G`
        genes for a given donor
    metric : {'mean', 'median'} or callable, optional
        Mechanism by which to reduce donor-level expression data into a single
        dataframe. If a callable, should be able to accept an `N`-dimensional
        input and the `axis` keyword argument and return an `N-1`-dimensional
        output. Default: 'mean'

    Returns
    -------
    expression : (R, G) pandas.DataFrame
        Microarray expression for `R` regions in `atlas` for `G` genes,
        aggregated across donors.
    """

    metric = utils.check_metric(metric)

    return pd.concat(expression).groupby('label').aggregate(np.mean)
