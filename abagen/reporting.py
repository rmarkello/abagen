# -*- coding: utf-8 -*-
"""
Functions for generating workflow methods reports

Note: all text contained within this module is released under a `CC-0 license
<https://creativecommons.org/publicdomain/zero/1.0/>`_.
"""

import logging

import numpy as np

from . import __version__
from .datasets import check_donors, fetch_donor_info
from .images import coerce_atlas_to_dict
from .utils import first_entry

LGR = logging.getLogger('abagen')
METRICS = {
    np.mean: 'mean',
    np.median: 'median'
}
REFERENCES = dict(
    A2019N=('Arnatkevic̆iūtė, A., Fulcher, B. D., & Fornito, A. (2019). '
            'A practical guide to linking brain-wide gene expression and '
            'neuroimaging data. Neuroimage, 189, 353-367.'),
    F2016P=('Fulcher, B. D., & Fornito, A. (2016). A transcriptional '
            'signature of hub connectivity in the mouse connectome. '
            'Proceedings of the National Academy of Sciences, 113(5), '
            '1435-1440.'),
    F2013J=('Fulcher, B. D., Little, M. A., & Jones, N. S. (2013). '
            'Highly comparative time-series analysis: the empirical '
            'structure of time series and their methods. Journal of '
            'the Royal Society Interface, 10(83), 20130048.'),
    H2012N=('Hawrylycz, M. J., Lein, E. S., Guillozet-Bongaarts, A. L., '
            'Shen, E. H., Ng, L., Miller, J. A., ... & Jones, A. R. '
            '(2012). An anatomically comprehensive atlas of the adult '
            'human brain transcriptome. Nature, 489(7416), 391-399.'),
    H2015N=('Hawrylycz, M., Miller, J. A., Menon, V., Feng, D., '
            'Dolbeare, T., Guillozet-Bongaarts, A. L., ... & Lein, E. '
            '(2015). Canonical genetic signatures of the adult human '
            'brain. Nature Neuroscience, 18(12), 1832.'),
    M2021B=('Markello, R.D., Arnatkevic̆iūtė, A., Poline, J-B., Fulcher, B. '
            'D., Fornito, A., & Misic, B. (2021). Standardizing workflows in '
            'imaging transcriptomics with the abagen toolbox. Biorxiv.'),
    P2017G=('Parkes, L., Fulcher, B. D., Yücel, M., & Fornito, A. '
            '(2017). Transcriptional signatures of connectomic '
            'subregions of the human striatum. Genes, Brain and '
            'Behavior, 16(7), 647-663.'),
    Q2002N=('Quackenbush, J. (2002). Microarray data normalization and '
            'transformation. Nature Genetics, 32(4), 496-501.'),
    R2018N=('Romero-Garcia, R., Whitaker, K. J., Váša, F., Seidlitz, '
            'J., Shinn, M., Fonagy, P., ... & NSPN Consortium. (2018). '
            'Structural covariance networks are coupled to expression '
            'of genes enriched in supragranular layers of the human '
            'cortex. NeuroImage, 171, 256-267.')
)


class Report:
    """ Generates report of methods for :func:`abagen.get_expression_data()`

    Refer to the doc-string of the workflow for an overview of paramter options
    """

    def __init__(self, atlas, atlas_info=None, *, ibf_threshold=0.5,
                 probe_selection='diff_stability', donor_probes='aggregate',
                 lr_mirror=None, missing=None, tolerance=2, sample_norm='srs',
                 gene_norm='srs', norm_matched=True, norm_structures=False,
                 region_agg='donors', agg_metric='mean', corrected_mni=True,
                 reannotated=True, donors='all', return_donors=False,
                 data_dir=None, counts=None, n_probes=None, n_genes=None):

        efflevel = LGR.getEffectiveLevel()
        LGR.setLevel(100)

        atlas, self.group_atlas = \
            coerce_atlas_to_dict(atlas, donors, atlas_info=atlas_info,
                                 data_dir=data_dir)
        self.atlas = first_entry(atlas)
        self.atlas_info = self.atlas.atlas_info
        self.ibf_threshold = ibf_threshold
        self.probe_selection = probe_selection
        self.donor_probes = donor_probes
        self.lr_mirror = lr_mirror
        self.missing = missing
        self.tolerance = tolerance
        self.sample_norm = sample_norm
        self.gene_norm = gene_norm
        self.norm_matched = norm_matched
        self.norm_structures = norm_structures
        self.region_agg = region_agg
        self.agg_metric = METRICS.get(agg_metric, agg_metric)
        self.corrected_mni = corrected_mni
        self.reannotated = reannotated
        self.donors = donors
        self.return_donors = return_donors
        self.counts = counts
        self.n_probes = n_probes
        self.n_genes = n_genes
        self.body = self.gen_report()

        LGR.setLevel(efflevel)

    def gen_report(self):
        """ Generates main text of report
        """

        report = ''
        report += """
        Regional microarry expression data were obtained from {n_donors}
        post-mortem brains ({n_female} female, ages {min}--{max}, {mean:.2f}
        +/- {std:.2f}) provided by the Allen Human Brain Atlas (AHBA,
        https://human.brain-map.org; [H2012N]). Data were processed with the
        abagen toolbox (version {vers}; https://github.com/rmarkello/abagen;
        [M2021B])
        """.format(**_get_donor_demographics(self.donors), vers=__version__)

        if self.atlas.volumetric and self.group_atlas:
            report += """
            using a {n_region}-region volumetric atlas in MNI space.<br>
            """.format(n_region=len(self.atlas.labels))
        elif not self.atlas.volumetric and self.group_atlas:
            report += """
            using a {n_region}-region surface-based atlas in MNI space.<br>
            """.format(n_region=len(self.atlas.labels))
        elif self.atlas.volumetric and not self.group_atlas:
            report += """
            using a {n_region}-region volumetric atlas, independently aligned
            to each donor's native MRI space.<br>
            """.format(n_region=len(self.atlas.labels))
        else:
            report += ".<br>"

        if self.reannotated:
            report += """
            First, microarray probes were reannotated using data provided by
            [A2019N]; probes not matched to a valid Entrez ID were discarded.
            """
        else:
            report += """
            First, microarray probes not matched to a valid Entrez ID were
            discarded.
            """

        if self.ibf_threshold > 0:
            report += """
            Next, probes were filtered based on their expression intensity
            relative to background noise [Q2002N], such that probes with
            intensity less than the background in >={threshold:.2f}% of samples
            across donors were discarded
            """.format(threshold=self.ibf_threshold * 100)
            if self.n_probes is not None:
                report += """
                , yielding {n_probes:,} probes
                """.format(n_probes=self.n_probes)
            report += "."

        if self.probe_selection == 'average':
            report += """
            When multiple probes indexed the expression of the same gene we
            calculated the mean expression across probes.
            """
        elif self.probe_selection == 'diff_stability':
            report += r"""
            When multiple probes indexed the expression of the same gene, we
            selected and used the probe with the most consistent pattern of
            regional variation across donors (i.e., differential stability;
            [H2015N]), calculated with:<br>

            $$ \Delta_{{S}}(p) = \frac{{1}}{{\binom{{N}}{{2}}}} \,
            \sum_{{i=1}}^{{N-1}} \sum_{{j=i+1}}^{{N}}
            \rho[B_{{i}}(p), B_{{j}}(p)] $$<br>

            where $ \rho $ is Spearman's rank correlation of the expression of
            a single probe, p, across regions in two donors $B_{{i}}$ and
            $B_{{j}}$, and N is the total number of donors. Here, regions
            correspond to the structural designations provided in the ontology
            from the AHBA.
            """
        elif self.probe_selection == "pc_loading":
            report += """
            When multiple probes indexed the expression of the same gene we
            selected and used the probe with the highest loading on the first
            principal component taken from a decomposition of probe expression
            across samples from all donors [P2017G].
            """
        elif self.probe_selection == "max_intensity":
            report += """
            When multiple probes indexed the expression of the same gene we
            selected and used the probe with the highest mean intensity across
            samples.
            """
        elif self.probe_selection == "max_variance":
            report += """
            When multiple probes indexed the expression of the same gene we
            selected and used the probe with the highest variance across
            samples.
            """
        elif self.probe_selection == "corr_intensity":
            report += """
            When multiple probes indexed the expression of the same gene we
            selected and used the expression profile of a single representative
            probe. For genes where only two probes were available we selected
            the probe with the highest mean intensity across samples; where
            three or more probes were available, we calculated the correlation
            of probe expression across samples and selected the probe with the
            highest average correlation.
            """
        elif self.probe_selection == "corr_variance":
            report += """
            When multiple probes indexed the expression of the same gene we
            selected and used the expression profile of a single representative
            probe. For genes where only two probes were available we selected
            the probe with the highest variance across samples; where three or
            more probes were available, we calculated the correlation of probe
            expression across samples and selected the probe with the highest
            average correlation.
            """
        elif self.probe_selection == "rnaseq":
            report += """
            When multiple probes indexed the expression of the same gene we
            selected and used the probe with the most consistent pattern of
            regional expression to RNA-seq data (available for two donors in
            the AHBA). That is, we calculated the Spearman’s rank correlation
            between each probes' microarray expression and RNA-seq expression
            data of the corresponding gene, and selected the probe with the
            highest correspondence. Here, regions correspond to the structural
            designations provided in the ontology from the AHBA.
            """

        if (self.donor_probes == "aggregate" and self.probe_selection not in
                ['average', 'diff_stability', 'rnaseq']):
            report += """
            The selection of probes was performed using sample expression data
            aggregated across all donors.<br>
            """
        elif (self.donor_probes == "aggregate" and self.probe_selection in
                ['average', 'diff_stability', 'rnaseq']):
            report += "<br>"
        elif self.donor_probes == "independent":
            report += """
            The selection of probes was performed independently for each donor,
            such that the probes chosen to represent each gene could differ
            across donors.<br>
            """
        elif self.donor_probes == "common":
            report += """
            The selection of probes was performed independently for each donor,
            and the probe most commonly selected across all donors was chosen
            to represent the expression of each gene for all donors.<br>
            """

        if self.corrected_mni and self.group_atlas:
            report += """
            The MNI coordinates of tissue samples were updated to those
            generated via non-linear registration using the Advanced
            Normalization Tools (ANTs; https://github.com/chrisfilo/alleninf).
            """

        if self.lr_mirror == 'bidirectional':
            report += """
            To increase spatial coverage, tissue samples were mirrored
            bilaterally across the left and right hemispheres [R2018N].
            """
        elif self.lr_mirror == 'leftright':
            report += """
            To increase spatial coverage, tissue samples in the left hemisphere
            were mirrored into the right hemisphere [R2018N].
            """
        elif self.lr_mirror == 'rightleft':
            report += """
            To increase spatial coverage, tissue samples in the right
            hemisphere were mirrored into the left hemisphere [R2018N].
            """

        if self.tolerance == 0 and not self.atlas.volumetric:
            report += """
            Samples were assigned to brain regions by minimizing the Euclidean
            distance between the {space} coordinates of each sample and the
            nearest surface vertex. Samples where the Euclidean distance to the
            nearest vertex was greater than the mean distance for all samples
            belonging to that donor were excluded.
            """.format(space='MNI' if self.group_atlas else 'native voxel')
        elif self.tolerance == 0 and self.atlas.volumetric:
            report += """
            Samples were assigned to brain regions in the provided atlas only
            if their {space} coordinates were directly within a voxel belonging
            to a parcel.
            """.format(space='MNI' if self.group_atlas else 'native voxel')
        elif self.tolerance > 0 and not self.atlas.volumetric:
            report += """
            Samples were assigned to brain regions by minimizing the Euclidean
            distance between the {space} coordinates of each sample and the
            nearest surface vertex. Samples where the Euclidean distance to the
            nearest vertex was more than {tolerance} standard deviations above
            the mean distance for all samples belonging to that donor were
            excluded.
            """.format(space='MNI' if self.group_atlas else 'native voxel',
                       tolerance=self.tolerance)
        elif self.tolerance > 0 and self.atlas.volumetric:
            report += """
            Samples were assigned to brain regions in the provided atlas if
            their {space} coordinates were within {tolerance} mm of a given
            parcel.
            """.format(space='MNI' if self.group_atlas else 'native voxel',
                       tolerance=self.tolerance)
        elif self.tolerance < 0 and not self.atlas.volumetric:
            report += """
            Samples were assigned to brain regions by minimizing the Euclidean
            distance between the {space} coordinates of each sample and the
            nearest surface vertex. Samples where the Euclidean distance to the
            nearest vertex was more than {tolerance}mm were excluded.
            """.format(space='MNI' if self.group_atlas else 'native voxel',
                       tolerance=abs(self.tolerance))

        if self.atlas_info is not None:
            report += """
            To reduce the potential for misassignment, sample-to-region
            matching was constrained by hemisphere and gross structural
            divisions (i.e., cortex, subcortex/brainstem, and cerebellum, such
            that e.g., a sample in the left cortex could only be assigned to an
            atlas parcel in the left cortex; [A2019N]).
            """

        if self.missing == 'centroids':
            report += """
            If a brain region was not assigned a sample from any donor based on
            the above procedure, the tissue sample closest to the centroid of
            that parcel was identified independently for each donor. The
            average of these samples was taken across donors, weighted by the
            distance between the parcel centroid and the sample, to obtain an
            estimate of the parcellated expression values for the missing
            region.
            """
            if self.counts is not None:
                n_missing = np.sum(self.counts.sum(axis=1) == 0)
                if n_missing > 0:
                    report += """
                    This procedure was performed for {n_missing} regions that
                    were not assigned tissue samples.
                    """.format(n_missing=n_missing)
        elif self.missing == 'interpolate':
            report += """
            If a brain region was not assigned a tissue sample based on the
            above procedure, every {node} in the region was mapped to the
            nearest tissue sample from the donor in order to generate a dense,
            interpolated expression map. The average of these expression values
            was taken across all {nodes} in the region, weighted by the
            distance between each {node} and the sample mapped to it, in order
            to obtain an estimate of the parcellated expression values for the
            missing region.
            """.format(node='voxel' if self.atlas.volumetric else 'vertex',
                       nodes='voxels' if self.atlas.volumetric else 'vertices')

        if self.norm_matched:
            report += """
            All tissue samples not assigned to a brain region in the provided
            atlas were discarded.
            """
        report += "<br>"

        if self.sample_norm is not None:
            report += """
            Inter-subject variation was addressed {}
            """.format(_get_norm_procedure(self.sample_norm, 'sample_norm'))

        if self.gene_norm is None:
            pass
        elif self.gene_norm == self.sample_norm:
            report += """
            Gene expression values were then normalized across tissue samples
            using an identical procedure.
            """
        elif (self.gene_norm != self.sample_norm
                and self.sample_norm is not None):
            report += """
            Inter-subject variation in gene expression was then addressed {}
            """.format(_get_norm_procedure(self.gene_norm, 'gene_norm'))
        elif self.gene_norm != self.sample_norm and self.sample_norm is None:
            report += """
            Inter-subject variation was addressed {}
            """.format(_get_norm_procedure(self.gene_norm, 'gene_norm'))

        if not self.norm_matched and not self.norm_structures:
            report += """
            All available tissue samples were used in the normalization process
            regardless of whether they were assigned to a brain region.
            """
        elif not self.norm_matched and self.norm_structures:
            report += """
            All available tissue samples were used in the normalization process
            regardless of whether they were matched to a brain region; however,
            normalization was performed separately for samples in distinct
            structural classes (i.e., cortex, subcortex/brainstem, cerebellum).
            """
        elif self.norm_matched and self.norm_structures:
            report += """
            Normalization was performed separately for samples in distinct
            structural classes (i.e., cortex, subcortex/brainstem, cerebellum).
            """

        if not self.norm_matched and self.gene_norm is not None:
            report += """
            Tissue samples not matched to a brain region were discarded after
            normalization.
            """

        if self.region_agg == 'donors' and self.agg_metric == 'mean':
            report += """
            Samples assigned to the same brain region were averaged separately
            for each donor{donors}, yielding a regional expression matrix
            {n_donors}
            """.format(donors=' and then across donors' if
                              not self.return_donors else '',
                       n_donors=' for each donor ' if
                                self.return_donors else '')
        elif self.region_agg == 'samples' and self.agg_metric == 'mean':
            report += """
            Samples assigned to the same brain region were averaged across all
            donors, yielding a regional expression matrix
            """
        elif self.region_agg == 'donors' and self.agg_metric == 'median':
            report += """
            The median value of samples assigned to the same brain region was
            computed separately for each donor{donors}, yielding a regional
            expression matrix{n_donors}
            """.format(donors=' and then across donors' if
                              not self.return_donors else '',
                       n_donors=' for each donor' if
                                self.return_donors else '')
        elif self.region_agg == 'samples' and self.agg_metric == 'median':
            report += """
            The median value of samples assigned to the same brain region was
            computed across donors, yielding a single regional expression
            matrix
            """
        if self.n_genes is not None:
            report += """
            with {n_region} rows, corresponding to brain regions, and
            {n_genes:,} columns, corresponding to the retained genes
            """.format(n_region=len(self.atlas.labels),
                       n_genes=self.n_genes)
        report += "\b."

        report += str(_add_references(report))
        return _sanitize_text(report)


def _sanitize_text(text):
    """ Cleans up line and paragraph breaks in `text`
    """

    text = ' '.join(text.replace('\n', ' ').split())
    return text.replace('<br> ', '\n\n') \
               .replace('<p>', '\n')


def _get_donor_demographics(donors):
    """
    Gets demographic info about the requested `donors`

    Parameters
    ----------
    donors : list, optional
        List of donors. Can be either donor numbers or UID. If not specified
        will use all available donors. Note that donors '9861' and '10021'
        have samples from both left + right hemispheres; all other donors have
        samples from the left hemisphere only. Default: 'all'

    Returns
    -------
    info : dict
        With keys ['n_donor', 'n_female','min', 'max', 'mean', 'std']. The
        values represent features of the age distribution of requested donors.
    """

    donors = [int(i) for i in check_donors(donors)]
    info = fetch_donor_info().set_index('uid').loc[donors]
    age_info = info['age'].describe()

    dinfo = dict(n_donors=len(info), n_female=sum(info['sex'] == 'F'))
    dinfo.update(age_info.loc[['min', 'max', 'mean', 'std']].to_dict())

    return dinfo


def _get_norm_procedure(norm, parameter):
    """
    Returns methods description of the provided `norm`

    Parameters
    ----------
    norm : str
        Normalization procedure

    Returns
    -------
    procedure : str
        Reporting procedures
    """

    if parameter not in ('sample_norm', 'gene_norm'):
        raise ValueError(f'Invalid norm parameter {parameter}')

    if parameter == 'sample_norm':
        mod = ('tissue sample expression values across genes')
        suff = ('tissue sample across genes')
    else:
        mod = ('gene expression values across tissue samples')
        suff = ('gene across tissue samples')

    sigmoid = r"""
    independently for each donor using a sigmoid function [F2016P]:<br>

    $$ x_{{norm}} = \frac{{1}}{{1 + \exp(-\frac{{(x - \overline{{x}})}}
    {{\sigma_{{x}}}})}} $$<br>

    where $\bar{{x}}$ is the arithmetic mean and $\sigma_{{x}}$ is the
    sample standard deviation of the expression of a single
    """
    robust_sigmoid = r"""
    using a robust sigmoid function [F2013J]:<br>

    $$ x_{{norm}} = \frac{{1}}{{1 + \exp(-\frac{{(x-\langle x \rangle)}}
    {{\text{{IQR}}_{{x}}}})}} $$<br>

    where $\langle x \rangle$ is the median and $\text{{IQR}}_{{x}}$
    is the normalized interquartile range of the expression of a single
    """
    rescale = r"""
    Normalized expression values were then rescaled to the unit interval: <br>

    $$ x_{{scaled}} = \frac{{x_{{norm}} - \min(x_{{norm}})}}
    {{\max(x_{{norm}}) - \min(x_{{norm}})}} $$<br>
    """

    procedure = ''
    if norm in ['center', 'demean']:
        procedure += """
        by demeaning {mod} independently for each donor.
        """.format(mod=mod)
    elif norm == 'zscore':
        procedure += """
        by mean- and variance-normalizing (i.e., z-scoring) {mod} independently
        for each donor.
        """.format(mod=mod)
    elif norm == 'minmax':
        procedure += r"""
        by rescaling {mod} to the unit interval independently for each donor:
        <br>

        $$ x_{{{{scaled}}}} = \frac{{{{x - \min(x)}}}}{{{{\max(x) - \min(x)}}}}
        $$<br>
        """.format(mod=mod)
    elif norm in ['sigmoid', 'sig']:
        procedure += r"""
        by normalizing {mod} {sigmoid} {suff}.
        """.format(mod=mod, sigmoid=sigmoid, suff=suff)
    elif norm in ['scaled_sigmoid', 'scaled_sig']:
        procedure += r"""
        by normalizing {mod} {sigmoid} {suff}. {rescale}
        """.format(mod=mod, sigmoid=sigmoid, suff=suff, rescale=rescale)
    elif norm in ['scaled_sigmoid_quantiles', 'scaled_sig_qnt']:
        procedure += r"""
        by normalizing {mod} {sigmoid} {suff}, calculated using only data in
        the 5–95th percentile range to downweight the impact of outliers.
        {rescale}
        """.format(mod=mod, sigmoid=sigmoid, suff=suff, rescale=rescale)
    elif norm in ['robust_sigmoid', 'rsig', 'rs']:
        procedure += r"""
        by normalizing {mod} {robust_sigmoid} {suff}.
        """.format(mod=mod, robust_sigmoid=robust_sigmoid, suff=suff)
    elif norm in ['scaled_robust_sigmoid', 'scaled_rsig', 'srs']:
        procedure += r"""
        by normalizing {mod} {robust_sigmoid} {suff}. {rescale}
        """.format(mod=mod, robust_sigmoid=robust_sigmoid, suff=suff,
                   rescale=rescale)
    elif norm in ['mixed_sigmoid', 'mixed_sig']:
        procedure += r"""
        by normalizing {mod} using a mixed sigmoid function [F2013J]:<br>

        $$ x_{{{{norm}}}} = \left\{{{{\begin{{{{array}}}}{{{{r r}}}}
        \frac{{{{1}}}}{{{{1 + \exp(-\frac{{{{(x-\overline{{{{x}}}})}}}}
        {{{{\sigma_{{{{x}}}}}}}})}}}} ,& \text{{{{IQR}}}}_{{{{x}}}} = 0
        \frac{{{{1}}}}{{{{1 + \exp(-\frac{{{{(x-\langle x \rangle)}}}}
        {{{{\text{{{{IQR}}}}_{{{{x}}}}}}}})}}}} ,& \text{{{{IQR}}}}_{{{{x}}}}
        \neq 0 \end{{{{array}}}}\right. $$<br>

        where $\bar{{{{x}}}}$ is the arithmetic mean, $\sigma_{{{{x}}}}$ is the
        sample standard deviation, $\langle x \rangle$ is the median, and
        $\text{{{{IQR}}}}_{{{{x}}}}$ is the normalized interquartile range of
        the expression value of a single {suff}. {rescale}
        """.format(mod=mod, suff=suff, rescale=rescale)

    return procedure


def _add_references(report):
    """
    Detects references in `report` and generates list

    Parameters
    ----------
    report : str
        Report body

    Returns
    -------
    references : str
        List of references to be appended to `report`
    """

    refreport = ''
    for ref, cite in REFERENCES.items():
        if ref in report:
            refreport += f'[{ref}]: {cite}<p>'

    if len(refreport) > 0:
        refreport = '<br> REFERENCES<p>----------<p>' + refreport

    if refreport.endswith('<p> '):
        refreport = refreport[:-4]

    return refreport
