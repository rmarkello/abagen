.. _usage_reporting:

Generating reporting methods
============================

There are a lot of options and parameters to choose from when processing data
with :func:`abagen.get_expression_data` and while we've attempted to select
reasonable defaults, we don't want to limit your options—you are free to pick
and choose any combination of inputs to process the AHBA data! That said, we
also wanted to make it easy for you to report *exactly* what was done to the
data based on the parameters you choose, so to that end we have added the
``return_report`` parameter to :func:`abagen.get_expression_data`.

When ``return_report=True``, the workflow will return an extra output. That is,
in addition to the default regional microarray expression dataframe, a string
will be returned that describes, in detail, all the processing that was done to
the AHBA in the process of generating the expression matrix. We have tried to
write this in such a way that you can simply copy-and-paste the provided text
into the methods section of a paper, though you are of course free to edit it
as you see fit (though if you feel edits are necessary please let us know and
we can modify the generation more permanently!).

Example report
--------------

A report can be generated with:

.. code-block::

    >>> expression, report = abagen.get_expression_data(atlas['image'], atlas['info'],
    ...                                                 return_report=True)

Alternatively, you can use the ``abagen.reporting`` module to generate a report
directly without having to re-run the entire pipeline. (Note that the ``Report``
class accepts (nearly) all the same parameters as the ``get_expression_data()``
workflow.)

.. doctest::

    >>> from abagen import reporting
    >>> generator = reporting.Report(atlas['image'], atlas['info'])
    >>> report = generator.gen_report()

The returned ``report`` (with default parameters) will look something like this
example:

    *Regional microarry expression data were obtained from 6 post-mortem brains
    (1 female, ages 24.0--57.0, 42.50 +/- 13.38) provided by the Allen Human
    Brain Atlas (AHBA, https://human.brain-map.org; [H2012N]). Data were
    processed with the abagen toolbox (version 1.0;
    https://github.com/rmarkello/abagen) using a 83-region volumetric atlas
    in MNI space.*

    *First, microarray probes were reannotated using data provided by [A2019N];
    probes not matched to a valid Entrez ID were discarded. Next, probes were
    filtered based on their expression intensity relative to background noise
    [Q2002N], such that probes with intensity less than the background in
    >=50.00% of samples across donors were discarded, yielding 31,569 probes.
    When multiple probes indexed the expression of the same gene, we selected
    and used the probe with the most consistent pattern of regional variation
    across donors (i.e., differential stability; [H2015N]), calculated with:*

    :math:`\Delta_{S}(p) = \frac{1}{\binom{N}{2}} \, \sum_{i=1}^{N-1} \sum_{j=i+1}^{N} \rho[B_{i}(p), B_{j}(p)]`

    *where* :math:`\rho` *is Spearman's rank correlation of the expression of a
    single probe, p, across regions in two donor brains* :math:`B_{i}` *and*
    :math:`B_{j}`, *and N is the total number of donors. Here, regions
    correspond to the structural designations provided in the ontology from the
    AHBA.*

    *The MNI coordinates of tissue samples were updated to those generated via
    non-linear registration using the Advanced Normalization Tools (ANTs;
    https://github.com/chrisfilo/alleninf). Samples were assigned to brain
    regions in the provided atlas if their MNI coordinates were within 2 mm of
    a given parcel. To reduce the potential for misassignment, sample-to-region
    matching was constrained by hemisphere and gross structural divisions
    (i.e., cortex, subcortex/brainstem, and cerebellum, such that e.g., a
    sample in the left cortex could only be assigned to an atlas parcel in the
    left cortex; [A2019N]). All tissue samples not assigned to a brain region
    in the provided atlas were discarded.*

    *Inter-subject variation was addressed by normalizing tissue sample
    expression values across genes using a robust sigmoid function [F2013J]:*

    :math:`x_{norm} = \frac{1}{1 + \exp(-\frac{(x-\langle x \rangle)} {\text{IQR}_{x}})}`

    *where* :math:`\langle x \rangle` *is the median and* :math:`\text{IQR}_{x}`
    *is the normalized interquartile range of the expression of a single tissue
    sample across genes. Normalized expression values were then rescaled to the
    unit interval:*

    :math:`x_{scaled} = \frac{x_{norm} - \min(x_{norm})} {\max(x_{norm}) - \min(x_{norm})}`

    *Gene expression values were then normalized across tissue samples using an
    identical procedure. Samples assigned to the same brain region were
    averaged separately for each donor and then across donors, yielding a
    regional expression matrix with 83 rows, corresponding to brain regions,
    and 15,633 columns, corresponding to the retained genes.*

    *REFERENCES* |br|
    -\-\-\-\-\-\-\-\-\- |br|
    *[A2019N]: Arnatkevic̆iūtė, A., Fulcher, B. D., & Fornito, A. (2019). A
    practical guide to linking brain-wide gene expression and neuroimaging
    data. Neuroimage, 189, 353-367.* |br|
    *[F2013J]: Fulcher, B. D., Little, M. A., & Jones, N. S. (2013). Highly
    comparative time-series analysis: the empirical structure of time series
    and their methods. Journal of the Royal Society Interface, 10(83),
    20130048.* |br|
    *[H2012N]: Hawrylycz, M. J., Lein, E. S., Guillozet-Bongaarts, A. L., Shen,
    E. H., Ng, L., Miller, J. A., ... & Jones, A. R. (2012). An anatomically
    comprehensive atlas of the adult human brain transcriptome. Nature,
    489(7416), 391-399.* |br|
    *[H2015N]: Hawrylycz, M., Miller, J. A., Menon, V., Feng, D., Dolbeare,
    T., Guillozet-Bongaarts, A. L., ... & Lein, E. (2015). Canonical genetic
    signatures of the adult human brain. Nature Neuroscience, 18(12), 1832.*
    |br|
    *[Q2002N]: Quackenbush, J. (2002). Microarray data normalization and
    transformation. Nature Genetics, 32(4), 496-501.*

Note that due to text formatting limitations in Python, relevant equations used
for e.g., normalizing the expression data will be provided in LaTeX format
(i.e., surrounded by `$$` characters and with TeX math commands).

.. important::
    Please note that we explicitly release all text in the ``abagen.reporting``
    `module <ref_reporting>`_ (used to generate the above-referenced reports)
    under a `CC0 license <https://creativecommons.org/publicdomain/zero/1.0/>`_
    such that it can be used in manuscripts without modification.

.. |br| raw:: html

   <br>
