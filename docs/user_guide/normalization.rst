.. _usage_normalization:

Data normalization options
==========================

The microarray expression data provided by the AHBA has been subjected to some
`normalization procedures <help.brain-map.org/display/humanbrain/
Documentation>`_ designed to mitigate potential differences in expression
values between donors due to "batch effects." Despite these procedures, there
are still some notable differences between donors present in the downloadable
data.

By default, :func:`abagen.get_expression_data` aggregates expression data
across donors (though this can be prevented via the ``return_donors``
parameter). Prior to aggregation, the function performs a within-donor
normalization procedure to attempt to mitigate donor-specific effects; however,
there are a number of ways to achieve this.

Currently, ``abagen`` supports five options for normalizing data:

    1. :ref:`usage_norm_center`,
    2. :ref:`usage_norm_zscore`,
    3. :ref:`usage_norm_minmax`,
    4. :ref:`usage_norm_rs`, and
    5. :ref:`usage_norm_srs`

All the options have been used at various points throughout the published
record, so while there is no "right" choice we do encourage using the default
option (:ref:`scaled robust sigmoid <usage_norm_srs>`) due to recent work by
Arnatkevičiūte et al., 2019 showing that it is---as the name might
suggest---robust to outlier effects commonly observed in microarray data.

We describe all the methods in detail here; these can be implemented by passing
the ``sample_norm`` and ``gene_norm`` keyword arguments to
:func:`abagen.get_expression_data`. For a selection of references to published
works that have used these different methods please see the documentation of
:func:`abagen.normalize_expression`.

.. _usage_norm_sampledonor:

``sample_norm`` vs ``gene_norm``
--------------------------------

Microarray expression data can be normalized in two directions:

    1. Each sample can be normalized across all genes, or
    2. Each gene can be normalized across all samples

These different forms of normalization are controlled by two parameters in the
:func:`abagen.get_expression_data` function: ``sample_norm`` and ``gene_norm``.
Note that normalization of each sample across all genes occurs before
normalization of each gene across all samples.

Both parameters can accept the same five arguments (detailed below), and both
are turned on by default.

.. _usage_norm_center:

Centering
---------

.. code-block:: python

    >>> abagen.get_expression_data(atlas['image'], sample_norm='center', gene_norm='center')

Microarray values are centered with:

.. math::

    x_{norm} = x_{y} - \bar{x}

where :math:`\bar{x}` is the mean of the microarray expression values.

.. _usage_norm_zscore:

Z-score
-------

.. code-block:: python

    >>> abagen.get_expression_data(atlas['image'], sample_norm='zscore', gene_norm='zscore')

Microarray values are normalized using a basic z-score function:

.. math::

    x_{norm} = \frac{x_{y} - \bar{x}}
                    {\sigma_{x}}

where :math:`\bar{x}` is the mean and :math:`\sigma_{x}` is the sample standard
deviation of the microarray expression values.

.. _usage_norm_minmax:

Min-max
-------

.. code-block:: python

    >>> abagen.get_expression_data(atlas['image'], sample_norm='minmax', gene_norm='minmax')

Microarray values are rescaled to the unit interval with:

.. math::

   x_{norm} = \frac{x_{y} - \text{min}(x)}
                   {\text{max}(x) - \text{min}(x)}

.. _usage_norm_rs:

Robust sigmoid
--------------

.. code-block:: python

    >>> abagen.get_expression_data(atlas['image'], sample_norm='rs', gene_norm='rs')

Microarray values are normalized using a robust sigmoid function:

.. math::

   x_{y} = \frac{1}
                {1 + \exp \left( \frac{-(x_{y} - \langle x \rangle)}
                                      {\text{IQR}_{x}}
                          \right)}

where :math:`\langle x \rangle` is the median and :math:`\text{IQR}_{x}` is the
normalized interquartile range of the microarray expression values given as:

.. math::

   \DeclareMathOperator\erf{erf}
   \text{IQR}_{x} = \frac{Q_{3} - Q{1}}
                         {2 \cdot \sqrt{2} \cdot \erf^{-1}\left(\frac{1}{2}\right)}
            \approx \frac{Q_{3} - Q_{1}}
                         {1.35}

.. _usage_norm_srs:

Scaled robust sigmoid
---------------------

.. code-block:: python

    >>> abagen.get_expression_data(atlas['image'], sample_norm='rs', gene_norm='srs')

Microarray values are processed with the :ref:`robust sigmoid <usage_norm_rs>`
function and then rescaled to the unit interval with the :ref:`min-max
<usage_norm_minmax>` function.

.. _usage_norm_none:

No normalization
----------------

.. code-block:: python

    >>> abagen.get_expression_data(atlas['image'], sample_norm=None, gene_norm=None)

Providing ``None`` to the ``sample_norm`` and ``gene_norm`` parameters will
prevent any normalization procedure from being performed on the data. Use this
with caution!
