.. _usage_donor_norm:

Donor normalization options
===========================

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
there are a number of ways to do that.

Currently, ``abagen`` supports three options for normalizing donor data:

    1. :ref:`usage_donors_srs`,
    2. :ref:`usage_donors_zscore`, and
    3. :ref:`usage_donors_batch`

All the options have been used at various points throughout the published
record, so while there is no "right" choice we do encourage using the default
option (:ref:`scaled robust sigmoid <usage_donors_srs>`) due to recent work by
Arnatkevičiūte et al., 2019 showing that it is---as the name would
suggest---robust to outlier effects commonly observed in microarray data.

We describe all the methods in detail here; these can be implemented by passing
the ``donor_norm`` keyword argument to :func:`abagen.get_expression_data`. For
a selection of references to published works that have used these different
methods please see the documentation of :func:`abagen.normalize_expression`.

.. _usage_donors_srs:

Scaled robust sigmoid
---------------------

.. code-block:: python

    >>> abagen.get_expression_data(atlas['image'], donor_norm='srs')

Microarray values are separately normalized for each gene across brain regions
using a robust sigmoid function:

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

Once the expression values have been normalized they are rescaled to the unit
interval with:

.. math::

   x_{norm} = \frac{x_{y} - \text{min}(x)}
                   {\text{max}(x) - \text{min}(x)}

Normalization is performed separately for each donor.

.. _usage_donors_zscore:

Z-score
-------

.. code-block:: python

    >>> abagen.get_expression_data(atlas['image'], donor_norm='zscore')

Microarray values are separately normalized for each gene across brain regions
using a basic z-score function:

.. math::

    x_{norm} = \frac{x_{y} - \bar{x}}
                    {\sigma_{x}}

where :math:`\bar{x}` is the mean and :math:`\sigma_{x}` is the sample standard
deviation of the microarray expression values.

Normalization is performed separately for each donor.

.. _usage_donors_batch:

Batch correction
----------------

.. code-block:: python

    >>> abagen.get_expression_data(atlas['image'], donor_norm='batch')

Region by gene expression matrices for each donor are vertically concatenated
(across donors) and donor-specific indicator variables are fit to the resulting
expression data matrix with a simple linear regression. Beta estimates for
donor effects are estimated independently for each gene:

.. math::

    x = \beta_{0} + \beta_{1} I_{1} + \beta_{1} I_{1} + \ldots + \beta_{n} I_{n} \epsilon

where :math:`\beta_{0}` is the intercept and :math:`I_{1}` is the indicator
variable for a given donor. Concatenated expression data are residualized based
on the regression fit and then unstacked into individual donor expression
matrices:

.. math::

   x_{norm} = x - (\beta_{1} I_{1} + \beta_{2} I_{2} + \ldots + \beta_{n} I_{n})

Note that the linear model fit includes the intercept but the intercept is not
removed during the residualization process.

Normalization is performed simultaneously for all donors.

No normalization
----------------

.. code-block:: python

    >>> abagen.get_expression_data(atlas['image'], donor_norm=None)

Providing ``None`` to the ``donor_norm`` parameter will prevent any
normalization procedure from being performed on the data. If you use this it
is **strongly** encouraged that you also specify ``return_donors=True`` so that
expression data are not aggregated across donors.
