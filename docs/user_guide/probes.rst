.. _usage_probe_selection:

Probe selection options
=======================

The probes used to measure microarray expression levels in the AHBA data are
often redundant; that is, there are frequently several probes indexing the same
gene. Since the output of the :func:`abagen.get_expression_data` workflow is a
region by gene dataframe, at some point we need to transition from indexing
probe expression level to indexing gene expression levels. Effectively, this
means we need to condense or select from the redundant probes for each gene;
unfortunately, there are a number of ways to do that.

Currently, ``abagen`` supports seven options for this probe to gene conversion.
All the options have been used at various points throughout the published
record, so while there is no "right" answer we do encourage using the default
option (`Differential stability`_) due to recent work by Arnatkevičiūte et al.,
2019 showing that it provides the highest fidelity to RNA sequencing data.

Nonetheless, we describe all the methods in detail here; these can be
implemented by passing the ``probe_selection`` keyword argument to
:func:`abagen.get_expression_data`. For a selection of references to published
works that have used these different methods please see the documentation of
:func:`abagen.probes.collapse_probes`.

Average
-------
Takes the average expression values for all probes indexing the same gene.

.. image:: imgs/average.png
   :align: center

.. code-block:: python

    >>> abagen.get_expression_data(atlas['image'], probe_selection='average')

Max intensity
-------------
Selects the probe with the highest average expression across all samples (where
samples are concatenated across donors).

.. image:: imgs/max_intensity.png
   :align: center

.. code-block:: python

    >>> abagen.get_expression_data(atlas['image'], probe_selection='max_intensity')


Max variance
------------
Selects the probe with the highest variance in expression across all samples
(where samples are concatenated across donors).

.. image:: imgs/max_variance.png
   :align: center

.. code-block:: python

    >>> abagen.get_expression_data(atlas['image'], probe_selection='max_variance')


Principal component loading
---------------------------
Selects the probe with the highest loading on the first principal component
derived from the probe microarray expression across all samples (where samples
are concatenated across donors).

.. image:: imgs/pc_loading.png
   :align: center

.. code-block:: python

    >>> abagen.get_expression_data(atlas['image'], probe_selection='pc_loading')


Correlation
-----------
When there are more than two probes indexing the same gene, selects the probe
with the highest average correlation to other probes across all samples (where
samples are concatenated across donors).

.. image:: imgs/correlation.png
   :align: center

When there are exactly two probes the correlation procedure cannot be used, and
so you can fall back to either the `Max intensity`_ (``corr_intensity``) or
the `Max variance`_ (``corr_variance``) criteria.

.. code-block:: python

    >>> abagen.get_expression_data(atlas['image'], probe_selection='corr_intensity')
    >>> abagen.get_expression_data(atlas['image'], probe_selection='corr_variance')


Differential stability
----------------------
Computes the Spearman correlation of microarray expression values for each
probe across brain regions for every **pair** of donors. Correlations are
averaged and the probe with the highest correlation is retained.

.. image:: imgs/diff_stability.png
   :align: center

.. code-block:: python

    >>> abagen.get_expression_data(atlas['image'], probe_selection='diff_stability')
