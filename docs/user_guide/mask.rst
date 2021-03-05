.. _usage_mask:

Using a binary mask
===================

.. _usage_mask_basic:

Basic usage
-----------

Sometimes, you're not interested in aggregating microarray expression samples
within regions of an atlas—you want the actual, sample-level data instead. In
this case, we provides the :func:`abagen.get_samples_in_mask` function.

To demonstrate how this works, we'll first make a brainmask for the left
parahippocampal gyrus using the region definition from the Desikan-Killiany
atlas:

.. doctest::

    >>> import nibabel as nib
    >>> import pandas as pd
    >>> atlas = abagen.fetch_desikan_killiany()
    >>> dk = nib.load(atlas['image'])
    >>> info = pd.read_csv(atlas['info'])
    >>> phg = int(info.query('label == "parahippocampal" & hemisphere == "L"')['id'])
    >>> img = dk.__class__(dk.dataobj[:] == phg, dk.affine, dk.header)

We can then use this mask to obtain all the microarray samples that fall within
its boundaries:

.. doctest::

    >>> expression, coords = abagen.get_samples_in_mask(mask=img)

:func:`abagen.get_samples_in_mask` returns two objects: (1) the the samples x
gene expression matrix (``exp``), and (2) an array of MNI coordinates for those
samples (``coords``). Because this is using :func:`abagen.get_expression_data`
under the hood, the returned expression data have been preprocessed (i.e.,
filtered, normalized) according to that workflow. As such, you can provide all
the same parameters and keyword arguments to :func:`abagen.get_samples_in_mask`
as you can to :func:`abagen.get_expression_data` (with the exception of
``atlas`` which is superseded by ``mask`` and ``region_agg``/``agg_metric``
which will be ignored). Refer to the :ref:`API documentation <api_ref>`) for
more details!

Since the returned expression dataframe is a samples x gene matrix (rather than
regions x gene), the index of the dataframe corresponds to the unique well ID
of the relevant sample (rather than the atlas region):

.. doctest::

    >>> print(expression)
    gene_symbol      A1BG  A1BG-AS1       A2M  ...       ZYX     ZZEF1      ZZZ3
    well_id                                    ...
    2850         0.654914  0.234039  0.283280  ...  0.020379  0.228080  0.000000
    998          0.428705  0.375819  0.457741  ...  0.254195  0.315383  0.502122
    990          0.400673  0.409852  0.561666  ...  0.270064  0.397740  0.522261
    ...               ...       ...       ...  ...       ...       ...       ...
    159226055    0.418706  0.751837  0.087808  ...  0.651541  0.410095  0.462773
    159226117    0.533079  0.773214  0.265615  ...  0.441826  0.389615  0.455249
    158158343    0.362038  0.553050  0.314730  ...  0.346605  0.261426  0.337738
    <BLANKLINE>
    [38 rows x 15633 columns]

This allows you to match up the samples with additional data provided by
the AHBA (e.g., ontological information) as desired.

.. _usage_mask_all:

Get ALL the samples
-------------------

If you want all of the available processed samples rather than only those
within a given mask you can call the function without providing an explicit
mask (this is the default when no ``mask`` parameter is passed):

.. doctest::
    :options: +SKIP

    >>> expression, coords = abagen.get_samples_in_mask(mask=None)

This will return all samples (after dropping those where the listed MNI
coordinates don't match the listed hemisphere designation, etc.).
