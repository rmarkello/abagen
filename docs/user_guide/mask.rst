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
    >>> atlas = abagen.fetch_desikan_killiany()
    >>> dk = nib.load(atlas['image'])
    >>> phg = dk.__class__(dk.dataobj[:] == 67, dk.affine, dk.header)

We can then use this mask to obtain all the microarray samples that fall within
its boundaries:

.. doctest::

    >>> exp, coords = abagen.get_samples_in_mask(mask=phg)

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

    >>> exp.head()
    gene_symbol      A1BG  A1BG-AS1       A2M  ...       ZYX     ZZEF1      ZZZ3
    well_id                                    ...
    2850         0.679381  0.069934  0.398107  ...  0.000000  0.000000  0.000000
    998          0.354435  0.441267  0.752046  ...  0.487668  0.091851  0.681403
    990          0.318782  0.531280  0.958299  ...  0.536189  0.187061  0.695825
    ...               ...       ...       ...  ...       ...       ...       ...
    141667159    1.000000  1.000000  0.476643  ...  0.840348  0.175459  0.445667
    159226157    1.000000  0.000000  1.000000  ...  0.000000  0.000000  0.000000
    159226117    0.000000  1.000000  0.000000  ...  1.000000  1.000000  1.000000
    <BLANKLINE>
    [40 rows x 15633 columns]

This allows you to match up the samples with additional data provided by
the AHBA (e.g., ontological information) as desired.

.. _usage_mask_all:

Get ALL the samples
-------------------

If you want all of the available processed samples rather than only those
within a given mask you can call the function without providing an explicit
mask:

.. doctest::
    :options: +SKIP

    >>> exp, coords = abagen.get_samples_in_mask(mask=None)

This will return all samples (after dropping those where the listed MNI
coordinates don't match the listed hemisphere designation, etc.).
