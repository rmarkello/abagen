.. _usage_parcellations:

Defining a parcellation
=======================

In order to process the microarray expression data, you'll need a parcellation.
A parcellation is an atlas image that contains regions or parcels denoted by
unique integer IDs. You can use whatever atlas you'd like as long as it is a
volumetric atlas in MNI space.

For demonstration purposes, ``abagen`` has a copy of the `Desikan-Killiany
atlas <https://surfer.nmr.mgh.harvard.edu/ftp/articles/desikan06-parcellation.
pdf>`_ that you can use:

.. code-block:: python

   >>> atlas = abagen.fetch_desikan_killiany()

Here, the returned object ``atlas`` is a dictionary with two keys: ``image``,
which is a filepath to the Nifti containing the atlas data, and ``info``, which
is a filepath to a CSV file containing auxilliary information about the
parcellation:

.. code-block:: python

    >>> print(atlas.keys())
    dict_keys(['image', 'info'])
    >>> atlas['image']
    '/local/path/to/atlas-desikankilliany.nii.gz'
    >>> atlas['info']
    '/local/path/to/atlas-desikankilliany.csv'

While only the atlas image (i.e., Nifti file) is required for processing the
microarray data, the CSV can also be very useful; we can use the CSV file to
constrain the matching of tissue samples to anatomical regions in our atlas.

If you want to supply your own CSV file with information about an atlas, you
must ensure it has the following columns:

  1. ``id``: an integer ID corresponding to the labels in the ``atlas`` image
  2. ``hemisphere``: a L/R hemispheric designation (i.e., 'L' or 'R')
  3. ``structure``: a broad structural class designation (i.e., one of
     'cortex', 'subcortex', 'cerebellum', 'brainstem', 'white matter', or
     'other')

For example, a valid CSV might look like this:

.. code-block:: python

    >>> import pandas as pd
    >>> pd.read_csv(atlas['info']).head()
       id                    label hemisphere structure
    0   1  lateralorbitofrontal_rh          R    cortex
    1   2         parsorbitalis_rh          R    cortex
    2   3           frontalpole_rh          R    cortex
    3   4   medialorbitofrontal_rh          R    cortex
    4   5      parstriangularis_rh          R    cortex

Notice that extra columns (i.e., ``label``) are okay as long as the three
required columns are present! If you want to confirm your file is formatted
correctly you can use :func:`abagen.utils.check_atlas_info()`.
