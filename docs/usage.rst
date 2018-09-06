.. _usage:

Example usage
=============

At a minimum, you will need an atlas image that you want to use to parcellate
the microarray expression data. The supplied atlas should be a Nifti image
where regions or parcels are denoted by unique integer IDs. You can use
whatever atlas you would like, but ``abagen`` ships with a volume-based version
of the `Desikan-Killiany atlas <https://surfer.nmr.mgh.harvard.edu/ftp/
articles/desikan06-parcellation.pdf>`_ that we'll use for demonstration
purposes:

.. code-block:: python

   >>> import abagen
   >>> atlas = abagen.fetch_desikan_killiany()

The returned object ``atlas`` is a ``dict`` with two keys: ``image``, which
points to the Nifti file containing the atlas data, and ``info``, which points
to a CSV file containing auxilliary information about the parcellation.

.. code-block:: python

    >>> atlas.image
    '/local/path/to/atlas-desikankilliany.nii.gz'
    >>> atlas.info
    '/local/path/to/atlas-desikankilliany.csv'

While we only **need** the image, the info CSV can also be very useful; we can
use the CSV file to constrain the matching of microarray expression samples to
anatomical regions.

If you want to supply your own CSV file with information about an atlas, you
must ensure it has the following columns:

  1. ``id``: ID corresponding to the integer labels in the ``atlas`` image
  2. ``hemisphere``: L/R hemispheric designation
  3. ``structure``: broad structural class (i.e., 'cortex', 'subcortex', or
     'cerebellum')

For example, a valid CSV might look like this:

.. code-block:: python

    >>> import pandas as pd
    >>> pd.read_csv(atlas.info)
       id                    label hemisphere structure
    0   1  lateralorbitofrontal_rh          R    cortex
    1   2         parsorbitalis_rh          R    cortex
    2   3           frontalpole_rh          R    cortex
    3   4   medialorbitofrontal_rh          R    cortex
    4   5      parstriangularis_rh          R    cortex
    ...

Notice that extra columns (i.e., ``label``) are okay, as long as the three
required columns as present.

Once you have an atlas and, if desired, a supplementary CSV file, you can
download the AHBA data:

.. code-block:: python

    >>> files = abagen.fetch_microarray(donors='all')

.. note::

    Downloading all the data (approximately 4GB) will take a long time!
    Thankfully, you only have to do this step once.

If you do not specify ``donors='all'``, microarray expression data from only
one donor will be downloaded. If you have already downloaded the microarray
expression from the Allen Brain Institute website, you can set the ``data_dir``
argument to use those files instead of re-downloading it:

.. code-block:: python

    >>> files = abagen.fetch_microarray(data_dir='/local/path/to/download', donors='all')


The returned object ``files`` is a dictionary pointing to the different data
files supplied by the Allen Institute. We can then use those files, along with
the atlas image and atlas info, to generate a region x gene expression array:

.. code-block:: python

    >>> expression = abagen.get_expression_data(files, atlas.image, atlas.info)

.. note::

    Wrangling all the raw microarray data is quite time-consuming! This call
    will take quite a while. You can speed it up by making sure that you've
    performed the :ref:`full_installation`!

Unfortunately, due to how samples were collected from the donor brains, it is
possible that some regions in the atlas may not be represented by any
expression data. If you require a full matrix with expression data for *every*
region, you can specify the following:

.. code-block:: python

    >>> expression = abagen.get_expression_data(files, atlas.image, atlas.info, exact=False)


By default, ``abagen`` will attempt to be as precise as possible in matching
microarray samples with brain regions. Specifying ``exact=False`` will, at the
cost of this precision, ensure that every brain region is matched to *at least*
one sample. You can investigate other options for modifying how the
``expression`` array is generated in the documentation by typing
``help(abagen.get_expression_data)``.
