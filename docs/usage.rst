.. _usage:

Using abagen
============

Downloading the AHBA data
-------------------------

In order to use ``abagen``, you'll need to download the `AHBA microarray data
<http://human.brain-map.org/static/download>`_. You can download it with the
following command:

.. code-block:: python

    >>> import abagen
    >>> files = abagen.fetch_microarray(donors='all')

.. note::

    Downloading the entire dataset (about 4GB) can take a long time!

This command will download data from all six available donors into a folder
called ``allenbrain`` in the current directory. If you have already downloaded
the data, you can provide the ``data_dir`` argument to specify where the files
have been stored:

.. code-block:: python

    >>> files = abagen.fetch_microarray(data_dir='/path/to/local/download', donors='all')

The returned object ``files`` is a dictionary with filepaths to the five
different file types in the AHBA dataset:

.. code-block:: python

    >>> print(files.keys())
    dict_keys(['microarray', 'ontology', 'pacall', 'probes', 'annotation'])

We can take a look at the data in these files using the ``abagen.io``
functions. There are ``io`` functions for each of the five file types; you can
get more information on the functions and the data contained in each file type
by looking at the reference :ref:`api_ref`. Notably, all ``io`` functions
return ``pandas.DataFrame`` objects for ease-of-use:

.. code-block:: python

    >>> abagen.io.read_annotation(files.annotation[0]).head()
               structure_id  slab_num  well_id  ...  mni_x mni_y mni_z
    sample_id                                   ...
    0                  4077        22      594  ...    5.9 -27.7  49.7
    1                  4323        11     2985  ...   29.2  17.0  -2.9
    2                  4323        18     2801  ...   28.2 -22.8  16.8
    3                  4440        18     2273  ...  -24.6 -24.6   1.3
    4                  4266        17     2785  ...   31.1 -31.3  -7.3

    [5 rows x 13 columns]
    >>> abagen.io.read_probs(files.probes[0]).head()
                          probe_name  gene_id    ...     entrez_id chromosome
    probe_id                                     ...
    1058685              A_23_P20713      729    ...         733.0          9
    1058684   CUST_15185_PI416261804      731    ...         735.0          5
    1058683             A_32_P203917      731    ...         735.0          5
    1058682             A_23_P138819      736    ...         740.0         11
    1058681             A_24_P232500      736    ...         740.0         11

    [5 rows x 6 columns]

Once you've downloaded the data you'll need to parcellate it.

Defining a parcellation
-----------------------

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
    >>> atlas.image
    '/local/path/to/atlas-desikankilliany.nii.gz'
    >>> atlas.info
    '/local/path/to/atlas-desikankilliany.csv'

While only the atlas image (i.e., Nifti file) is required for processing the
microarray data, the CSV can also be very useful; we can use the CSV file to
constrain the matching of tissue samples to anatomical regions in our atlas.

If you want to supply your own CSV file with information about an atlas, you
must ensure it has the following columns:

  1. ``id``: an integer ID corresponding to the labels in the ``atlas`` image
  2. ``hemisphere``: a L/R hemispheric designation (i.e., 'L' or 'R')
  3. ``structure``: a broad structural class designation (i.e., one of
     'cortex', 'subcortex', or 'cerebellum')

For example, a valid CSV might look like this:

.. code-block:: python

    >>> import pandas as pd
    >>> pd.read_csv(atlas.info).head()
       id                    label hemisphere structure
    0   1  lateralorbitofrontal_rh          R    cortex
    1   2         parsorbitalis_rh          R    cortex
    2   3           frontalpole_rh          R    cortex
    3   4   medialorbitofrontal_rh          R    cortex
    4   5      parstriangularis_rh          R    cortex

Notice that extra columns (i.e., ``label``) are okay, as long as the three
required columns are present!

Getting expression data
-----------------------

Now that the microarray data have been downloaded and we have a parcellation,
we can process the data. This is as simple as:

.. code-block:: python

    >>> expression = abagen.get_expression_data(files, atlas.image, atlas.info)

.. note::

    Wrangling all the raw microarray data can be quite time-consuming! If you'd
    like to speed up this step you can make sure you've performed the
    :ref:`io_installation`.

The ``expression`` object returned is a ``pandas.DataFrame``, where rows
correspond to region labels as defined in our atlas image, columns correspond
to genes, and entry values are normalized microarray expression data averaged
across donors:

.. code-block:: python

    >>> expression.head()
    gene_symbol    MRPL49    ZNHIT2     ...       A_32_P9207  A_32_P94122
    label                               ...
    1            0.407088  0.478699     ...         0.305448     0.470933
    2            0.391223  0.636014     ...         0.383983     0.585307
    3                 NaN       NaN     ...              NaN          NaN
    4            0.492941  0.373068     ...         0.364473     0.246995
    5            0.358736  0.241114     ...         0.250388     0.215016

    [5 rows x 20597 columns]

Unfortunately, due to how tissue samples were collected from the donor brains,
it is possible that some regions in an atlas may not be represented by any
expression data. As you can see above, the third row is filled with NaN values.
That region, corresponding to the right frontal pole in the Desikan-Killiany
atlas, was not matched to any tissue samples; this is likely due to the fact
that only two of the six donors had any tissue samples taken from the right
hemisphere.

If you require a full matrix with expression data for *every* region, you can
specify the following:

.. code-block:: python

    >>> expression = abagen.get_expression_data(files, atlas.image, atlas.info, exact=False)
    >>> expression.head()
    gene_symbol    MRPL49    ZNHIT2     ...       A_32_P9207  A_32_P94122
    label                               ...
    1            0.408125  0.488091     ...         0.307095     0.480003
    2            0.392768  0.644425     ...         0.386279     0.591653
    3            0.507654  0.000000     ...         0.257872     0.201342
    4            0.494560  0.388589     ...         0.366562     0.257212
    5            0.359371  0.258268     ...         0.251915     0.225289

    [5 rows x 20597 columns]

By default, ``get_expression_data()`` will attempt to be as precise as possible
in matching microarray samples with brain regions. Specifying ``exact=False``
will, at the cost of this precision, ensure that every brain region is matched
to *at least* one sample.

You can investigate other options for modifying how the ``expression`` array is
generated by looking at the :ref:`api_ref`.
