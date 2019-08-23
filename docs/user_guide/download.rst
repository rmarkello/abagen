.. _usage_download:

Downloading AHBA data
=====================

In order to use ``abagen``, you'll need to download the `AHBA microarray data
<https://human.brain-map.org/static/download>`_. You can download it with the
following command:

.. code-block:: python

    >>> import abagen
    >>> files = abagen.fetch_microarray(donors='all')

.. note::

    Downloading the entire dataset (about 4GB) can take a long time depending
    on your internet connection speed!

This command will download data from all six available donors into a folder
called ``allenbrain`` in the current directory. If you have already downloaded
the data, you can provide the ``data_dir`` argument to specify where the files
have been stored:

.. code-block:: python

    >>> files = abagen.fetch_microarray(donors='all', data_dir='/path/to/data')

The returned object ``files`` is a dictionary with filepaths to the five
different file types in the AHBA dataset:

.. code-block:: python

    >>> print(files.keys())
    dict_keys(['microarray', 'ontology', 'pacall', 'probes', 'annotation'])

We can take a look at the data in these files using the :mod:`abagen.io`
functions. There are IO functions for each of the five file types; you can get
more information on the functions and the data contained in each file type
by looking at the :ref:`api_ref`. Notably, all IO functions return
:class:`pandas.DataFrame` objects for ease-of-use:

.. code-block:: python

    >>> abagen.io.read_annotation(files['annotation'][0]).head()
               structure_id  slab_num  well_id  ...  mni_x mni_y mni_z
    sample_id                                   ...
    0                  4077        22      594  ...    5.9 -27.7  49.7
    1                  4323        11     2985  ...   29.2  17.0  -2.9
    2                  4323        18     2801  ...   28.2 -22.8  16.8
    3                  4440        18     2273  ...  -24.6 -24.6   1.3
    4                  4266        17     2785  ...   31.1 -31.3  -7.3

    [5 rows x 13 columns]

    >>> abagen.io.read_probes(files['probes'][0]).head()
                          probe_name  gene_id    ...     entrez_id chromosome
    probe_id                                     ...
    1058685              A_23_P20713      729    ...         733.0          9
    1058684   CUST_15185_PI416261804      731    ...         735.0          5
    1058683             A_32_P203917      731    ...         735.0          5
    1058682             A_23_P138819      736    ...         740.0         11
    1058681             A_24_P232500      736    ...         740.0         11

    [5 rows x 6 columns]

Once you've downloaded the data you'll need to parcellate it.
