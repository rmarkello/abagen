.. _usage_download:

The Allen Human Brain Atlas dataset
===================================

.. _usage_download_fetching:

Fetching the AHBA data
----------------------

In order to use ``abagen``, you'll need to download the `AHBA microarray data
<https://human.brain-map.org/static/download>`_. You can download it with the
following command:

.. doctest::

    >>> import abagen
    >>> files = abagen.fetch_microarray(donors='all')

.. note::

    Downloading the entire dataset (about 4GB) can take a long time depending
    on your internet connection speed! If you don't want to download all the
    donors you can provide the subject IDs of the donors you want as a list
    (e.g., ``['9861', '10021']``) instead of passing ``'all'``.

This command will download data from the specified donors into a folder called
``allenbrain`` in the ``$HOME/abagen-data`` directory. If you have already
downloaded the data you can provide the ``data_dir`` argument to specify where
the files have been stored:

.. doctest::

    >>> files = abagen.fetch_microarray(donors=['12876', '15496'], data_dir='/path/to/my/data/')  # doctest: +SKIP

Alternatively, ``abagen`` will check the directory specified by the
environmental variable ``$ABAGEN_DATA`` and use that as the download location
if the dataset does not already exist there.

If you provide a path to ``data_dir`` (or specify a path with ``$ABAGEN_DATA``)
the directory specified should have the following structure:

.. code-block::

    /path/to/my/data/
    ├── normalized_microarray_donor10021/
    │   ├── MicroarrayExpression.csv
    │   ├── Ontology.csv
    │   ├── PACall.csv
    │   ├── Probes.csv
    │   └── SampleAnnot.csv
    ├── normalized_microarray_donor12876/
    ├── normalized_microarray_donor14380/
    ├── normalized_microarray_donor15496/
    ├── normalized_microarray_donor15697/
    └── normalized_microarray_donor9861/

(Note the directory does not have to be named ``allenbrain`` for this to work.)

.. _usage_download_loading:

Loading the AHBA data
---------------------

The ``files`` object returned by :func:`abagen.fetch_microarray` is a
dictionary with filepaths to the five different file types in the AHBA
microarray dataset:

.. doctest::

    >>> print(files.keys())
    dict_keys(['microarray', 'ontology', 'pacall', 'probes', 'annotation'])

You can load the data in these files using the :mod:`abagen.io` functions.
There are IO functions for each of the five file types; you can get more
information on the functions and the data contained in each file type by
looking at the ref:`api_ref`. Notably, all IO functions return
``pandas.DataFrame`` objects for ease-of-use.

For example, you can load the annotation file for the first donor with:

.. doctest::

    >>> annotation = abagen.io.read_annotation(files['annotation'][0])
    >>> print(annotation)
               structure_id  slab_num  well_id  ... mni_x mni_y mni_z
    sample_id                                   ...
    0                  4077        22      594  ...   5.9 -27.7  49.7
    1                  4323        11     2985  ...  29.2  17.0  -2.9
    2                  4323        18     2801  ...  28.2 -22.8  16.8
    ...                 ...       ...      ...  ...   ...   ...   ...
    943                4758        67     1074  ...   7.9 -72.3 -40.6
    944                4760        67     1058  ...   8.3 -57.4 -59.0
    945                4761        67     1145  ...   9.6 -46.7 -47.6
    <BLANKLINE>
    [946 rows x 13 columns]


And you can do the same for, e.g., the probe file with:

.. doctest::

    >>> probes = abagen.io.read_probes(files['probes'][0])
    >>> print(probes)
                          probe_name  gene_id  ... entrez_id chromosome
    probe_id                                   ...
    1058685              A_23_P20713      729  ...     733.0          9
    1058684   CUST_15185_PI416261804      731  ...     735.0          5
    1058683             A_32_P203917      731  ...     735.0          5
    ...                          ...      ...  ...       ...        ...
    1071209             A_32_P885445  1012197  ...       NaN        NaN
    1071210               A_32_P9207  1012198  ...       NaN        NaN
    1071211              A_32_P94122  1012199  ...       NaN        NaN
    <BLANKLINE>
    [58692 rows x 6 columns]


The other IO functions work similarly for the remaining filetypes.
