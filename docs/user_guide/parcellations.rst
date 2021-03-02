.. _usage_parcellations:

Defining a parcellation
=======================

.. _usage_parcellations_acceptable:

Acceptable parcellations
------------------------

In order to process the microarray expression data from AHBA you'll need a
parcellation (or "atlas"). Here, we define a parcellation (atlas) as any image
in `MNI space`_ that contains regions or parcels denoted by unique integer IDs.
You can use whatever parcellation/atlas you'd like as long as it is a
volumetric image in MNI (ICBM152) space.

.. note::

    We hope to support surface parcellations in the not-so-distant future! See
    this `GitHub issue <https://github.com/rmarkello/abagen/issues/50>`__ for
    more information.

For demonstration purposes, ``abagen`` has a copy of the `Desikan-Killiany
atlas <https://surfer.nmr.mgh.harvard.edu/ftp/articles/desikan06-parcellation.
pdf>`_ that you can use:

.. doctest::

   >>> atlas = abagen.fetch_desikan_killiany()

Here, the returned object ``atlas`` is a dictionary with two keys: ``image``,
which is filepath to a Nifti image containing the atlas data, and ``info``,
which is a filepath to a CSV file containing auxilliary information about the
parcellation:

.. doctest::

    >>> print(atlas['image'])  # doctest: +ELLIPSIS
    /.../data/atlas-desikankilliany.nii.gz
    >>> print(atlas['info'])  # doctest: +ELLIPSIS
    /.../data/atlas-desikankilliany.csv.gz

.. _usage_parcellations_additional:

Providing additional parcellation info
--------------------------------------

While only the image (i.e., Nifti file) is required for processing the
microarray data, the CSV file with information on the parcellation scheme can
also be very useful. In particular, ``abagen`` can use the CSV file to
constrain the matching of tissue samples to anatomical regions in the atlas
image.

If you want to supply your own CSV file with information about an atlas you
must ensure it has the following columns:

  1. ``id``: an integer ID corresponding to the labels in the ``atlas`` image
  2. ``hemisphere``: a L/R hemispheric designation (i.e., 'L' or 'R')
  3. ``structure``: a broad structural class designation (i.e., one of
     'cortex', 'subcortex', 'cerebellum', 'brainstem', 'white matter', or
     'other')

For example, a valid CSV might look like this:

.. doctest::

    >>> import pandas as pd
    >>> atlas_info = pd.read_csv(atlas['info'])
    >>> print(atlas_info)
        id                    label hemisphere  structure
    0    1  lateralorbitofrontal_rh          R     cortex
    1    2         parsorbitalis_rh          R     cortex
    2    3           frontalpole_rh          R     cortex
    ..  ..                      ...        ...        ...
    80  81           hippocampus_lh          L  subcortex
    81  82              amygdala_lh          L  subcortex
    82  83             brainstem_lh          L  subcortex
    <BLANKLINE>
    [83 rows x 4 columns]


Notice that extra columns (i.e., ``label``) are okay as long as the three
required columns are present! If you want to confirm your file is formatted
correctly you can use :func:`abagen.utils.check_atlas_info`:

.. doctest::

    >>> from abagen import utils
    >>> utils.check_atlas_info(atlas['image'], atlas['info'], validate=True)

If something is amiss with the file this function will raise an error and try
to give some information about what you should check for.

.. important::

    You might be asking: **"why should I provide this extra information for**
    **my parcellation?"** Providing this CSV file will ensure that microarray
    samples designated as belonging to a given hemisphere/structure by the AHBA
    ontology are not matched to regions in the ``atlas`` image with different
    hemispheric/structural designations. That is, if the AHBA ontology
    specifies that a tissue sample comes from the left hemisphere subcortex, it
    will only ever be matched to regions in ``atlas`` belonging to the left
    hemisphere subcortex.

    While this seems trivial, it is **very important** because there are
    numerous tissue samples which occur on the boundaries of hemispheres and
    structural classes (i.e., cortex/subcortex). In many instances, these
    samples won't fall directly within a region of the ``atlas``, at which
    point ``abagen`` will attempt to match them to nearby regions. Without the
    hemisphere/structure information provided by this CSV file there is a high
    likelihood of misassigning samples, leading to biased or skewed expression
    data.

Individualized parcellations
----------------------------

Instead of providing a single parcellation image in MNI space that will be used
for all donors, you can instead provide a parcellation image for each donor in
the space of their "raw" (or native) T1w image. ``abagen`` ships with versions
of the Desikan-Killiany parcellation defined in donor-native space:

.. doctest::

    >>> atlas = abagen.fetch_desikan_killiany(native=True)
    >>> print(atlas['image'].keys())
    dict_keys(['9861', '10021', '12876', '14380', '15496', '15697'])
    >>> print(atlas['image']['9861'])  # doctest: +ELLIPSIS
    /.../data/native_dk/9861/atlas-desikankilliany.nii.gz

Note here that ``atlas['image']`` is a dictionary, where the keys are donor IDs
and the corresponding values are paths to the parcellation for each donor. The
primary workflows in ``abagen`` that accept a single atlas (i.e.,
:func:`abagen.get_expression_data` and :func:`abagen.get_samples_in_mask`) will
also accept a dictionary of this format.

.. note::
    The donor-native versions of the DK parcellation shipped with ``abagen``
    were generated by Arnatkevičiūte et al., 2018, *NeuroImage*, and are
    provided under the CC BY 4.0 license. They have been re-numbered to match
    the MNI-space DK parcellation provided with earlier versions of ``abagen``.

.. _MNI space: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1088516/
