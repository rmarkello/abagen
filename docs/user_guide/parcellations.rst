.. _usage_parcellations:

Defining a parcellation
=======================

.. _usage_parcellations_acceptable:

Acceptable parcellations
------------------------

In order to process the microarray expression data from AHBA you'll need a
parcellation (or "atlas"). Here, we define a parcellation (atlas) as either (1)
a NIFTI image in `MNI space`_, or (2) a tuple of GIFTI images in `fsaverage
space`_ (and with fsaverage5 resolution!). In both cases, parcels in the atlas
should be denoted by unique integer IDs (distinct across hemispheres). The
primary workflows in ``abagen`` are designed to readily accept any
parcellations / atlases in this format; however, if you want to use a different
format please refer to :ref:`usage_parcellations_special`.

For demonstration purposes, ``abagen`` has a copy of the `Desikan-Killiany
atlas <https://surfer.nmr.mgh.harvard.edu/ftp/articles/desikan06-parcellation.
pdf>`_ that you can use. Here, we load the volumetric atlas by default:

.. doctest::

    >>> import abagen
    >>> atlas = abagen.fetch_desikan_killiany()

The returned object ``atlas`` is a dictionary with two keys: ``image``, which
is filepath to a NIFTI image containing the atlas data, and ``info``, which is
a filepath to a CSV file containing extra information about the parcellation:

.. doctest::

    >>> print(atlas['image'])  # doctest: +ELLIPSIS
    /.../data/atlas-desikankilliany.nii.gz
    >>> print(atlas['info'])  # doctest: +ELLIPSIS
    /.../data/atlas-desikankilliany.csv

You can load the surface version of the atlas by providing the ``surface``
parameter:

.. doctest::

    >>> atlas = abagen.fetch_desikan_killiany(surface=True)
    >>> print(atlas['image'])  # doctest: +ELLIPSIS
    ('/.../data/atlas-desikankilliany-lh.label.gii.gz', '/.../data/atlas-desikankilliany-rh.label.gii.gz')

.. NOTE: Add this in v1.0 and remove this in v1.3
.. .. warning::

    As of `v1.0` the ordering of the parcels in the Desikan-Killiany atlas
    provided with ``abagen`` has changed! This means that any results using
    this atlas will be different than those previously obtained with identical
    code. If you need the older version you can navigate to GitHub and obtain
    it from the history of the ``abagen`` repository. We strongly encourage you
    to use the new ordering, however, which is more consistent with traditional
    (i.e., FreeSurfer) conventions.

.. _usage_parcellations_additional:

Providing additional parcellation info
--------------------------------------

While only the image (i.e., NIFTI or GIFTIs) is required for processing the
microarray data, the CSV file with information on the parcellation scheme can
also be very useful. In particular, ``abagen`` can use the CSV file to
constrain the matching of tissue samples to anatomical regions in the atlas
image.

.. note::

    If you are using a surface atlas and your GIFTI files have valid `label
    tables <https://www.nitrc.org/projects/gifti/>`_ then ``abagen`` will
    automatically create a pandas.DataFrame with all the relevant information
    described below. However, you can always provide a separate CSV file if
    you are unsure and this will override any label tables present in the
    GIFTI files.

If you want to supply your own CSV file with information about an atlas you
must ensure it has (at least) the following columns:

  1. ``id``: an integer ID corresponding to the labels in the ``atlas`` image
  2. ``hemisphere``: a left/right/bilateral hemispheric designation (i.e., 'L',
     'R', or 'B')
  3. ``structure``: a broad structural class designation (i.e., one of
     'cortex', 'subcortex/brainstem', 'cerebellum', 'white matter', or
     'other')

For example, a valid CSV might look like this:

.. doctest::

    >>> import pandas as pd
    >>> atlas_info = pd.read_csv(atlas['info'])
    >>> print(atlas_info)
        id                    label hemisphere            structure
    0    1                 bankssts          L               cortex
    1    2  caudalanteriorcingulate          L               cortex
    2    3      caudalmiddlefrontal          L               cortex
    ..  ..                      ...        ...                  ...
    80  81              hippocampus          R  subcortex/brainstem
    81  82                 amygdala          R  subcortex/brainstem
    82  83                brainstem          B  subcortex/brainstem
    <BLANKLINE>
    [83 rows x 4 columns]

Notice that extra columns (i.e., ``label``) are okay as long as the three
required columns are present! If you want to confirm your file is formatted
correctly you can use :func:`abagen.images.check_atlas`:

.. note::

    Do not run :func:`abagen.images.check_atlas` below if you are following
    the tutorial for your first trial of abagen. It will convert the atlas
    from ``atlas['image'], atlas['info']`` into :obj:`abagen.AtlasTree` object.
    If you ran it accidentally, use ``atlas`` directly in the next steps,
    instead of ``atlas['image'], atlas['info']``. For more, see the section
    on :ref:`usage_parcellations_special`.

.. doctest::

    >>> from abagen import images
    >>> atlas = abagen.fetch_desikan_killiany()
    >>> # atlas = images.check_atlas(atlas['image'], atlas['info']);

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

.. _usage_parcellation_individual:

Individualized parcellations
----------------------------

Instead of providing a single parcellation image that will be used for all
donors, you can instead provide a parcellation image for each donor in the
space of their "raw" (or native) T1w image. ``abagen`` ships with versions
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

We also provide donor-specific surface atlases (derived from the FreeSurfer
outputs that can be fetched with :func:`abagen.datasets.fetch_freesurfer`).
These atlases are also shipped with ``abagen`` and can be loaded with:

.. doctest::

    >>> atlas = abagen.fetch_desikan_killiany(native=True, surface=True)
    >>> print(atlas['image'].keys())
    dict_keys(['9861', '10021', '12876', '14380', '15496', '15697'])
    >>> print(atlas['image']['9861'])  # doctest: +ELLIPSIS
    ('/.../9861/atlas-desikankilliany-lh.label.gii.gz', '/.../9861/atlas-desikankilliany-rh.label.gii.gz')

Note that if you are using your own donor-specific surface atlases they must,
by default, be based on the geometry of the FreeSurfer surfaces provided with
:func:`abagen.datasets.fetch_freesurfer`. If you wish to use surface atlases
based on different geometry please refer to :ref:`usage_parcellations_special`,
below.

Finally, when in doubt we recommend simply using a standard-space, group-level
atlas; however, we are actively investigating whether native-space atlases
provide any measurable benefits to the ``abagen`` workflows.

.. note::

    The donor-native volumetric versions of the DK parcellation shipped with
    ``abagen`` were generated by Arnatkevičiūte et al., 2018, *NeuroImage*, and
    are provided under the CC BY 4.0 license. The donor-native surface versions
    of the DK parcellation were generated by Romero-Garcia et al., 2017,
    *NeuroImage*, and are also provided under the CC BY 4.0 license.

.. _usage_parcellations_special:

Non-standard parcellations
--------------------------

If you'd like to use a non-standard atlas in the primary ``abagen`` workflows
that may be possible---with some caveats. That is, the constraining factor here
is the coordinates of the tissue samples from the AHBA: they are available in
(1) the native space of each donor's MRI, or (2) MNI152 space, and we strongly
encourage you to use one of these options (rather than e.g., attempting to
register the coordinates to a new space). If you provide a group-level atlas
the toolbox will default to using the MNI152 coordinates; if you provide
donor-specific atlases then the tooblox will use the native coordinates. Thus,
by default, ``abagen`` prefers you use one of the atlas conformations described
above.

However, if you have an atlas in a different space or resolution you can
(potentially) use it in the primary ``abagen`` workflows. To do this you will
need to create a :obj:`abagen.AtlasTree` object. All atlases provided are
internally coerced to `AtlasTree` instances, which is then used to assign
microarray tissue samples to parcels in the atlas.

Take, for example, a surface atlas in fsaverage6 resolution (by default,
surface atlases are assumed to be fsaverage5 resolution). In this case, you
simply need to supply the relevant geometry files for the atlas and specify
the space of the atlas:

.. code::

    >>> from abagen import images
    >>> atlas = ('/.../fsaverage6-lh.label.gii', '/.../fsaverage6-rh.label.gii')
    >>> surf = ('/.../fsaverage6-lh.surf.gii', '/.../fsaverage6-lh.surf.gii')
    >>> atlas = images.check_atlas(atlas, geometry=surf, space='fsaverage6')

The same procedure can be used for an atlas using fsLR geometry:

.. code::

    >>> from abagen import images
    >>> atlas = ('/.../fslr32k-lh.label.gii', '/.../fslr32k-rh.label.gii')
    >>> surf = ('/.../fslr32k-lh.surf.gii', '/.../fslr32k-lh.surf.gii')
    >>> atlas = images.check_atlas(atlas, geometry=surf, space='fslr')

.. _MNI space: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1088516/
.. _fsaverage space: https://surfer.nmr.mgh.harvard.edu/
