.. _api_ref:

.. currentmodule:: abagen

Reference API
=============

.. contents:: **List of modules**
   :local:

.. _ref_workflows:

:mod:`abagen.allen` - Primary workflows
---------------------------------------
.. automodule:: abagen.allen
   :no-members:
   :no-inherited-members:

.. currentmodule:: abagen

.. autosummary::
   :template: function.rst
   :toctree: generated/

   abagen.get_expression_data
   abagen.get_samples_in_mask
   abagen.get_interpolated_map

.. _ref_datasets:

:mod:`abagen.datasets` - Fetching AHBA datasets
-----------------------------------------------
.. automodule:: abagen.datasets
   :no-members:
   :no-inherited-members:

.. currentmodule:: abagen

.. autosummary::
   :template: function.rst
   :toctree: generated/

   abagen.fetch_desikan_killiany
   abagen.fetch_donor_info
   abagen.fetch_freesurfer
   abagen.fetch_gene_group
   abagen.fetch_microarray
   abagen.fetch_raw_mri
   abagen.fetch_rnaseq
   abagen.datasets.fetch_fsaverage5
   abagen.datasets.fetch_fsnative

.. _ref_images:

:mod:`abagen.images` - Image processing functions
-------------------------------------------------
.. automodule:: abagen.images
   :no-members:
   :no-inherited-members:

.. currentmodule:: abagen

.. autosummary::
   :template: function.rst
   :toctree: generated/

   abagen.leftify_atlas
   abagen.check_atlas
   abagen.annot_to_gifti
   abagen.relabel_gifti

.. _ref_correct:

:mod:`abagen.correct` - Post-processing corrections
---------------------------------------------------
.. automodule:: abagen.correct
   :no-members:
   :no-inherited-members:

.. currentmodule:: abagen

.. autosummary::
   :template: function.rst
   :toctree: generated/

   abagen.remove_distance
   abagen.keep_stable_genes
   abagen.normalize_expression

.. _ref_matching:

:mod:`abagen.matching` - Functions for matching samples
-------------------------------------------------------
.. automodule:: abagen.matching
   :no-members:
   :no-inherited-members:

.. currentmodule:: abagen

.. autosummary::
   :template: class.rst
   :toctree: generated/

   abagen.AtlasTree

.. currentmodule:: abagen

.. autosummary::
   :template: function.rst
   :toctree: generated/

   abagen.matching.get_centroids
   abagen.matching.closest_centroid

.. _ref_reporting:

:mod:`abagen.reporting` - Functions for generating reports
----------------------------------------------------------
.. automodule:: abagen.reporting
   :no-members:
   :no-inherited-members:

.. currentmodule:: abagen

.. autosummary::
   :template: class.rst
   :toctree: generated/

   abagen.Report

.. _ref_io:

:mod:`abagen.io` - Loading AHBA data files
------------------------------------------
.. automodule:: abagen.io
   :no-members:
   :no-inherited-members:

.. currentmodule:: abagen.io

.. autosummary::
   :template: function.rst
   :toctree: generated/

   abagen.io.read_annotation
   abagen.io.read_microarray
   abagen.io.read_ontology
   abagen.io.read_pacall
   abagen.io.read_probes
   abagen.io.read_genes
   abagen.io.read_tpm
   abagen.io.read_counts

.. _ref_mouse:

:mod:`abagen.mouse` - Working with the Allen Mouse Brain Atlas
--------------------------------------------------------------
.. automodule:: abagen.mouse
   :no-members:
   :no-inherited-members:

.. currentmodule:: abagen.mouse

.. autosummary::
   :template: function.rst
   :toctree: generated/

   abagen.mouse.available_gene_info
   abagen.mouse.available_structure_info
   abagen.mouse.available_unionization_info

   abagen.mouse.get_gene_info
   abagen.mouse.get_structure_info
   abagen.mouse.get_structure_coordinates
   abagen.mouse.get_unionization_from_gene

   abagen.mouse.fetch_allenref_genes
   abagen.mouse.fetch_allenref_structures
   abagen.mouse.fetch_rubinov2015_structures
