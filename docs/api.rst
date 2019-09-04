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

   abagen.fetch_microarray
   abagen.fetch_raw_mri
   abagen.fetch_desikan_killiany
   abagen.fetch_gene_group

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

:mod:`abagen.utils` - Utility functions
---------------------------------------
.. automodule:: abagen.utils
   :no-members:
   :no-inherited-members:

.. currentmodule:: abagen.utils

.. autosummary::
   :template: function.rst
   :toctree: generated/

   abagen.utils.check_atlas_info
   abagen.utils.check_metric
