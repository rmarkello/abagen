abagen: A toolbox for the Allen Brain Atlas genetics data
=========================================================

This package provides a Python interface for working with the `Allen Human
Brain Atlas <http://human.brain-map.org/>`_ (AHBA) microarray expression data.

.. image:: https://travis-ci.org/rmarkello/abagen.svg?branch=master
   :target: https://travis-ci.org/rmarkello/abagen
.. image:: https://codecov.io/gh/rmarkello/abagen/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/rmarkello/abagen
.. image:: https://readthedocs.org/projects/abagen/badge/?version=latest
   :target: http://abagen.readthedocs.io/en/latest
.. image:: https://img.shields.io/badge/License-BSD%203--Clause-blue.svg
   :target: https://opensource.org/licenses/BSD-3-Clause

.. _overview:

Overview
--------

In 2013, the Allen Institute for Brain Science released the Allen Human Brain
Atlas, `a dataset <http://human.brain-map.org/>`_ containing microarray
expression data collected from six human brains. This dataset has offered an
unprecedented opportunity to examine the genetic underpinnings of the human
brain, and has already yielded novel insight into e.g., `adolescent brain
development <http://www.pnas.org/content/113/32/9105.long>`_ and `functional
brain organization <http://science.sciencemag.org/content/348/6240/1241.long>`_.

However, in order to be effectively used in most analyses, the AHBA microarray
expression data often needs to be (1) collapsed into regions of interest (e.g.,
parcels or networks), and (2) combined across donors. While this may
potentially seem trivial, there are numerous analytic choices in these steps
that can dramatically influence the resulting data and any downstream analyses.
Indeed, Arnatkevičiūte et al., 2018 ([1]_) provided a thorough treatment of
this in a `recent manuscript <https://www.biorxiv.org/content/early/2018/07/30/
380089>`_, demonstrating how the techniques and code used to prepare the raw
AHBA data have varied widely across published reports.

The current Python package, ``abagen``, aims to provide a reproducible pipeline
for processing and preparing the AHBA microarray expression data for analysis.
If you'd like more information about the package, including how to install it
and some example instructions on its use, check out our `documentation <https:
//abagen.readthedocs.io>`_!

.. _development:

Development and getting involved
--------------------------------

This package has been largely developed in the spare time of a single graduate
student (`@rmarkello <https://github.com/rmarkello>`_) with help from some
incredible `contributors <https://github.com/rmarkello/abagen/graphs/
contributors>`_. While it would be |sparkles| amazing |sparkles| if anyone else
finds it helpful, given the limited time constraints of graduate school, the
current package is not currently accepting requests for new features.

However, if you're interested in getting involved in the project, we're
thrilled to welcome new contributors! You shouldstart by reading our
`contributing guidelines <https://github.com/rmarkello/abagen/blob/master/
CONTRIBUTING.md>`_ and `code of conduct <https://github.com/rmarkello/abagen/
blob/master/CODE_OF_CONDUCT.md>`_. Once you're done with that, take a look at
our `issues <https://github.com/rmarkello/abagen/issues>`_ to see if there's
anything you might like to work on. Alternatively, if you've found a bug, are
experiencing a problem, or have a question, create a new issue with some
information about it!

.. _acknowledgments:

Acknowledgments
---------------

While this package was initially created in early 2018, many of the current
functions in the project were inspired by the workflow laid out in
Arnatkevičiūte et al., 2018. As such, if you use this code it would be good
to (1) provide a link back to the ``abagen`` repository with the version of the
code used, and (2) cite their paper:

.. [1] Arnatkeviciute, A., Fulcher, B. D., & Fornito, A. (2018). A practical
   guide to linking brain-wide gene expression and neuroimaging data. bioRxiv,
   380089.

Additionally, please add the following publication to all Acknowledgments and References sections per The Allen Institute's `citation policies <https://alleninstitute.org/legal/citation-policy/>`_:

.. [2] Hawrylycz, M.J. et al. (2012) An anatomically comprehensive atlas of the adult human transcriptome, Nature 489: 391-399. doi:10.1038/nature11405

.. _licensing:

License Information
-------------------

This codebase is licensed under the 3-clause BSD license. The full license can
be found in the `LICENSE <https://github.com/rmarkello/abagen/blob/master/
LICENSE>`_ file in the ``abagen`` distribution.

Reannotated gene information located at ``abagen/data/reannotated.csv.gz`` is
taken from [1]_ and is separately licensed under the `CC BY 4.0 <https://
creativecommons.org/licenses/by/4.0/legalcode>`_; these data can also be found
on `figshare <https://figshare.com/s/441295fe494375aa0c13>`_.

Corrected MNI coordinates used to match AHBA tissues samples to MNI space
located at ``abagen/data/corrected_mni_coordinates.csv`` are taken from the
`alleninf package <https://github.com/chrisfilo/alleninf>`_, provided under
the `3-clause BSD license <https://opensource.org/licenses/BSD-3-Clause>`_.

All microarray expression data is copyrighted under `non-commercial reuse
policies <http://alleninstitute.org/legal/terms-use/>`_ by the Allen Institute
for Brain Science (© 2010 Allen Institute for Brain Science. Allen Human Brain
Atlas. Available from: human.brain-map.org).

All trademarks referenced herein are property of their respective holders.

.. |sparkles| replace:: ✨
