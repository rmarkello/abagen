# -*- coding: utf-8 -*-
"""
Defines information about the abagen distribution
"""

long_description = """
abagen: A toolbox for the Allen Brain Atlas genetics data
=========================================================

This package provides a Python interface for fetching and working with the
`Allen Human Brain Atlas`_ (AHBA) microarray expression data.

.. image:: https://circleci.com/gh/rmarkello/abagen.svg?style=shield
   :target: https://circleci.com/gh/rmarkello/abagen
.. image:: https://codecov.io/gh/rmarkello/abagen/branch/main/graph/badge.svg
   :target: https://codecov.io/gh/rmarkello/abagen
.. image:: https://readthedocs.org/projects/abagen/badge/?version=latest
   :target: https://abagen.readthedocs.io/en/stable
.. image:: https://img.shields.io/badge/License-BSD%203--Clause-blue.svg
   :target: https://opensource.org/licenses/BSD-3-Clause
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3451463.svg
   :target: https://doi.org/10.5281/zenodo.3451463

.. _readme_overview:

Overview
--------

In 2013, the Allen Institute for Brain Science released the `Allen Human Brain
Atlas`_, a dataset containing microarray expression data collected from six
human brains (Hawrylycz et al., 2012) . This dataset has offered an
unprecedented opportunity to examine the genetic underpinnings of the human
brain, and has already yielded novel insight into e.g., `adolescent brain
development <https://www.pnas.org/content/113/32/9105.long>`__ and `functional
brain organization <https://science.sciencemag.org/content/348/6240/
1241.long>`__.

However, in order to be effectively used in most analyses, the AHBA microarray
expression data often needs to be (1) collapsed into regions of interest (e.g.,
parcels or networks), and (2) combined across donors. While this may
potentially seem trivial, there are a number of analytic choices in these steps
that can dramatically influence the resulting data and any downstream analyses.
Arnatkevičiūte et al., 2019 provided a thorough treatment of this in a `recent
manuscript <https://www.sciencedirect.com/science/article/pii/
S1053811919300114>`__, demonstrating how the techniques and code used to
prepare the raw AHBA data have varied widely across published reports.

The current Python package, ``abagen``, aims to provide reproducible workflows
for processing and preparing the AHBA microarray expression data for analysis.

.. _readme_requirements:

Installation requirements
-------------------------

Currently, ``abagen`` works with Python 3.6+ and requires a few dependencies:

    - nibabel
    - numpy (>=1.14.0)
    - pandas (>=0.25.0), and
    - scipy

There are some additional (optional) dependencies you can install to speed up
some functions:

    - fastparquet, and
    - python-snappy

These latter packages are primarily used to facilitate loading the (rather
large!) microarray expression dataframes provided by the Allen Institute,

For detailed information on how to install ``abagen``, including these
dependencies, refer to our `installation instructions`_.

.. _readme_usage:

Quickstart
----------

At it's core, using ``abagen`` is as simple as:

.. code-block:: python

    >>> import abagen
    >>> expression = abagen.get_expression_data('myatlas.nii.gz')  # doctest: +SKIP

where ``'myatlas.nii.gz'`` points to a brain parcellation file.

This function can also be called from the command line with:

.. code-block:: bash

    $ abagen --output-file expression.csv myatlas.nii.gz

For more detailed instructions on how to use ``abagen`` please refer to our
`user guide`_!

.. _readme_development:

Development and getting involved
--------------------------------

If you've found a bug, are experiencing a problem, or have a question about
using the package, please head on over to our `GitHub issues`_ and make a new
issue with some information about it! Someone will try and get back to you
as quickly as possible, though please note that the primary developer for
``abagen`` (@rmarkello) is a graduate student so responses make take some time!

If you're interested in getting involved in the project: welcome |sparkles|!
We're thrilled to welcome new contributors. You should start by reading our
`code of conduct`_; all activity on ``abagen`` should adhere to the CoC. After
that, take a look at our `contributing guidelines`_ so you're familiar with the
processes we (generally) try to follow when making changes to the repository!
Once you're ready to jump in head on over to our issues to see if there's
anything you might like to work on.

.. _readme_acknowledgments:

Citing ``abagen``
-----------------

For up-to-date instructions on how to cite abagen please refer to our
`documentation <https://abagen.readthedocs.io/en/stable/citing.html>`_.

.. _readme_licensing:

License Information
-------------------

This codebase is licensed under the `3-clause BSD license`_. The full license
can be found in the `LICENSE`_ file in the ``abagen`` distribution.

Reannotated gene information located at ``abagen/data/reannotated.csv.gz`` and
individualized donor parcellations for the Desikan-Killiany atlas located at
``abagen/data/native_dk`` are taken from Arnatkevičiūte et al., 2018 and are
separately licensed under the `CC BY 4.0`_; these data can also be found on
`figshare <https://figshare.com/s/441295fe494375aa0c13>`__.

Corrected MNI coordinates used to match AHBA tissues samples to MNI space
located at ``abagen/data/corrected_mni_coordinates.csv`` are taken from the
`alleninf package`_, provided under the `3-clause BSD license`_.

All microarray expression data is copyrighted under `non-commercial reuse
policies`_ by the Allen Institute for Brain Science (© 2010 Allen Institute for
Brain Science. Allen Human Brain Atlas. Available from: `Allen Human Brain
Atlas`_).

All trademarks referenced herein are property of their respective holders.

.. |sparkles| replace:: ✨
.. |warning| replace:: 🚨
.. _3-clause BSD license: https://opensource.org/licenses/BSD-3-Clause
.. _Allen Human Brain Atlas: https://human.brain-map.org/
.. _alleninf package: https://github.com/chrisfilo/alleninf
.. _CC BY 4.0: https://creativecommons.org/licenses/by/4.0/legalcode
.. _code of conduct: https://github.com/rmarkello/abagen/blob/main/CODE_OF_CONDUCT.md
.. _contributing guidelines: https://github.com/rmarkello/abagen/blob/main/CONTRIBUTING.md
.. _contributors: https://github.com/rmarkello/abagen/graphs/contributors
.. _user guide: https://abagen.readthedocs.io/en/stable/usage.html
.. _GitHub issues: https://github.com/rmarkello/abagen/issues
.. _installation instructions: https://abagen.readthedocs.io/en/stable/installation.html
.. _LICENSE: https://github.com/rmarkello/abagen/blob/main/LICENSE
.. _non-commercial reuse policies: https://alleninstitute.org/legal/terms-use/
"""  # noqa
