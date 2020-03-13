---
title: 'abagen: A Python toolbox for the Allen Brain Atlas'
tags:
  - Python
  - neuroimaging
  - transcriptomics
authors:
  - name: Ross D. Markello
    orcid: 0000-0003-1057-1336
    affiliation: 1
  - name: Golia Shafiei
    orcid: 0000-0002-2036-5571
    affiliation: 1
  - name: Ying-Qiu Zheng
    orcid: 0000-0003-1236-0700
    affiliation: 2
  - name: Vincent Bazinet
    affiliation: 1
  - name: Anthony Gifuni
    affilation: 1
  - name: Ellen Wang
    affiliation: 1
  - name: Romke Hannema
    affiliation: 1
  - name: Bratislav Misic
    orcid: 0000-0003-0307-2862
    affiliation: 1
affiliations:
  - name: McGill University, Montreal, Canada
    index: 1
  - name: University of Oxford, Oxford, United Kingdom
    index: 2
date: 17 February 2020
bibliography: paper.bib
---

# Summary

In 2010, the Allen Institute for Brain Science released the [Allen Human Brain
Atlas](https://human.brain-map.org/) (AHBA), a dataset containing microarray
expression data collected from post-mortem human brains [@hawrylycz:2012]. This
dataset has offered an unprecedented opportunity to examine the genetic
scaffolding of the human brain, and has already yielded novel insights into
e.g., adolescent brain development [@whitakervertes:2016] and the functional
organization of the brain [@richiardi:2015].

In order to be effectively used in most research analyses, data from the Allen
Human Brain Atlas needs to be (1) pre-processed and normalized, (2) aggregated
into distinct brain regions, and (3) combined across donors. Notably, there are
a number of analytic choices in these steps that can dramatically influence the
resulting data and any downstream analyses. @arnatkeviciute:2019 provided a
thorough overview of this, demonstrating how the techniques and code used to
prepare the raw data from the AHBA have varied widely across the published
literature.

``abagen`` provides reproducible workflows for downloading, preparing, and
processing data from the AHBA, accepting any brain atlas—a file that specifies
the division of the brain into distinct regions—and returning ready-to-analyze
microarray expression data.

For a more detailed explanation about how abagen works and instructions on how
to install the software refer to our ReadTheDocs documentation
([https://abagen.readthedocs.io](https://abagen.readthedocs.io)).

# Acknowledgement

This research was undertaken thanks in part tofunding from the Canada First
Research Excellence Fund, awarded to McGill University for the Healthy Brains
for Healthy Lives initiative. BM acknowledges support from the Natural Sciences
and Engineering Research Council of Canada (NSERC Discovery Grant RGPIN
\#017-04265), from the Canada Research Chairs Program and from the Fonds de
Recherche du Quebec - Sante (Chercheur Boursier).

# References
