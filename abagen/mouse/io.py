# -*- coding: utf-8 -*-
"""
functions to read default mouse gene list and structure list
"""
import pandas as pd


def read_all_genes():
    """
    read all available gene acronyms
    that have section data sets available
    in mouse atlasa in the form of acronyms
    :return: list
        a list of gene acronyms
    """
    all_gene_acronyms = pd.read_csv(
        "abagen/data/all_genes_available.csv"
    )
    return all_gene_acronyms['acronym'].values


def read_all_structures():
    """
    read default structure (ROI) acronyms as documented in
    Rubinov et al, 2015
    :return: list
        a list of structure acronyms
    """
    roilabels = pd.read_csv(
        "abagen/data/roilabels-rubinov2015pnas.csv"
    )
    return roilabels['roiacronyms'].values
