# -*- coding: utf-8 -*-
"""
functions to read default mouse gene list and structure list
"""
import pandas as pd


def read_all_genes(entry_type='acronym'):
    """
    read all available gene acronyms
    that have section data sets available
    in mouse atlas in the form of acronyms
    :param entry_type: str, optional
        the type of gene identifier. default:'acronym'
        supported:
            'id', return the list of gene IDs (integers)
            'acronym', return the list of gene acronyms (str)
            'name', return the list of gene names (str)
    :return: list
        a list of gene acronyms
    """
    all_gene_acronyms = pd.read_csv(
        "abagen/data/all_genes_available.csv"
    )
    return all_gene_acronyms[entry_type].values


def read_all_structures(entry_type='acronym'):
    """
    read default structure (ROI) acronyms as documented in
    Rubinov et al, 2015
    :param entry_type: str, optional
        the type of structure identifier. default:'acronym'
        supported:
            'id', return the list of structure IDs (integers)
            'acronym', return the list of structure acronyms (str)
            'name', return the list of structure names (str)

    :return: list
        a list of structure acronyms
    """
    roilabels = pd.read_csv(
        "abagen/data/roilabels-rubinov2015pnas.csv"
    )
    return roilabels[entry_type].values
