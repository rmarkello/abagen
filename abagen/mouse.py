# -*- coding: utf-8 -*-

import requests
from xml.etree import ElementTree as ET
import numpy as np
from scipy.stats import zscore
import pandas as pd

API_QUERY_STRING = 'http://api.brain-map.org/api/v2/data/SectionDataSet/query.xml?id='
API_OPTION_STRING = '&include=structure_unionizes%28structure%29'


def get_mouse_expression_from_single_experiment(experiment_id, structure_list=None):
    """
    fetches mouse gene expression data of a single experiment,
    either saggital or coronal, according to the experiment id

    Parameters
        ----------
    experiment_id: an integer, specifying the experiment id
    structure_list: a list of strings, optional
        the list of ROIs in the form of acronyms.
        default: ROIs as included in Rubinov et al, 2015

    Returns
    -------
    gene_epr: array_like
        a (N, ) numpy array of regional gene expressions with N=len(structure_list),
        corresponding to the ROIs in structure_list

    Raises
    ------
    ValueError
        If `atlas_info` does not have sufficient information
    """
    if structure_list is None:
        # read default ROI list
        # (see Rubinov et al, 2015 for the criteria of choosing the ROIs)
        structure_list = pd.read_csv()

    query_url = API_QUERY_STRING + str(experiment_id) + API_OPTION_STRING

