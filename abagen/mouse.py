# -*- coding: utf-8 -*-

import requests
from xml.etree import ElementTree as ET
import numpy as np
from scipy.stats import zscore
import pandas as pd

API_QUERY_STRING = 'http://api.brain-map.org/api/v2/data/SectionDataSet/query.xml?id='
API_OPTION_STRING = '&include=structure_unionizes%28structure%29'


def get_mouse_expression_from_single_experiment(experiment_id, roi_list=None):
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
    if roi_list is None:
        # read default ROI list
        # (see Rubinov et al, 2015 for the criteria of choosing the ROIs)
        roilabels = pd.read_csv("abagen/data/roilabels-rubinov2015pnas.csv")
        roi_list = roilabels['roiacronyms']

    all_structure = []
    all_epr = np.empty((0, 1))

    # make the query
    query_url = API_QUERY_STRING + str(experiment_id) + API_OPTION_STRING
    r = requests.get(query_url)
    root = ET.fromstring(r.content)

    for item in root.findall(
        'section-data-sets/section-data-set/structure-unionizes/structure-unionize'
    ):
        # append new epr value
        for subitem in item.findall('expression-energy'):
            all_epr = np.append(all_epr, float(subitem.text))
        # append new structure label
        for subitem in item.findall('structure/acronym'):
            all_structure.append(subitem.text)

    # extract the regions in roi_list
    data_in_the_list = [
        (item, all_epr[k]) for k, item in enumerate(all_structure) if item in roi_list
    ]
    structure_in_the_list = [item[0] for item in data_in_the_list]
    epr_in_the_list = np.array(
        [item[1] for item in data_in_the_list]
    )

    

