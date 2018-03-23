# -*- coding: utf-8 -*-
"""
Functions for loading the various files associated with the Allen Brain Atlas
human microarray datasets.

Primarily here because I always forget which ones have headers and which ones
don't, so I figure this is a bit more foolproof.
"""

import pandas as pd


def read_microarray(fname):
    if not isinstance(fname, str):
        return fname
    return pd.read_csv(fname, header=None, index_col=0)


def read_ontology(fname):
    if not isinstance(fname, str):
        return fname
    return pd.read_csv(fname)


def read_pacall(fname):
    if not isinstance(fname, str):
        return fname
    return pd.read_csv(fname, header=None, index_col=0)


def read_probes(fname):
    if not isinstance(fname, str):
        return fname
    return pd.read_csv(fname, index_col=0)


def read_sampleannot(fname):
    if not isinstance(fname, str):
        return fname
    return pd.read_csv(fname)
