# -*- coding: utf-8 -*-
"""
Tests for abagen.process module
"""

import pandas as pd

from abagen import process


def test_load_alleninf_coords():
    coords = process._load_alleninf_coords()
    assert isinstance(coords, pd.DataFrame)
    assert coords.index.name == 'well_id'
    assert list(coords.columns) == ['mni_x', 'mni_y', 'mni_z']
    assert coords.shape == (3702, 3)
