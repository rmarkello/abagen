# -*- coding: utf-8 -*-
"""
Tests for abagen.samples module
"""

import pandas as pd

from abagen import samples


def test_replace_mni_coords(testfiles):
    for fn in testfiles['annotation']:
        out = samples._replace_mni_coords(fn)
        assert isinstance(out, pd.DataFrame)
        assert all(f in out.columns for f in ['mni_x', 'mni_y', 'mni_z'])
