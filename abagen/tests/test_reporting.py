# -*- coding: utf-8 -*-
"""
Tests for abagen.reporting module
"""

from abagen import reporting


def test_sanitize_text():
    assert reporting._sanitize_text('normal text') == 'normal text'
    assert reporting._sanitize_text('header\n\n<br> this is\n a   test ') == \
        'header \n\nthis is a test'


def test_get_donor_demographics():
    out = reporting._get_donor_demographics('all')
    keys = ('n_donors', 'n_female', 'min', 'max', 'mean', 'std')
    assert all(key in out for key in keys)
    assert out['n_donors'] == 6 and out['n_female'] == 1

    out = reporting._get_donor_demographics(['9861', '10021'])
    assert all(key in out for key in keys)
    assert out['n_donors'] == 2 and out['n_female'] == 0


def test_get_norm_procedure():
    pass


def test_add_references():
    assert reporting._add_references('') == ''
    assert reporting._add_references('no references') == ''
    assert reporting._add_references('[MISSNO]') == ''
    assert reporting._add_references('[A2019N]') != ''


def test_reports(atlas, surface):
    # not sure how else to check this beyond a simple smoke test
    # i'm certainly not gonna recode all the logic here...
    volreport = reporting.Report(atlas['image'], atlas['info'])
    assert isinstance(volreport, reporting.Report)
    assert volreport.body != ''
    assert "83-region volumetric atlas in MNI space" in volreport.body

    # check that surface returns different output than volumetric atlas
    surfreport = reporting.Report(surface['image'], surface['info'])
    assert isinstance(surfreport, reporting.Report)
    assert surfreport.body != ''
    assert "68-region surface-based atlas in MNI space" in surfreport.body
    assert surfreport.body != volreport.body

    # check a couple of specific parameter choices :man_shrugging:
    report = reporting.Report(atlas['image'], atlas['info'],
                              lr_mirror='bidirectional', missing='centroids',
                              sample_norm='srs', gene_norm='zscore').body
    assert 'tissue samples were mirrored' in report
    assert 'tissue sample closest to the centroid of' in report
    assert 'variation in gene expression was then addressed' in report
