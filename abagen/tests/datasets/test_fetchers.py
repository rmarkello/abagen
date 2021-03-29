# -*- coding: utf-8 -*-
"""
Tests for abagen.datasets.fetchers module
"""

from pathlib import Path

import pandas as pd
import pytest

from abagen.datasets import fetchers


def test_fetch_microarray():
    # test different `donors` inputs
    f1 = fetchers.fetch_microarray(donors=['12876'])
    f2 = fetchers.fetch_microarray(donors='12876')
    f3 = fetchers.fetch_microarray(donors='H0351.1009')
    f4 = fetchers.fetch_microarray(donors=None)

    # test n_proc
    fetchers.fetch_microarray(donors=['12876', '15496'], n_proc=2)

    assert f1 == f2 == f3 == f4
    assert len(f1) == 1  # only one donor
    for k in ['microarray', 'annotation', 'pacall', 'probes', 'ontology']:
        assert k in f1['12876']

    # check downloading incorrect donor
    with pytest.raises(ValueError):
        fetchers.fetch_microarray(donors='notadonor')

    with pytest.raises(ValueError):
        fetchers.fetch_microarray(donors=['notadonor'])


def test_fetch_raw_mri():
    # test different `donors` inputs
    f1 = fetchers.fetch_raw_mri(donors=['12876'])
    f2 = fetchers.fetch_raw_mri(donors='12876')
    f3 = fetchers.fetch_raw_mri(donors='H0351.1009')
    f4 = fetchers.fetch_raw_mri(donors=None)

    assert f1 == f2 == f3 == f4
    assert len(f1) == 1
    for k in ['t1w', 't2w']:
        assert k in f1['12876']

    # check downloading incorrect donor
    with pytest.raises(ValueError):
        fetchers.fetch_raw_mri(donors='notadonor')

    with pytest.raises(ValueError):
        fetchers.fetch_raw_mri(donors=['notadonor'])


def test_fetch_rnaseq():
    f1 = fetchers.fetch_rnaseq(donors=['9861'])
    f2 = fetchers.fetch_rnaseq(donors='9861')
    f3 = fetchers.fetch_rnaseq(donors='H0351.2001')
    f4 = fetchers.fetch_rnaseq(donors=None)

    assert f1 == f2 == f3 == f4
    assert len(f1) == 1
    for k in ['genes', 'ontology', 'counts', 'tpm', 'annotation']:
        assert k in f1['9861']

    with pytest.raises(ValueError):
        fetchers.fetch_rnaseq(donors='notadonor')

    with pytest.raises(ValueError):
        fetchers.fetch_rnaseq(donors=['9861', 'notadonor'])


@pytest.mark.parametrize('group, expected', [
    ('brain', 2413),
    ('neuron', 2530),
    ('oligodendrocyte', 1769),
    ('synaptome', 1886),
    ('layers', 46)
])
def test_get_gene_group(group, expected):
    assert len(fetchers.fetch_gene_group(group)) == expected

    with pytest.raises(ValueError):
        fetchers.fetch_gene_group('notagroup')


def test_fetch_desikan_killiany():
    atlas = fetchers.fetch_desikan_killiany()
    assert sorted(atlas) == ['image', 'info']
    img, info = Path(atlas['image']), Path(atlas['info'])
    assert img.name == 'atlas-desikankilliany.nii.gz' and img.exists()
    assert info.name == 'atlas-desikankilliany.csv' and info.exists()

    atlas = fetchers.fetch_desikan_killiany(native=True)
    donors = ['10021', '12876', '14380', '15496', '15697', '9861']
    assert sorted(atlas) == ['image', 'info']
    assert sorted(atlas['image']) == donors
    for donor in donors:
        img = Path(atlas['image'][donor])
        assert img.name == 'atlas-desikankilliany.nii.gz' and img.exists()
    info = Path(atlas['info'])
    assert info.name == 'atlas-desikankilliany.csv' and info.exists()

    atlas = fetchers.fetch_desikan_killiany(surface=True)
    assert sorted(atlas) == ['image', 'info']
    assert len(atlas['image']) == 2
    for hemi, img in zip(('lh', 'rh'), atlas['image']):
        img = Path(img)
        assert img.name == f'atlas-desikankilliany-{hemi}.label.gii.gz'
        assert img.exists()
    assert Path(atlas['info']).name == 'atlas-desikankilliany.csv'

    atlas = fetchers.fetch_desikan_killiany(native=True, surface=True)
    donors = ['10021', '12876', '14380', '15496', '15697', '9861']
    assert sorted(atlas) == ['image', 'info']
    assert sorted(atlas['image']) == donors
    for donor in donors:
        imgs = atlas['image'][donor]
        assert len(imgs) == 2
        for img in imgs:
            assert img.endswith('.gii.gz')
            assert Path(img).exists()
    info = Path(atlas['info'])
    assert info.name == 'atlas-desikankilliany.csv' and info.exists()


def test_fetch_donor_info():
    cols = [
        'donor', 'uid', 'age', 'sex', 'ethnicity', 'medical_conditions',
        'post_mortem_interval_hours'
    ]
    donor_info = fetchers.fetch_donor_info()
    assert isinstance(donor_info, pd.DataFrame)
    assert donor_info.shape == (6, 7)
    assert all(donor_info.columns == cols)


def test_fetch_freesurfer():
    # test different `donors` inputs
    f1 = fetchers.fetch_freesurfer(donors=['12876'])
    f2 = fetchers.fetch_freesurfer(donors='12876')
    f3 = fetchers.fetch_freesurfer(donors='H0351.1009')
    f4 = fetchers.fetch_freesurfer(donors=None)

    assert f1 == f2 == f3 == f4

    fpath = Path(f1['12876'])
    assert fpath.exists() and fpath.is_dir()

    for sd in ('label', 'mri', 'stats', 'surf'):
        sdpath = fpath / sd
        assert sdpath.exists() and sdpath.is_dir()

    # check downloading incorrect donor
    with pytest.raises(ValueError):
        fetchers.fetch_freesurfer(donors='notadonor')

    with pytest.raises(ValueError):
        fetchers.fetch_freesurfer(donors=['notadonor'])


def test_fetch_fsaverage5():
    fs5 = fetchers.fetch_fsaverage5()
    assert len(fs5) == 2
    for hemi in ('lh', 'rh'):
        assert hasattr(fs5, hemi)
        hemi = getattr(fs5, hemi)
        for attr, exp in zip(('vertices', 'faces'), (10242, 20480)):
            assert hasattr(hemi, attr)
            assert getattr(hemi, attr).shape == (exp, 3)

    fs5 = fetchers.fetch_fsaverage5(load=False)
    for hemi in ('lh', 'rh'):
        assert hasattr(fs5, hemi)
        hemi = getattr(fs5, hemi)
        assert Path(hemi).is_file()


def test_fetch_fsnative():
    fsn = fetchers.fetch_fsnative(donors=['12876'])
    fetchers.fetch_fsnative(donors='12876')
    fetchers.fetch_fsnative(donors='H0351.1009')
    fetchers.fetch_fsnative(donors=None)

    assert len(fsn) == 2
    for hemi in ('lh', 'rh'):
        assert hasattr(fsn, hemi)
        hemi = getattr(fsn, hemi)
        for attr in ('vertices', 'faces'):
            assert hasattr(hemi, attr)

    fsn = fetchers.fetch_fsnative(donors=['12876'], load=False)
    for hemi in ('lh', 'rh'):
        assert hasattr(fsn, hemi)
        hemi = getattr(fsn, hemi)
        assert Path(hemi).is_file()
