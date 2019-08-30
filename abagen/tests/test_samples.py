# -*- coding: utf-8 -*-
"""
Tests for abagen.samples module
"""

import numpy as np
import pandas as pd
import pytest

from abagen import samples


def test_update_mni_coords(testfiles):
    # mni values are getting replaced so who cares about original
    # but ids are important and need to be real!
    x = y = z = [-10, 20, 30, 40]
    ids = [594, 2985, 1058, 1145]
    annotation = pd.DataFrame(dict(mni_x=x, mni_y=y, mni_z=z, well_id=ids))
    out = samples.update_mni_coords(annotation)

    # confirm that no samples were lost / reordered during the update process
    # and that we still have all our columns
    assert np.all(out['well_id'] == annotation['well_id'])
    assert np.all(out.columns == annotation.columns)

    # but DO confirm that _something_ changes about the dataframes (i.e., our
    # bogus coordinates should be different)
    with pytest.raises(AssertionError):
        pd.testing.assert_frame_equal(out, annotation)
    assert np.all(out['mni_x'] != annotation['mni_x'])
    assert np.all(out['mni_y'] != annotation['mni_y'])
    assert np.all(out['mni_z'] != annotation['mni_z'])

    # if we provide invalid well_ids we should receive an error!
    annotation['well_id'] = [594, 2985, 1058, 99999999999]
    with pytest.raises(KeyError):
        samples.update_mni_coords(annotation)

    # make sure we don't get any errors on real data, though
    for annotation in testfiles['annotation']:
        samples.update_mni_coords(annotation)


@pytest.mark.parametrize('path, expected', [
    ('/4005/4006/4007/4275/4276/4277/4278/12899/4286/', 'subcortex'),
    ('/4005/4006/4007/4275/4327/4341/4342/4344/', 'subcortex'),
    ('/4005/4006/4007/4008/4084/4103/4111/4112/4113/', 'cortex'),
    ('/4005/4006/4833/4696/4697/12930/12931/12933/4751/', 'cerebellum'),
    ('/4005/4006/9512/9676/9677/9680/9681/', 'brainstem'),
    ('/4005/4006/4833/9131/9132/9133/9492/', 'brainstem'),
    ('/4005/9218/9298/12959/265505622/', 'white matter'),
    ('/4005/9218/9219/9220/9227/', 'white matter'),
    ('/4005/9352/9418/9419/9708/', 'other'),
    ('/4005/9352/9353/9400/9402/', 'other'),
    ('/4005/059123', None),
    ('thisisnotapath', None),  # TODO: should this error?
])
def test__get_struct(path, expected):
    out = samples._get_struct(path)
    assert out == expected if expected is not None else out is None


def test_drop_mismatch_samples(testfiles):
    sids = ['4251', '4323', '059123']
    hemi = ['L', 'R', np.nan]
    paths = [
        '/4005/4006/4007/4008/4219/4249/12896/4251/',
        '/4005/4006/4007/4275/4321/4323/',
        '/4005/059123/',
    ]
    ontology = pd.DataFrame(dict(id=sids, hemisphere=hemi,
                                 structure_id_path=paths))

    # create annotation file with two mismatched samples (idxs 1 & 3)
    mni_x = [-10, -20, 30, 40, 0]
    sids = ['4251', '4323', '4323', '4251', '059123']
    annotation = pd.DataFrame(dict(mni_x=mni_x, structure_id=sids))

    # here's what we expect
    expected = pd.DataFrame(dict(mni_x=[-10, 30, 0],
                                 structure_id=['4251', '4323', '059123'],
                                 hemisphere=['L', 'R', np.nan],
                                 structure=['cortex', 'subcortex', None]),
                            index=[0, 2, 4])

    # do we get what we expect?
    out = samples.drop_mismatch_samples(annotation, ontology)
    pd.testing.assert_frame_equal(out, expected)

    # make sure we don't get errors on real data, either
    for an, on in zip(testfiles['annotation'], testfiles['ontology']):
        samples.drop_mismatch_samples(an, on)


@pytest.mark.xfail
def test__assign_sample(testfiles, atlas):
    assert False


@pytest.mark.xfail
def test__check_label(testfiles, atlas):
    assert False


@pytest.mark.xfail
def test_label_samples(testfiles, atlas):
    assert False


@pytest.mark.xfail
def test_mirror_samples(testfiles):
    assert False


@pytest.mark.xfail
def test__mirror_samples(testfiles):
    assert False


@pytest.mark.xfail
def test__mirror_ontology(testfiles):
    assert False
