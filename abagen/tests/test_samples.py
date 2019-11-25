# -*- coding: utf-8 -*-
"""
Tests for abagen.samples module
"""

import numpy as np
import pandas as pd
import pytest

from abagen import samples


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  generate fake data (based largely on real data) so we know what to expect  #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


@pytest.fixture(scope='module')
def ontology():
    """ Fake ontology dataframe
    """
    sid = [4251, 4260, 4322, 4323, 9422]
    hemi = ['L', 'R', 'L', 'R', np.nan]
    acronym = ['S', 'S', 'Cl', 'Cl', 'CC']
    path = [
        '/4005/4006/4007/4008/4219/4249/12896/4251/',
        '/4005/4006/4007/4008/4219/4249/12896/4260/',
        '/4005/4006/4007/4275/4321/4322/',
        '/4005/4006/4007/4275/4321/4323/',
        '/4005/9352/9418/9422/',
    ]
    name = [
        'subiculum, left',
        'subiculum, right',
        'claustrum, left',
        'claustrum, right',
        'central canal',
    ]
    return pd.DataFrame(dict(id=sid, hemisphere=hemi, name=name,
                             acronym=acronym, structure_id_path=path))


@pytest.fixture(scope='module')
def mm_annotation():
    """ Fake annotation dataframe with some samples mislabelled
    """
    mni_x = [-10, -20, 30, 40, 0]
    sid = [4251, 4323, 4323, 4251, 9422]
    sacr = ['S', 'Cl', 'Cl', 'S', 'CC']
    sname = [
        'subiculum, left',
        'claustrum, right',
        'claustrum, right',
        'subiculum, left',
        'central canal'
    ]
    ind = pd.Series(range(len(sid)), name='sample_id')
    return pd.DataFrame(dict(mni_x=mni_x, structure_id=sid,
                             structure_acronym=sacr, structure_name=sname),
                        index=ind)


@pytest.fixture(scope='module')
def annotation(mm_annotation):
    """ Fake annotation dataframe
    """
    out = mm_annotation.loc[[0, 2, 4]].reset_index(drop=True)
    out.index.name = 'sample_id'
    return out


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# test all the functions on our generated fake data so we know what to expect #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


def test_update_mni_coords():
    # xyz coordinates are getting replaced so who cares about the original
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
    ('/4005/4006/4833', None),
    ('thisisnotapath', None),  # TODO: should this error?
])
def test__get_struct(path, expected):
    out = samples._get_struct(path)
    assert out == expected if expected is not None else out is None


def test_drop_mismatch_samples(mm_annotation, ontology):
    # here's what we expect (i.e., indices 1 & 3 are dropped and the structure
    # for the remaining samples is correctly extracted from the paths)
    expected = pd.DataFrame(dict(hemisphere=['L', 'R', np.nan],
                                 mni_x=[-10, 30, 0],
                                 structure_acronym=['S', 'Cl', 'CC'],
                                 structure=['subcortex', 'subcortex', 'other'],
                                 structure_id=[4251, 4323, 9422],
                                 structure_name=['subiculum, left',
                                                 'claustrum, right',
                                                 'central canal']),
                            index=[0, 2, 4])

    # do we get what we expect? (ignore ordering of columns / index)
    out = samples.drop_mismatch_samples(mm_annotation, ontology)
    pd.testing.assert_frame_equal(out, expected, check_like=True)


def test__mirror_ontology(annotation, ontology):
    aexp = pd.DataFrame(dict(mni_x=[-10, 30, 0],
                             structure_acronym=['S', 'Cl', 'CC'],
                             structure_id=[4260, 4322, 9422],
                             structure_name=['subiculum, right',
                                             'claustrum, left',
                                             'central canal']),
                        index=pd.Series(range(3), name='sample_id'))

    # this function doesn't touch mni_x -- it just assumes that all the
    # hemisphere designation are incorrect and updates the structure id, name,
    # and acronyms accordingly
    a = samples._mirror_ontology(annotation, ontology)
    pd.testing.assert_frame_equal(a, aexp, check_like=True)


def test_mirror_samples(annotation, ontology):
    # we're changing quite a bit of stuff in the annotation dataframe
    aexp = pd.DataFrame(dict(mni_x=[-10, 30, 0, 10, -30],
                             structure_acronym=['S', 'Cl', 'CC', 'S', 'Cl'],
                             structure_id=[4251, 4323, 9422, 4260, 4322],
                             structure_name=['subiculum, left',
                                             'claustrum, right',
                                             'central canal',
                                             'subiculum, right',
                                             'claustrum, left']),
                        index=[0, 1, 2, 0, 1])

    # but let's confirm all the outputs are as-expected
    a = samples.mirror_samples(annotation, ontology)
    pd.testing.assert_frame_equal(a, aexp, check_like=True)


@pytest.mark.xfail
def test__assign_sample(atlas):
    atlas = atlas['image']
    assert samples._assign_sample([[0, 0, 0]], atlas, tolerance=0) == 0
    assert samples._assign_sample([[26, 96, 66]], atlas, tolerance=0) == 71
    assert False


@pytest.mark.xfail
def test__check_label():
    assert False


@pytest.mark.xfail
def test_label_samples():
    assert False


def test_groupby_index():
    # default usage (no params)
    microarray = pd.DataFrame([[0., 1.], [1., 2.], [5., 6.], [0., 1.]],
                              index=pd.Series([1, 1, 1, 2], name='label'))
    expected = pd.DataFrame([[2., 3.], [0., 1.]],
                            index=pd.Series([1, 2], name='label'))
    out = samples.groupby_index(microarray)
    pd.testing.assert_frame_equal(out, expected, check_like=True)

    # supplying `labels` appends NaN rows to output
    expected = pd.DataFrame([[2., 3.], [0., 1.], [np.nan, np.nan]],
                            index=pd.Series([1, 2, 3], name='label'))
    out = samples.groupby_index(microarray, labels=[1, 2, 3])
    pd.testing.assert_frame_equal(out, expected, check_like=True)

    # changing `metric` works
    expected = pd.DataFrame([[1., 2.], [0., 1.]],
                            index=pd.Series([1, 2], name='label'))
    out = samples.groupby_index(microarray, metric='median')
    pd.testing.assert_frame_equal(out, expected, check_like=True)


def test_aggregate_samples():
    m1 = pd.DataFrame([[0., 1.], [1., 2.], [5., 6.], [0., 1.]],
                      index=pd.Series([1, 1, 1, 2], name='label'))
    m2 = pd.DataFrame([[0., 1.], [1., 2.], [5., 6.], [10., 11.]],
                      index=pd.Series([1, 1, 2, 3], name='label'))

    # region_agg='donors'
    expected = pd.DataFrame([[1.25, 2.25],
                             [2.5, 3.5],
                             [10., 11.],
                             [np.nan, np.nan]],
                            index=pd.Series([1, 2, 3, 4], name='label'))
    out = samples.aggregate_samples([m1, m2], labels=[1, 2, 3, 4],
                                    region_agg='donors', agg_metric='mean')
    pd.testing.assert_frame_equal(out, expected, check_like=True)

    # region_agg = 'samples'
    expected = pd.DataFrame([[1.4, 2.4],
                             [2.5, 3.5],
                             [10., 11.],
                             [np.nan, np.nan]],
                            index=pd.Series([1, 2, 3, 4], name='label'))
    out = samples.aggregate_samples([m1, m2], labels=[1, 2, 3, 4],
                                    region_agg='samples', agg_metric='mean')
    pd.testing.assert_frame_equal(out, expected, check_like=True)

    # return_donors=True, agg_metric='median'
    expected = [
        pd.DataFrame([[1., 2.], [0., 1.], [np.nan, np.nan], [np.nan, np.nan]],
                     index=pd.Series([1, 2, 3, 4], name='label')),
        pd.DataFrame([[0.5, 1.5], [5, 6], [10, 11], [np.nan, np.nan]],
                     index=pd.Series([1, 2, 3, 4], name='label'))
    ]
    out = samples.aggregate_samples([m1, m2], labels=[1, 2, 3, 4],
                                    return_donors=True, agg_metric='median')
    for e, o in zip(expected, out):
        pd.testing.assert_frame_equal(o, e, check_like=True)

    # check that poorly normalized genes are removed
    m1.loc[2, 1], m2.loc[2, 1] = np.nan, np.nan
    expected = pd.DataFrame([[1.25],
                             [2.5],
                             [10],
                             [np.nan]],
                            index=pd.Series([1, 2, 3, 4], name='label'))
    out = samples.aggregate_samples([m1, m2], labels=[1, 2, 3, 4],
                                    region_agg='donors', agg_metric='mean')
    pd.testing.assert_frame_equal(out, expected, check_like=True)

    # invalid method for region_agg
    with pytest.raises(ValueError):
        samples.aggregate_samples([m1, m2], region_agg='notamethod')
    # region_agg='samples' incompatible with return_donors=True
    with pytest.raises(ValueError):
        samples.aggregate_samples([m1, m2], region_agg='samples',
                                  return_donors=True)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# test all the above functions again on real data to make sure we don't choke #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


def test_update_mni_coords_real(testfiles):
    for annotation in testfiles['annotation']:
        samples.update_mni_coords(annotation)


def test_label_samples_real(testfiles, atlas):
    out = samples.label_samples(testfiles['annotation'][0], atlas['image'])
    assert isinstance(out, pd.DataFrame)
    assert out.index.name == 'sample_id'
    assert out.columns == ['label']


def test_drop_mismatch_samples_real(testfiles):
    for an, on in zip(testfiles['annotation'], testfiles['ontology']):
        samples.drop_mismatch_samples(an, on)


def test_mirror_samples_real(testfiles):
    orig = [363, 470]
    for an, on, o in zip(testfiles['annotation'], testfiles['ontology'], orig):
        out = samples.mirror_samples(an, on)
        # there should be more than the original # of samples but less than or
        # equal to 2x that number (we can't MORE than duplicate)
        assert len(out) > o and len(out) <= o * 2


def test__mirror_ontology_real(testfiles):
    orig = [363, 470]
    for a, o, l in zip(testfiles['annotation'], testfiles['ontology'], orig):
        annot = samples._mirror_ontology(a, o)
        assert len(annot) == l
