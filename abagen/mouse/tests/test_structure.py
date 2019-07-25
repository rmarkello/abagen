# -*- coding: utf-8 -*-

import random
import numpy as np
import pytest
from abagen.mouse import io, structure

RANDOM = random.Random(1234)
STRUCTS = io.fetch_allenref_structures()
COORD_ID = {
    'sagittal': (7800, 3400, 1050),
    'coronal': (7800, 3400, 10350)
}
COORD_ACRONYM = {
    'sagittal': (5740, 2060, 2730),
    'coronal': (5740, 2060, 8670)
}


def test_get_structure_info():
    samples = STRUCTS.loc[RANDOM.sample(range(len(STRUCTS)), 10)]

    acronyms = structure.get_structure_info(id=samples['id'],
                                            attributes='acronym')
    assert sorted(acronyms['acronym']) == sorted(samples['acronym'])

    names = structure.get_structure_info(id=samples['id'],
                                         attributes='name')
    assert sorted(names['name']) == sorted(samples['name'])

    with pytest.raises(ValueError):
        structure.get_structure_info(acronym=samples['acronym'],
                                     attributes='invalid')

    # multiple attributes
    info = structure.get_structure_info(id=samples['id'],
                                        attributes=['name', 'acronym'])
    assert sorted(acronyms['acronym']) == sorted(samples['acronym'])
    assert sorted(info['name']) == sorted(samples['name'])

    # all attributes
    info = structure.get_structure_info(id=samples['id'])
    assert len(info.columns) == 16

    # exception: invalid gene identifiers
    with pytest.raises(ValueError):
        structure.get_structure_info(id=-100000)
    with pytest.raises(ValueError):
        structure.get_structure_info(acronym='notastructure')


@pytest.mark.parametrize('space', ['coronal', 'sagittal'])
def test_get_structure_coordinates(space):
    info = structure.get_structure_coordinates(id=1018,
                                               reference_space=space)
    assert tuple(np.asarray(info[['x', 'y', 'z']])[0]) == COORD_ID[space]

    info = structure.get_structure_coordinates(acronym='SSp',
                                               reference_space=space)
    assert tuple(np.asarray(info[['x', 'y', 'z']])[0]) == COORD_ACRONYM[space]

    # exception: structure is invalid
    with pytest.raises(ValueError):
        structure.get_structure_coordinates(id=-1000000)
    with pytest.raises(ValueError):
        structure.get_structure_coordinates(acronym='notastructure')
