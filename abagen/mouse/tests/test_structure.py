import pytest
from ..structure import (check_structure_validity,
                         get_structure_info,
                         get_structure_coordinates)
import random
from ..io import read_all_structures

# number of tests to make
TEST_COUNT = 10
STRUCTURES_LIST_ACRONYM = read_all_structures(entry_type='acronym')
STRUCTURES_LIST_ID = read_all_structures(entry_type='id')
STRUCTURES_LIST_NAME = read_all_structures(entry_type='name')
RANDOM_STRING = 'random_string'
RANDOM_ID = 000000000

TEST_COOR_STRUCTURE_ID = 1018
COOR_ID = [(7800, 3400, 1050), (7800, 3400, 10350)]

TEST_COOR_STRUCTURE_ACRONYM = 'SSp'
COOR_ACRONYM = [(5740, 2060, 2730), (5740, 2060, 8670)]

# sample 10 structures
TEST_SAMPLES = random.sample(
    range(len(STRUCTURES_LIST_ACRONYM)), TEST_COUNT
)


def test_check_structure_validity():
    for sample in TEST_SAMPLES:
        # acronym is given
        validity, root = check_structure_validity(
            acronym=STRUCTURES_LIST_ACRONYM[sample]
        )
        assert validity is True
        # structure id is given
        validity, root = check_structure_validity(
            structure_id=STRUCTURES_LIST_ID[sample]
        )
        assert validity is True
        # structure name is given
        validity, root = check_structure_validity(
            name=STRUCTURES_LIST_NAME[sample]
        )
        assert validity is True

    # exceptions: structure id is invalid
    with pytest.raises(ValueError):
        check_structure_validity(structure_id=RANDOM_ID)
    # exceptions: structure acronym is invalid
    with pytest.raises(ValueError):
        check_structure_validity(acronym=RANDOM_STRING)
    # exceptions: structure id is invalid
    with pytest.raises(ValueError):
        check_structure_validity(name=RANDOM_STRING)


def test_get_structure_info():
    for sample in TEST_SAMPLES:
        # single attribute
        structure_info = get_structure_info(
            structure_id=STRUCTURES_LIST_ID[sample],
            attributes='acronym'
        )
        # structure_info is str
        assert structure_info == STRUCTURES_LIST_ACRONYM[sample]
        structure_info = get_structure_info(
            acronym=STRUCTURES_LIST_ACRONYM[sample],
            attributes='id'
        )
        # structure_infor is int
        assert structure_info == STRUCTURES_LIST_ID[sample]
        structure_info = get_structure_info(
            structure_id=STRUCTURES_LIST_ID[sample],
            attributes='name'
        )
        # structure_info is str
        assert structure_info == STRUCTURES_LIST_NAME[sample]
        # no such attribute, return an empty structure_info
        structure_info = get_structure_info(
            structure_id=STRUCTURES_LIST_ID[sample],
            attributes=RANDOM_STRING
        )
        # structure_info is str
        assert not structure_info
        structure_info = get_structure_info(
            acronym=STRUCTURES_LIST_ACRONYM[sample],
            attributes=RANDOM_STRING
        )
        # structure_info is str
        assert not structure_info
        # multiple attributes
        # attributes = 'all'
        structure_info = get_structure_info(
            structure_id=STRUCTURES_LIST_ID[sample],
        )
        assert structure_info['acronym'] == STRUCTURES_LIST_ACRONYM[sample]
        assert structure_info['name'] == STRUCTURES_LIST_NAME[sample]
        structure_info = get_structure_info(
            acronym=STRUCTURES_LIST_ACRONYM[sample],
            attributes=['id', 'name', RANDOM_STRING]
        )
        assert RANDOM_STRING not in structure_info
        assert structure_info['id'] == STRUCTURES_LIST_ID[sample]
        assert structure_info['name'] == STRUCTURES_LIST_NAME[sample]

    # exceptions: structure is invalid
    with pytest.raises(ValueError):
        get_structure_info(
            structure_id=RANDOM_ID
        )
    with pytest.raises(ValueError):
        get_structure_info(
            acronym=RANDOM_STRING,
            attributes=['id', 'name']
        )


def test_get_structure_center():
    coor = get_structure_coordinates(
        structure_id=TEST_COOR_STRUCTURE_ID,
    )
    # coor is a list
    assert coor[0] == COOR_ID[0] and coor[1] == COOR_ID[1]

    coor = get_structure_coordinates(
        acronym='SSp',
    )
    assert coor[0] == COOR_ACRONYM[0] and coor[1] == COOR_ACRONYM[1]

    # exception: structure is invalid
    with pytest.raises(ValueError):
        get_structure_coordinates(
            structure_id=RANDOM_ID,
        )
    with pytest.raises(ValueError):
        get_structure_coordinates(
            acronym=RANDOM_STRING
        )
