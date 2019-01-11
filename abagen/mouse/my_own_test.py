from ..structure import check_structure_validity
from ..io import read_all_structures
import random

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

if __name__ == "__main__":
    for sample in TEST_SAMPLES:
        validity, root1, root2 = check_structure_validity(
            acronym=STRUCTURES_LIST_ACRONYM[sample]
        )
        print(validity)