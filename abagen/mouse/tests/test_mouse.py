# -*- coding: utf-8 -*-

import numpy as np
import pytest
from abagen.mouse import mouse

STRUCTURES = [182305713, 182305709, 182305705]
ATTRIBUTES = ['expression_energy', 'expression_density', 'sum_pixels']
EXPERIMENTS = {
    986: {
        'expression_energy': np.array([7.73432, 7.28206, 3.82741]),
        'expression_density': np.array([0.0603072, 0.0553335, 0.0298628]),
        'sum_pixels': np.array([419628.0, 2238000.0, 3629050.0])
    },
    69782969: {
        'expression-energy': np.array([1.00879, 1.40198, 2.34988]),
        'expression-density': np.array([0.00650952, 0.00896709, 0.015409]),
        'sum-pixels': np.array([195863.0, 1051230.0, 1440430.0])
    }
}
GBA_UNIONIZATION = {
    'expression_energy': np.array([0.236301, 0.261266, 0.416281]),
    'expression_density': np.array([0.00207073, 0.00232453, 0.00372069]),
    'sum_pixels': np.array([419628.0, 2238000.0, 3636670.0]),
    'voxel_energy_cv': np.array([0.923913, 0.995631, 0.81265])
}


@pytest.mark.parametrize(('genes', 'direction', 'expected'), [
    (dict(acronym='Snca'), 'sagittal', [988, 990, 79904550]),
    (dict(acronym='Snca'), 'coronal', [986, 79908848]),
    (dict(acronym='Gba'), 'sagittal', [1612]),
    (dict(acronym='Gba'), 'coronal', [1611]),
    (dict(acronym='Elf4'), 'sagittal', [73834415, 77464840]),
    (dict(id=84193), 'sagittal', [70238925, 71147924, 71213117]),
    (dict(id=84193), 'coronal', [74047443]),
    (dict(id=18608), 'sagittal', [69289721]),
])
def test_get_experiments_from_gene(genes, direction, expected):
    with pytest.raises(ValueError):
        mouse._get_experiments_from_gene(acronym='notagene',
                                         slicing_direction='coronal')
    with pytest.raises(ValueError):
        mouse._get_experiments_from_gene(id=84193,
                                         slicing_direction='notadirection')
    with pytest.raises(ValueError):
        mouse._get_experiments_from_gene(id=-1000000)

    # get experiments from provided input
    exp = mouse._get_experiments_from_gene(**genes,
                                           slicing_direction=direction)
    assert len(exp) == len(expected)
    assert all(i == j for (i, j) in zip(sorted(exp), sorted(expected)))


@pytest.mark.parametrize(('experiment', 'attributes'), [
    (986, None), (986, ATTRIBUTES), (986, 'all'),
    (69782969, None), (69782969, ATTRIBUTES), (69782969, 'all')
])
def test_get_unionization_from_experiment(experiment, attributes):
    with pytest.raises(ValueError):
        mouse._get_unionization_from_experiment(-100, structures=STRUCTURES)

    # get data from provided experiment
    data = mouse._get_unionization_from_experiment(experiment,
                                                   structures=STRUCTURES,
                                                   attributes=attributes)
    data.index = data.index.droplevel('gene_id')

    if attributes is None:
        attributes = ['expression_density']
    elif attributes == 'all':
        attributes = mouse._UNIONIZATION_ATTRIBUTES

    assert len(data.columns) == len(attributes)
    for attr in set(EXPERIMENTS[experiment].keys()).intersection(attributes):
        assert np.allclose(np.asarray(data.loc[STRUCTURES, attr]),
                           EXPERIMENTS[experiment][attr])


@pytest.mark.parametrize(('attribute'), [
    'expression_energy', 'expression_density', 'sum_pixels', 'voxel_energy_cv'
])
def test_get_unionization_from_gene(attribute):
    with pytest.raises(ValueError):
        mouse.get_unionization_from_gene(acronym='Gba',
                                         slicing_direction='notadirection',
                                         structures=STRUCTURES)
    with pytest.raises(ValueError):
        mouse.get_unionization_from_gene(id=18608,
                                         slicing_direction='coronal',
                                         structures=STRUCTURES)

    # get data for provided gene
    data = mouse.get_unionization_from_gene(acronym='Gba',
                                            slicing_direction='coronal',
                                            structures=STRUCTURES,
                                            attributes=attribute)
    data.index = data.index.droplevel('gene_id')

    assert np.allclose(data.loc[STRUCTURES, attribute],
                       GBA_UNIONIZATION[attribute])
