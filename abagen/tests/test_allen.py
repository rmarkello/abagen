import pandas as pd
from nilearn._utils import check_niimg
from nilearn.image import new_img_like
import pytest
from abagen import allen
from abagen.datasets import fetch_desikan_killiany

ATLAS = fetch_desikan_killiany()


def test_label_samples(testfiles):
    out = allen.label_samples(testfiles.annotation[0], ATLAS.image)
    assert isinstance(out, pd.DataFrame)
    assert out.index.name == 'sample_id'
    assert out.columns == ['label']


def test_vanilla_get_expression_data(testfiles):
    out = allen.get_expression_data(testfiles, ATLAS.image)
    assert isinstance(out, pd.DataFrame)
    assert out.index.name == 'label'
    assert out.columns.name == 'gene_symbol'

    with pytest.raises(KeyError):
        allen.get_expression_data({'microarray': [1, 2, 3]}, ATLAS.image)


def test_extra_get_expression_data(testfiles):
    for opts in [{'atlas_info': ATLAS.info},
                 {'exact': False},
                 {'reannotated': False},
                 {'atlas_info': ATLAS.info, 'exact': False}]:
        out = allen.get_expression_data(testfiles, ATLAS.image, **opts)
        assert isinstance(out, pd.DataFrame)
        assert out.index.name == 'label'
        assert out.columns.name == 'gene_symbol'


def test_missing_labels(testfiles):
    # remove some labels from atlas image so numbers are non-sequential
    remove = [10, 20, 60]
    # subset atlas image
    atlas = check_niimg(ATLAS.image).get_data()
    for i in remove:
        atlas[atlas == i] = 0
    atlas = new_img_like(ATLAS.image, atlas)
    # subset atlas info
    atlas_info = pd.read_csv(ATLAS.info)
    atlas_info = atlas_info[~atlas_info.id.isin(remove)]
    # test get expression
    out, counts = allen.get_expression_data(testfiles, atlas, atlas_info,
                                            exact=False, return_counts=True)
    assert isinstance(out, pd.DataFrame)
    assert out.index.name == 'label'
    assert out.columns.name == 'gene_symbol'
    assert len(out) == len(atlas_info)

    assert isinstance(counts, pd.DataFrame)
    assert counts.shape == (len(atlas_info), len(testfiles.probes))
