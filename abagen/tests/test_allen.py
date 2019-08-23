import pandas as pd
from abagen import allen
from abagen.datasets import fetch_desikan_killiany
from abagen.utils import check_img

ATLAS = fetch_desikan_killiany()


def test_label_samples(testfiles):
    out = allen.label_samples(testfiles['annotation'][0], ATLAS['image'])
    assert isinstance(out, pd.DataFrame)
    assert out.index.name == 'sample_id'
    assert out.columns == ['label']


def test_vanilla_get_expression_data(testdir, testfiles):
    out = allen.get_expression_data(ATLAS['image'], data_dir=testdir,
                                    donors=['12876', '15496'])
    assert isinstance(out, pd.DataFrame)
    assert out.index.name == 'label'
    assert out.columns.name == 'gene_symbol'


def test_extra_get_expression_data(testdir, testfiles):
    for opts in [{'atlas_info': ATLAS['info']},
                 {'exact': False},
                 {'reannotated': False},
                 {'atlas_info': ATLAS['info'], 'exact': False}]:
        out = allen.get_expression_data(ATLAS['image'], data_dir=testdir,
                                        donors=['12876', '15496'], **opts)
        assert isinstance(out, pd.DataFrame)
        assert out.index.name == 'label'
        assert out.columns.name == 'gene_symbol'


def test_missing_labels(testdir, testfiles):
    # remove some labels from atlas image so numbers are non-sequential
    remove = [10, 20, 60]
    # subset atlas image
    atlas = check_img(ATLAS['image']).get_data()
    for i in remove:
        atlas[atlas == i] = 0
    atlas = ATLAS['image'].__class__(atlas, ATLAS['image'].affine)
    # subset atlas info
    atlas_info = pd.read_csv(ATLAS['info'])
    atlas_info = atlas_info[~atlas_info.id.isin(remove)]
    # test get expression
    out, counts = allen.get_expression_data(atlas, atlas_info,
                                            exact=False, return_counts=True,
                                            data_dir=testdir,
                                            donors=['12876', '15496'])
    assert isinstance(out, pd.DataFrame)
    assert out.index.name == 'label'
    assert out.columns.name == 'gene_symbol'
    assert len(out) == len(atlas_info)

    assert isinstance(counts, pd.DataFrame)
    assert counts.shape == (len(atlas_info), len(testfiles['probes']))
