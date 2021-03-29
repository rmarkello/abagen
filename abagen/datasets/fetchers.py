# -*- coding: utf-8 -*-
"""
Functions for downloading the Allen Brain Atlas human microarray dataset.
"""

from collections import namedtuple
from functools import partial
import multiprocessing as mp
import os
from pkg_resources import resource_filename

import nibabel as nib
import pandas as pd

from .. import io
from ..utils import load_gifti, first_entry
from .utils import _get_dataset_dir, _fetch_files

WELL_KNOWN_IDS = nib.volumeutils.Recoder(
    (('9861', 'H0351.2001', '178238387', '157722636', '157722638'),
     ('10021', 'H0351.2002', '178238373', '157723301', '157723303'),
     ('12876', 'H0351.1009', '178238359', '157722290', '157722292'),
     ('15496', 'H0351.1015', '178238266', '162021642', '162021644'),
     ('14380', 'H0351.1012', '178238316', '157721937', '157721939'),
     ('15697', 'H0351.1016', '178236545', '157682966', '157682968')),
    fields=('subj', 'uid', 'url', 't1w', 't2w')
)
VALID_DONORS = sorted(WELL_KNOWN_IDS.value_set('subj')
                      | WELL_KNOWN_IDS.value_set('uid'))
RESOURCE = partial(resource_filename, 'abagen')


def check_donors(donors, default='12876', valid=VALID_DONORS):
    """
    Checks that provided `donors` are valid

    Parameters
    ----------
    donors : list of str
        List of donors to download; can be either donor number or UID. Can also
        specify 'all' to download all available donors. If 'None' is provided
        then `default` will be used.
    default : str, optional
        Default donor to use if `donors` is None. Default: '12876'
    valid : list of str, optional
        List of valid donnor numbers and UIDs. Default: :obj:`VALID_DONORS`

    Returns
    -------
    donors : list of str
        Donor subject IDs
    """

    if donors is None:
        donors = [default]
    elif donors == 'all':
        donors = valid
    elif isinstance(donors, str):
        donors = [donors]

    donors = list(donors).copy()
    for n, sub_id in enumerate(donors):
        if sub_id not in valid:
            raise ValueError('Invalid subject id: {0}. Subjects must in: {1}.'
                             .format(sub_id, valid))
        donors[n] = WELL_KNOWN_IDS[sub_id]  # convert to ID system
    donors = sorted(set(donors), key=lambda x: int(x))

    return donors


def fetch_microarray(data_dir=None, donors=None, resume=True, verbose=1,
                     convert=True, n_proc=1):
    """
    Downloads the Allen Human Brain Atlas microarray expression dataset

    Parameters
    ----------
    data_dir : str, optional
        Directory where data should be downloaded and unpacked. Default: $HOME/
        abagen-data
    donors : list, optional
        List of donors to download; can be either donor number or UID. Can also
        specify 'all' to download all available donors. Default: 12876
    resume : bool, optional
        Whether to resume download of a partly-downloaded file. Default: True
    verbose : int, optional
        Verbosity level (0 means no message). Default: 1
    convert : bool, optional
        Whether to convert downloaded CSV files into parquet format for faster
        loading in the future; only available if ``fastparquet`` and ``python-
        snappy`` are installed. Default: True
    n_proc : int, optional
        Number of processes to parallelize download if multiple donors are
        specified. Default: 1

    Returns
    -------
    data : dict
        Two-level nested dictionary, where top-level keys are donor IDs and
        second-level keys are ['microarray', 'ontology', 'pacall', 'probes',
        'annotation'], where corresponding values are lists of filepaths to
        downloaded CSV files.

    References
    ----------
    Hawrylycz, M. J., Lein, E. S., Guillozet-Bongaarts, A. L., Shen, E. H., Ng,
    L., Miller, J. A., ... & Abajian, C. (2012). An anatomically comprehensive
    atlas of the adult human brain transcriptome. Nature, 489(7416), 391.
    """

    url = "https://human.brain-map.org/api/v2/well_known_file_download/{}"

    dataset_name = 'microarray'
    data_dir = _get_dataset_dir(dataset_name, data_dir=data_dir,
                                verbose=verbose)

    sub_files = ('MicroarrayExpression.csv',
                 'Ontology.csv', 'PACall.csv',
                 'Probes.csv', 'SampleAnnot.csv')
    n_files = len(sub_files)
    donors = check_donors(donors)

    if n_proc < 0:
        n_proc = mp.cpu_count() + n_proc + 1

    files = [
        [(os.path.join('normalized_microarray_donor{}'.format(sub), fname),
            url.format(WELL_KNOWN_IDS.url[sub]),
            dict(uncompress=True,
                 move=os.path.join('normalized_microarray_donor{}'.format(sub),
                                   'donor{}.zip'.format(sub))))
         for fname in sub_files]
        for sub in donors
    ]

    if n_proc > 1:
        with mp.Pool(n_proc) as pool:
            results = [pool.apply_async(_fetch_files,
                                        (data_dir, f),
                                        dict(resume=resume, verbose=verbose))
                       for f in files]
            # flatten outputs into single list
            files = [fn for res in results for fn in res.get()]
    else:
        # flatten list of lists into single list
        files = [fn for f in files for fn in f]
        files = _fetch_files(data_dir, files, resume=resume, verbose=verbose)

    # if we want to convert files to parquet format it's good to do that now
    # this step is _already_ super long, so an extra 1-2 minutes is negligible
    if convert and io.use_parq:
        for fn in files[0::n_files] + files[2::n_files]:
            io._make_parquet(fn, convert_only=True)

    keys = ['microarray', 'ontology', 'pacall', 'probes', 'annotation']
    return {
        donor: dict(zip(keys, files[k:k + n_files]))
        for k, donor in zip(range(0, len(files), n_files), donors)
    }


def fetch_rnaseq(data_dir=None, donors=None, resume=True, verbose=1):
    """
    Downloads RNA-sequencing data from the Allen Human Brain Atlas

    Parameters
    ----------
    data_dir : str, optional
        Directory where data should be downloaded and unpacked. Default:
        current directory
    donors : list, optional
        List of donors to download; can be either donor number or UID. Can also
        specify 'all' to download all available donors (two). Default: 9861
    resume : bool, optional
        Whether to resume download of a partly-downloaded file. Default: True
    verbose : int, optional
        Verbosity level (0 means no message). Default: 1

    Returns
    -------
    data : dict
        Two-level nested dictionary, where top-level keys are donor IDs and
        second-level keys are ['counts', 'tpm', 'ontology', 'genes',
        'annotation'], where corresponding values are lists of filepaths to
        downloaded CSV files.

    References
    ----------
    Hawrylycz, M. J., Lein, E. S., Guillozet-Bongaarts, A. L., Shen, E. H., Ng,
    L., Miller, J. A., ... & Abajian, C. (2012). An anatomically comprehensive
    atlas of the adult human brain transcriptome. Nature, 489(7416), 391.
    """

    url = "https://human.brain-map.org/api/v2/well_known_file_download/{}"
    well_known_ids = {
        '9861': '278447594',
        '10021': '278448166'
    }

    dataset_name = 'rnaseq'
    data_dir = _get_dataset_dir(dataset_name, data_dir=data_dir,
                                verbose=verbose)

    sub_files = ('Genes.csv', 'Ontology.csv',
                 'RNAseqCounts.csv', 'RNAseqTPM.csv', 'SampleAnnot.csv')
    n_files = len(sub_files)
    valid = ['9861', '10021', 'H0351.2001', 'H0351.2002']
    donors = sorted(set(check_donors(donors, default=valid[0])) & set(valid),
                    key=lambda x: int(x))

    files = [
        [(os.path.join('rnaseq_donor{}'.format(sub), fname),
            url.format(well_known_ids[sub]),
            dict(uncompress=True,
                 move=os.path.join('rnaseq_donor{}'.format(sub),
                                   'donor{}.zip'.format(sub))))
         for fname in sub_files]
        for sub in donors
    ]

    files = [fn for f in files for fn in f]
    files = _fetch_files(data_dir, files, resume=resume, verbose=verbose)

    keys = ['genes', 'ontology', 'counts', 'tpm', 'annotation']
    return {
        donor: dict(zip(keys, files[k:k + n_files]))
        for k, donor in zip(range(0, len(files), n_files), donors)
    }


def fetch_raw_mri(data_dir=None, donors=None, resume=True, verbose=1):
    """
    Downloads the "raw" Allen Human Brain Atlas T1w/T2w MRI images

    Parameters
    ----------
    data_dir : str, optional
        Directory where data should be downloaded and unpacked. Default: $HOME/
        abagen-data
    donors : list, optional
        List of donors to download; can be either donor number or UID. Can also
        specify 'all' to download all available donors. Default: 12876
    resume : bool, optional
        Whether to resume download of a partly-downloaded file. Default: True
    verbose : int, optional
        Verbosity level (0 means no message). Default: 1

    Returns
    -------
    mris : dict
        Two-level nested dictionary, where top-level keys are donor IDs and
        second-level keys are ['t1w', 't2w'], where corresponding values are
        lists of filepaths to downloaded Nifti files
    """

    url = "https://human.brain-map.org/api/v2/well_known_file_download/{}"

    dataset_name = 'mri'
    data_dir = _get_dataset_dir(dataset_name, data_dir=data_dir,
                                verbose=verbose)

    sub_files = dict(t1w='T1.nii.gz', t2w='T2.nii.gz')
    n_files = len(sub_files)
    donors = check_donors(donors)

    files = [
        (os.path.join('mri_donor{}'.format(sub), fname),
         url.format(getattr(WELL_KNOWN_IDS, img)[sub]),
         dict(move=os.path.join('mri_donor{}'.format(sub),
                                fname)))
        for sub in donors
        for img, fname in sub_files.items()
    ]

    files = _fetch_files(data_dir, files, resume=resume, verbose=verbose)

    return {
        donor: dict(zip(sub_files.keys(), files[k:k + n_files]))
        for k, donor in zip(range(0, len(files), n_files), donors)
    }


def fetch_freesurfer(data_dir=None, donors=None, resume=True, verbose=1):
    """
    Downloads FreeSurfer reconstructions of the Allen Human Brain Atlas MRIs

    Parameters
    ----------
    data_dir : str, optional
        Directory where data should be downloaded and unpacked. Default: $HOME/
        abagen-data
    donors : list, optional
        List of donors to download; can be either donor number or UID. Can also
        specify 'all' to download all available donors. Default: 12876
    resume : bool, optional
        Whether to resume download of a partly-downloaded file. Default: True
    verbose : int, optional
        Verbosity level (0 means no message). Default: 1

    Returns
    -------
    freesurfer : dict
        Dictionary where keys are donor IDs and values are paths to FreeSurfer
        directories for requested `donors`

    References
    ----------
    Romero-Garcia, R., Whitaker, K., Vasa, F., Seidlitz, J., Shinn, M., Fonagy,
    P., Jones, P., et al. (2017). Data supporting NSPN publication "Structural
    covariance networks are coupled to expression of genes enriched in
    supragranular layers of the human cortex " [Dataset].
    https://doi.org/10.17863/CAM.11392
    """

    url = "https://www.repository.cam.ac.uk/bitstream/handle/1810/265272/" \
          "donor{}.zip"
    dataset_name = 'freesurfer'
    data_dir = _get_dataset_dir(dataset_name, data_dir=data_dir,
                                verbose=verbose)

    donors = check_donors(donors)
    files = [
        ('donor{}'.format(sub),
         url.format(sub),
         dict(uncompress=True,
              move=os.path.join('freesurfer.tar.gz')))
        for sub in donors
    ]

    files = _fetch_files(data_dir, files, resume=resume, verbose=verbose)

    return {
        donor: files[k]
        for k, donor in enumerate(donors)
    }


def fetch_desikan_killiany(native=False, surface=False, *args, **kwargs):
    """
    Fetches Desikan-Killiany atlas shipped with `abagen`

    Parameters
    ----------
    native : bool, optional
        Whether to return individualized atlases in donor native space.
        Default: False
    surface : bool, optional
        Whether to return surface instead of volumetric parcellation. This
        option is currently incompatible with ``native=True``; instead, refer
        to :func:`abagen.datasets.fetch_freesurfer` for donor-specific surface
        atlases. Default: False

    Returns
    -------
    atlas : dict
        Dictionary with keys ['image', 'info'] pointing to atlas image and
        information files. If ``native`` then 'image' is a dictionary where
        keys are donor IDs and values are image paths. If ``surface`` then
        'image' is a tuple of GIFTI files (.label.gii.gz)

    References
    ----------
    Desikan, R. S., Ségonne, F., Fischl, B., Quinn, B. T., Dickerson, B. C.,
    Blacker, D., ... & Albert, M. S. (2006). An automated labeling system
    for subdividing the human cerebral cortex on MRI scans into gyral based
    regions of interest. Neuroimage, 31(3), 968-980.

    Examples
    --------
    >>> import abagen
    >>> atlas = abagen.fetch_desikan_killiany()
    >>> print(atlas['image'])  # doctest: +ELLIPSIS
    /.../abagen/data/atlas-desikankilliany.nii.gz
    >>> print(atlas['info'])  # doctest: +ELLIPSIS
    /.../abagen/data/atlas-desikankilliany.csv

    When fetching native-space atlases, `atlas['image']` will be a dictionary
    where the keys are donor IDs and the values are paths to the donor-specific
    atlases:

    >>> atlas = abagen.fetch_desikan_killiany(native=True)
    >>> print(atlas['image'].keys())
    dict_keys(['9861', '10021', '12876', '14380', '15496', '15697'])
    >>> print(atlas['image']['9861'])  # doctest: +ELLIPSIS
    /.../abagen/data/native_dk/9861/atlas-desikankilliany.nii.gz
    """

    # grab resource filenames
    img = dict()
    for donor in check_donors('all'):
        fp = 'data' if not native else os.path.join('data', 'native_dk', donor)
        if surface:
            impath = tuple([
                RESOURCE(
                    os.path.join(fp, f'atlas-desikankilliany-{h}.label.gii.gz')
                )
                for h in ('lh', 'rh')
            ])
        else:
            impath = RESOURCE(os.path.join(fp, 'atlas-desikankilliany.nii.gz'))
        img[donor] = impath
    if not native:
        img = first_entry(img)
    info = RESOURCE('data/atlas-desikankilliany.csv')

    return dict(image=img, info=info)


def fetch_gene_group(group):
    """
    Return list of gene acronyms belonging to provided `group`

    Groups are defined as in [DS1]_

    Parameters
    ----------
    group : {'brain', 'neuron', 'oligodendrocyte', 'synaptome', 'layers'}
        Desired gene group

    Returns
    -------
    genes : list of str
        List of gene acronyms

    References
    ----------
    .. [DS1] Burt, J. B., Demirtaş, M., Eckner, W. J., Navejar, N. M., Ji, J.
       L., Martin, W. J., ... & Murray, J. D. (2018). Hierarchy of
       transcriptomic specialization across human cortex captured by
       structural neuroimaging topography. Nature neuroscience, 21(9), 1251.
    """

    groups = ['brain', 'neuron', 'oligodendrocyte', 'synaptome', 'layers']
    if group.lower() not in groups:
        raise ValueError('Provided group {} not one of the available gene '
                         'groups: {}'.format(group, groups))

    group = group.lower()
    fn = RESOURCE(os.path.join('data', 'burt2018_natneuro.csv.gz'))
    genes = pd.read_csv(fn).query('group == "{}"'.format(group))['acronym']

    return sorted(list(genes))


def fetch_donor_info():
    """
    Returns dataframe with donor demographic information

    Returns
    -------
    info : pandas.DataFrame
        With columns ['donor', 'age', 'sex', 'ethnicity', 'medical_conditions',
        'post_mortem_interval_hours'] detailing basic demographic info about
        donors
    """

    fn = RESOURCE(os.path.join('data', 'donor_info.csv'))
    donors = pd.read_csv(fn)

    return donors


Brain = namedtuple('Brain', ('lh', 'rh'))
Surface = namedtuple('Surface', ('vertices', 'faces'))


def fetch_fsaverage5(load=True):
    """
    Fetches and optionally loads fsaverage5 surface

    Parameters
    ----------
    load : bool, optional
        Whether to pre-load files. Default: True

    Returns
    -------
    brain : namedtuple ('lh', 'rh')
        If `load` is True, a namedtuple where each entry in the tuple is a
        hemisphere, represented as a namedtuple with fields ('vertices',
        'faces'). If `load` is False, a namedtuple where entries are filepaths.
    """

    hemispheres = []
    for hemi in ('lh', 'rh'):
        fn = RESOURCE(
            os.path.join('data', f'fsaverage5-pial-{hemi}.surf.gii.gz')
        )
        if load:
            hemispheres.append(Surface(*load_gifti(fn).agg_data()))
        else:
            hemispheres.append(fn)

    return Brain(*hemispheres)


def fetch_fsnative(donors, surf='pial', load=True, data_dir=None, resume=True,
                   verbose=1):
    """
    Fetches and optionally loads fsnative surface of `donor`

    Parameters
    ----------
    donors : str or list-of-str
        Donor(s) to download; can be either donor number or UID. Can also
        specify 'all' to download all available donors.
    surf : {'orig', 'white', 'pial', 'inflated', 'sphere'}, optional
        Which surface to load. Default: 'pial'
    load : bool, optional
        Whether to pre-load files. Default: True
    data_dir : str, optional
        Directory where data should be downloaded and unpacked. Default: $HOME/
        abagen-data
    resume : bool, optional
        Whether to resume download of a partly-downloaded file. Default: True
    verbose : int, optional
        Verbosity level (0 means no message). Default: 1

    Returns
    -------
    brain : namedtuple ('lh', 'rh')
        If `load` is True, a namedtuple where each entry in the tuple is a
        hemisphere, represented as a namedtuple with fields ('vertices',
        'faces'). If `load` is False, a namedtuple where entries are filepaths.
        If multiple donors are requested a dictionary is returned where keys
        are donor IDs.
    """

    donors = check_donors(donors)
    if len(donors) > 1:
        return {donor: fetch_fsnative(donor, surf, data_dir, resume, verbose)
                for donor in donors}

    donors = donors[0]
    fpath = fetch_freesurfer(donors=donors, data_dir=data_dir, resume=resume,
                             verbose=verbose)[donors]
    hemispheres = []
    for hemi in ('lh', 'rh'):
        fn = os.path.join(fpath, 'surf', f'{hemi}.{surf}')
        if load:
            hemispheres.append(Surface(*nib.freesurfer.read_geometry(fn)))
        else:
            hemispheres.append(fn)

    return Brain(*hemispheres)
