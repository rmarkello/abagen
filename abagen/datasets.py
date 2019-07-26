# -*- coding: utf-8 -*-
"""
Functions for downloading Allen Brain Atlas human microarray dataset

Modeled after ``nilearn.datasets``. Currently just downloads into current
working directory, but will likely be modified to download into a more
"standard" directory.
"""

import os
from pkg_resources import resource_filename
from nibabel.volumeutils import Recoder
from nilearn.datasets.utils import _fetch_files
import pandas as pd
from sklearn.utils import Bunch
from abagen import io

WELL_KNOWN_IDS = Recoder(
    (('9861', 'H0351.2001', '178238387', '157722636'),
     ('10021', 'H0351.2002', '178238373', '157723301'),
     ('12876', 'H0351.1009', '178238359', '157722290'),
     ('15496', 'H0351.1015', '178238266', '162021642'),
     ('14380', 'H0351.1012', '178238316', '157721937'),
     ('15697', 'H0351.1016', '178236545', '157682966')),
    fields=('subj', 'uid', 'url', 't1w',)
)

VALID_DONORS = sorted(WELL_KNOWN_IDS.value_set('subj') |
                      WELL_KNOWN_IDS.value_set('uid'))


def _get_dataset_dir(dataset_name, data_dir=None, verbose=1):
    """
    Gets path to `dataset_name`

    Parameters
    ----------
    dataset_name : str
        The name of the dataset in question
    data_dir : str, optional
        Path to use as data directory. If not specified, will check for
        environmental variables 'ABAGEN_DATA'; if that is not set, will use
        '~/abagen-data' instead. Default: None
    verbose : int, optional
        Verbosity level (0 means no message). Default: 1

    Returns
    -------
    data_dir : str
        Path to use as data directory

    References
    ----------
    Function lightly modified from `nilearn <https://github.com/nilearn/
    nilearn>`_, which is licensed under the BSD, referenced here:

    Copyright (c) 2007 - 2015 The nilearn developers.
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    a. Redistributions of source code must retain the above copyright notice,
        this list of conditions and the following disclaimer.
    b. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
    c. Neither the name of the nilearn developers nor the names of
        its contributors may be used to endorse or promote products
        derived from this software without specific prior written
        permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR
    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
    LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
    OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
    DAMAGE.
    """

    paths = [
            os.environ.get('ABAGEN_DATA'),
            os.path.join('~', 'abagen-data'),
            os.getcwd()
    ]
    if data_dir is not None:
        paths = [data_dir, os.path.dirname(data_dir)] + paths
    paths = [os.path.expanduser(d) for d in paths if d is not None]

    if verbose > 2:
        print('Dataset search paths: {}'.format(paths))

    for path in paths:
        if not os.path.basename(path) == dataset_name:
            path = os.path.join(path, dataset_name)
        if os.path.islink(path):
            link = os.readlink(path)
            if os.path.isabs(link):
                path = link
            path = os.path.join(os.path.dirname(path), link)
        if os.path.exists(path) and os.path.isdir(path):
            if verbose > 1:
                print('\nDataset found in {}\n'.format(path))
            return path

    errors = []
    for path in paths:
        if not os.path.basename(path) == dataset_name:
            path = os.path.join(path, dataset_name)
        if not os.path.exists(path):
            try:
                os.makedirs(path)
                if verbose > 0:
                    print('\nDataset created in {}\n'.format(path))
                return path
            except (NotADirectoryError, PermissionError) as exc:
                err_msg = getattr(exc, 'strerror', str(exc))
                errors.append('\n -{} ({})'.format(path, err_msg))

    raise OSError('Tried to store dataset {} in the following directories, '
                  'but: ' + ''.join(errors))


def fetch_microarray(data_dir=None, donors=['9861'], resume=True, verbose=1,
                     convert=True):
    """
    Downloads the Allen Human Brain Atlas microarray expression dataset

    Parameters
    ----------
    data_dir : str, optional
        Directory where data should be downloaded and unpacked. Default:
        current directory
    donors : list, optional
        List of donors to download; can be either donor number or UID. Can also
        specify 'all' to download all available donors. Default: 9861
    resume : bool, optional
        Whether to resume download of a partly-downloaded file. Default: True
    verbose : int, optional
        Verbosity level (0 means no message). Default: 1
    convert : bool, optional
        Whether to convert downloaded CSV files into parquet format for faster
        loading in the future; only available if ``fastparquet`` and ``python-
        snappy`` are installed. Default: True

    Returns
    -------
    data : :class:`sklearn.utils.Bunch`
        Dictionary-like object with keys ['microarray', 'ontology', 'pacall',
        'probes', 'annotation'], where corresponding values are lists of
        filepaths to downloaded CSV files.

    References
    ----------
    Hawrylycz, M. J., Lein, E. S., Guillozet-Bongaarts, A. L., Shen, E. H., Ng,
    L., Miller, J. A., ... & Abajian, C. (2012). An anatomically comprehensive
    atlas of the adult human brain transcriptome. Nature, 489(7416), 391.
    """

    url = "https://human.brain-map.org/api/v2/well_known_file_download/{}"

    dataset_name = 'allenbrain'
    data_dir = _get_dataset_dir(dataset_name, data_dir=data_dir,
                                verbose=verbose)

    sub_files = ('MicroarrayExpression.csv',
                 'Ontology.csv', 'PACall.csv',
                 'Probes.csv', 'SampleAnnot.csv')
    n_files = len(sub_files)

    if donors is not None and (isinstance(donors, (list, tuple))):
        for n, sub_id in enumerate(donors):
            if sub_id not in VALID_DONORS:
                raise ValueError('You provided invalid subject id {0} in a'
                                 'list. Subjects must be selected in {1}.'
                                 .format(sub_id, VALID_DONORS))
            donors[n] = WELL_KNOWN_IDS[sub_id]  # convert to ID system
    elif donors == 'all':
        donors = WELL_KNOWN_IDS.value_set('subj')
    else:
        donors = []
    donors = sorted(set(donors), key=lambda x: int(x))  # avoid duplicates

    files = [
         (os.path.join('normalized_microarray_donor{}'.format(sub), fname),
          url.format(WELL_KNOWN_IDS.url[sub]),
          dict(uncompress=True,
               move=os.path.join('normalized_microarray_donor{}'.format(sub),
                                 'donor{}.zip'.format(sub))))
         for sub in donors
         for fname in sub_files
    ]

    files = _fetch_files(data_dir, files, resume=resume, verbose=verbose)

    # if we want to convert files to parquet format it's good to do that now
    # this step is _already_ super long, so an extra 1-2 minutes is negligible
    if convert and io.use_parq:
        for fn in files[0::n_files] + files[2::n_files]:
            io._make_parquet(fn, convert_only=True)

    return Bunch(
        microarray=files[0::n_files],
        ontology=files[1::n_files],
        pacall=files[2::n_files],
        probes=files[3::n_files],
        annotation=files[4::n_files]
    )


def fetch_mri(data_dir=None, donors=['9861'], resume=True, verbose=1):
    """
    Downloads the Allen Human Brain Atlas T1w MRI images

    Parameters
    ----------
    data_dir : str, optional
        Directory where data should be downloaded and unpacked. Default:
        current directory
    donors : list, optional
        List of donors to download; can be either donor number or UID.
        Default: donor9861
    resume : bool, optional
        Whether to resume download of a partly-downloaded file. Default: True
    verbose : int, optional
        Verbosity level (0 means no message). Default: 1
    """

    raise NotImplementedError


def _fetch_alleninf_coords(*args, **kwargs):
    """
    Gets updated MNI coordinates for AHBA samples, as shipped with `alleninf`

    Returns
    -------
    coords : :class:`pandas.DataFrame`
        Updated MNI coordinates for all AHBA samples

    References
    ----------
    Updated MNI coordinates taken from https://github.com/chrisfilo/alleninf,
    which is licensed under the BSD-3 (reproduced here):

    Copyright (c) 2018, Krzysztof Gorgolewski
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

    * Neither the name of the copyright holder nor the names of its
    contributors may be used to endorse or promote products derived from
    this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
    LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.
    """

    coords = resource_filename('abagen', 'data/corrected_mni_coordinates.csv')
    coords = pd.read_csv(coords).rename(dict(corrected_mni_x='mni_x',
                                             corrected_mni_y='mni_y',
                                             corrected_mni_z='mni_z'),
                                        axis=1)
    return coords.set_index('well_id')


def fetch_desikan_killiany(*args, **kwargs):
    """
    Fetches Desikan-Killiany atlas shipped with `abagen`

    Returns
    -------
    atlas : :class:`sklearn.utils.Bunch`
        Dictionary-like object with attributes ['image', 'info'] pointing to
        atlas image (.nii.gz) and information (.csv) files

    References
    ----------
    Desikan, R. S., Ségonne, F., Fischl, B., Quinn, B. T., Dickerson, B. C.,
    Blacker, D., ... & Albert, M. S. (2006). An automated labeling system
    for subdividing the human cerebral cortex on MRI scans into gyral based
    regions of interest. Neuroimage, 31(3), 968-980.
    """
    # grab resource filenames
    image = resource_filename('abagen', 'data/atlas-desikankilliany.nii.gz')
    info = resource_filename('abagen', 'data/atlas-desikankilliany.csv')

    return Bunch(image=image, info=info)
