# -*- coding: utf-8 -*-

__version__ = '0.0.1'

NAME = 'abagen'
MAINTAINER = 'Ross Markello'
EMAIL = 'rossmarkello@gmail.com'
VERSION = __version__
LICENSE = 'BSD-3'
DESCRIPTION = """\
A toolbox for working with the Allen Brain Atlas human genetic data\
"""
LONG_DESCRIPTION = 'README.rst'
LONG_DESCRIPTION_CONTENT_TYPE = 'text/x-rst'
URL = 'https://github.com/rmarkello/{name}'.format(name=NAME)
DOWNLOAD_URL = ('http://github.com/rmarkello/{name}/archive/{ver}.tar.gz'
                .format(name=NAME, ver=__version__))

INSTALL_REQUIRES = [
    'nibabel',
    'nilearn',
    'numpy',
    'pandas',
    'requests',
    'scikit-learn',
    'scipy',
]

TESTS_REQUIRES = [
    'codecov',
    'pytest',
    'pytest-cov'
]

EXTRAS_REQUIRES = {
    'doc': [
        'sphinx>=1.2',
        'sphinx_rtd_theme'
    ],
    'io': [
        'fastparquet',
        'python-snappy'
    ],
    'tests': TESTS_REQUIRES
}

EXTRAS_REQUIRES['all'] = list(set([
    v for deps in EXTRAS_REQUIRES.values() for v in deps
]))


PACKAGE_DATA = {
    'abagen': [
        'data/*'
        ]
}

CLASSIFIERS = [
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: BSD License',
    'Programming Language :: Python :: 3.6',
]
