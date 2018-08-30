# -*- coding: utf-8 -*-

__version__ = '0.0.1'

NAME = 'abagen'
MAINTAINER = 'Ross Markello'
EMAIL = 'rossmarkello@gmail.com'
VERSION = __version__
LICENSE = 'MIT'
DESCRIPTION = """\
A toolbox for working with the Allen Brain Atlas human microarray datasets\
"""
LONG_DESCRIPTION = """\
"""
URL = 'https://github.com/rmarkello/{name}'.format(name=NAME)
DOWNLOAD_URL = ('http://github.com/rmarkello/{name}/archive/{ver}.tar.gz'
                .format(name=NAME, ver=__version__))

INSTALL_REQUIRES = [
    'nibabel',
    'nilearn',
    'numpy',
    'pandas',
    'scikit-learn',
    'scipy',
]

EXTRAS_REQUIRE = {
    'io': ['fastparquet', 'snappy']
}

TESTS_REQUIRE = [
    'pytest',
]

PACKAGE_DATA = {
}
