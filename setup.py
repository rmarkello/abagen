#!/usr/bin/env python
import sys

from setuptools import setup
import versioneer

SETUP_REQUIRES = ['setuptools >= 30.3.0']
SETUP_REQUIRES += ['wheel'] if 'bdist_wheel' in sys.argv else []

if __name__ == "__main__":
    setup(name='abagen',
          setup_requires=SETUP_REQUIRES,
          version=versioneer.get_version(),
          cmdclass=versioneer.get_cmdclass())
