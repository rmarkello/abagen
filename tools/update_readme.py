#!/usr/bin/env python
"""
Refresh README.rst file from long description

Should be run from abagen root (containing setup.py)

Authored by: @effigies and @matthew-brett for `nibabel`, licensed under MIT
"""

import os
import runpy

readme_lines = []
with open('README.rst', 'rt') as fobj:
    for line in fobj:
        readme_lines.append(line)
        if line.startswith('.. Following contents should be'):
            break
    else:
        raise ValueError('Expected comment not found')

rel = runpy.run_path(os.path.join('abagen', 'info.py'))

readme = ''.join(readme_lines) + '\n' + rel['long_description']

with open('README.rst', 'wt') as fobj:
    fobj.write(readme)

print('Done')
