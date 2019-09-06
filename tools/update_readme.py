#!/usr/bin/env python
"""
Refresh README.rst file from long description

Should be run from abagen root (containing setup.py)

Authored by: @effigies and @matthew-brett for `nibabel`, licensed under MIT
(https://github.com/nipy/nibabel/blob/master/COPYING)
"""

import os
import re
import runpy

readme_lines = []
with open('README.rst', 'rt') as fobj:
    for line in fobj:
        readme_lines.append(line)
        if line.startswith('.. Following contents should be'):
            break
    else:
        raise ValueError('Expected comment not found')

# strip doctest markers out of README
rel = runpy.run_path(os.path.join('abagen', 'info.py'))['long_description']
rel = re.sub(r'(\s)*\# doctest: +(\S)*\n', '\n', rel)

readme = ''.join(readme_lines) + '\n' + rel

with open('README.rst', 'wt') as fobj:
    fobj.write(readme)

print('Done')
