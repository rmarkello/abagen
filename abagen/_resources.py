# -*- coding: utf-8 -*-
"""
Locate files that ship inside the abagen package
"""

import os.path as op

_ROOT = op.dirname(op.abspath(__file__))


def resource_filename(package, resource):
    """
    Returns path to `resource` inside the installed `package`

    Replaces ``pkg_resources.resource_filename``, which is not available in
    recent versions of setuptools. Only the ``abagen`` package is supported.
    """

    if package != 'abagen':
        raise ValueError('Only the "abagen" package is supported.')
    return op.join(_ROOT, *resource.split("/"))
