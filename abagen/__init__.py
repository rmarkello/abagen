__all__ = [
    '__version__', '__doc__',
    'io', 'mouse', 'get_expression_data', 'keep_stable_genes',
    'normalize_expression', 'remove_distance', 'fetch_desikan_killiany',
    'fetch_gene_group', 'fetch_microarray', 'fetch_raw_mri'
]

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

from .info import long_description as __doc__

from . import io, mouse
from .allen import get_expression_data
from .correct import keep_stable_genes, normalize_expression, remove_distance
from .datasets import (fetch_desikan_killiany, fetch_gene_group,
                       fetch_microarray, fetch_raw_mri)
