__all__ = ['__version__', '__doc__', 'io', 'get_expression_data',
           'keep_stable_genes', 'remove_distance', 'aggregate_donors',
           'fetch_microarray', 'fetch_desikan_killiany', 'mouse']

from ._version import get_versions
version = get_versions()['version']
del get_versions

from .info import long_description as __doc__

from . import io, mouse
from .allen import get_expression_data
from .correct import keep_stable_genes, remove_distance
from .datasets import fetch_microarray, fetch_desikan_killiany
from .process import aggregate_donors
