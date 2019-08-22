__all__ = ['__version__', 'io', 'get_expression_data',
           'keep_stable_genes', 'remove_distance', 'aggregate_donors',
           'fetch_microarray', 'fetch_desikan_killiany', 'mouse']

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

from abagen.allen import get_expression_data
from abagen.correct import keep_stable_genes, remove_distance
from abagen.datasets import fetch_microarray, fetch_desikan_killiany
from abagen.process import aggregate_donors
from abagen import io
from abagen import mouse
