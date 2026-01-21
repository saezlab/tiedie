"""TieDIE: Tied Diffusion for Network Discovery

A network analysis algorithm that finds subnetworks connecting genomic
perturbations to transcriptional changes in large gene interaction networks.
"""

from .ppr import PprDiffuser
from .util import (
    parse_net,
    parse_heats,
    filter_linkers,
    normalize_heats,
    connected_subnets,
    get_network_nodes,
    find_linker_cutoff,
    map_ugraph_to_network,
)
from .kernel import Kernel
from .permute import NetBalancedPermuter
from ._metadata import __author__, __license__, __version__
from .kernel_scipy import ScipyKernel


__all__ = [
    '__version__',
    '__author__',
    '__license__',
    'Kernel',
    'ScipyKernel',
    'PprDiffuser',
    'NetBalancedPermuter',
    'parse_heats',
    'parse_net',
    'normalize_heats',
    'get_network_nodes',
    'filter_linkers',
    'find_linker_cutoff',
    'connected_subnets',
    'map_ugraph_to_network',
]
