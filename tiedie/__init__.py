"""
TieDIE: Tied Diffusion for Network Discovery

A network analysis algorithm that finds subnetworks connecting genomic
perturbations to transcriptional changes in large gene interaction networks.
"""

from ._metadata import __author__, __license__, __version__
from .kernel import Kernel
from .kernel_scipy import SciPYKernel
from .permute import NetBalancedPermuter
from .ppr import PPrDiffuser
from .util import (
    connectedSubnets,
    filterLinkers,
    findLinkerCutoff,
    getNetworkNodes,
    mapUGraphToNetwork,
    normalizeHeats,
    parseHeats,
    parseNet,
)

__all__ = [
    '__version__',
    '__author__',
    '__license__',
    'Kernel',
    'SciPYKernel',
    'PPrDiffuser',
    'NetBalancedPermuter',
    'parseHeats',
    'parseNet',
    'normalizeHeats',
    'getNetworkNodes',
    'filterLinkers',
    'findLinkerCutoff',
    'connectedSubnets',
    'mapUGraphToNetwork',
]
