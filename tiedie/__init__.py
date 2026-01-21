"""TieDIE: Tied Diffusion for Network Discovery

A network analysis algorithm that finds subnetworks connecting genomic
perturbations to transcriptional changes in large gene interaction networks.
"""

from .ppr import PPrDiffuser
from .util import (
    parseNet,
    parseHeats,
    filterLinkers,
    normalizeHeats,
    getNetworkNodes,
    connectedSubnets,
    findLinkerCutoff,
    mapUGraphToNetwork,
)
from .kernel import Kernel
from .permute import NetBalancedPermuter
from ._metadata import __author__, __license__, __version__
from .kernel_scipy import SciPYKernel


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
