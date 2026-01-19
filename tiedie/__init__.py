"""
TieDIE: Tied Diffusion for Network Discovery

A network analysis algorithm that finds subnetworks connecting genomic
perturbations to transcriptional changes in large gene interaction networks.
"""

from .kernel import Kernel
from .kernel_scipy import SciPYKernel
from .ppr import PPrDiffuser
from .permute import NetBalancedPermuter
from .util import (
    parseHeats,
    parseNet,
    normalizeHeats,
    getNetworkNodes,
    filterLinkers,
    findLinkerCutoff,
    connectedSubnets,
    mapUGraphToNetwork,
)

__version__ = '2.0.0'
__author__ = 'Evan Paull'
__email__ = 'epaull@soe.ucsc.edu'

__all__ = [
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
