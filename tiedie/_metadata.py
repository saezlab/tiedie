"""Package metadata (version, authors, etc)."""

__all__ = ['__version__', '__author__', '__license__']

import importlib.metadata

_FALLBACK_VERSION = '2.0.1'

try:
    __version__ = importlib.metadata.version('tiedie')
except importlib.metadata.PackageNotFoundError:
    # Package not installed (e.g. running from source checkout)
    __version__ = _FALLBACK_VERSION

__author__ = 'Evan O. Paull'
__license__ = 'GPL-3.0-or-later'
