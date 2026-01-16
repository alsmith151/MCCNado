from .mccnado import *
from .version import __version__


__doc__ = mccnado.__doc__
if hasattr(mccnado, "__all__"):
    __all__ = mccnado.__all__