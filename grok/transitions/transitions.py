
from astropy.io import registry
from astropy.table import Table
from .connect import TransitionsRead, TransitionsWrite


class Transitions(Table):
    """A class to represent atomic and molecular line transitions."""
    
    read = registry.UnifiedReadWriteMethod(TransitionsRead)
    write = registry.UnifiedReadWriteMethod(TransitionsWrite)