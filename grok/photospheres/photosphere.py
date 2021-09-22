
from astropy.io import registry
from astropy.table import Table
from .connect import PhotosphereRead, PhotosphereWrite



class Photosphere(Table):
    """A class to represent a model photosphere."""

    read = registry.UnifiedReadWriteMethod(PhotosphereRead)
    write = registry.UnifiedReadWriteMethod(PhotosphereWrite)
