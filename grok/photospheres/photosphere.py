
from astropy.io import registry
from astropy.table import Table
from .connect import PhotosphereRead, PhotosphereWrite



class Photosphere(Table):
    """A class to represent a model photosphere."""

    read = registry.UnifiedReadWriteMethod(PhotosphereRead)
    write = registry.UnifiedReadWriteMethod(PhotosphereWrite)


    @property
    def is_spherical_geometry(self):
        # TODO: we're assuming the radius units are in cm.
        radius = self.meta.get("radius", 0)
        # MARCs models give radius of 1 cm for plane-parallel models
        return radius > 1
    

    @property
    def is_plane_parallel_geometry(self):
        return not self.is_spherical_geometry