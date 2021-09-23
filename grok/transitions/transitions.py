
from astropy.io import registry
from astropy.table import (Row, Table)
from astropy import units as u
from .connect import TransitionsRead, TransitionsWrite

from .utils import (air_to_vacuum, vacuum_to_air, species_as_float, species)
    
class Transition(Row):

    def __getitem__(self, item):
        if item in self.colnames:
            return super().__getitem__(item)

        else:
            # There are some translations that we are willing to do on demand:
            # - lambda_air -> lambda
            # - lambda -> lambda_air
            # - species -> species_as_float
            # - species_as_float -> species

            functions = {
                "lambda": (air_to_vacuum, "lambda_air"),
                "lambda_air": (vacuum_to_air, "lambda"),
                "species_as_float": (species_as_float, "species"),
                "species": (species, "species_as_float"),
            }

            try:
                f, *keys = functions[item]
            except KeyError:
                raise
            else:
                return f(*(self[k] for k in keys))


class Transitions(Table):
    """A class to represent atomic and molecular line transitions."""
    
    read = registry.UnifiedReadWriteMethod(TransitionsRead)
    write = registry.UnifiedReadWriteMethod(TransitionsWrite)

    Row = Transition

    def __init__(self, *args, **kwargs):
        super(Transitions, self).__init__(*args, **kwargs)

        # Add units and descriptions for common names.
        descriptions = [
            ("lambda", u.Angstrom, "Wavelength (vacuum)"),
            ("lambda_air", u.Angstrom, "Wavelength (air)"),
            ("excitation_potential_lower", u.eV, "Lower excitation potential"),
            ("excitation_potential_upper", u.eV, "Upper excitation potential"),
            ("v_micro", u.km/u.s, "Microturbulent velocity"),
        ]

        for name, unit, description in descriptions:
            if name in self.dtype.names:
                self[name].unit = unit
                self[name].description = description
        
        return None

