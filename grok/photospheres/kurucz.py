
import re
import gzip
import numpy as np
from collections import OrderedDict
from astropy.io import registry

from .photosphere import Photosphere
from grok.utils import periodic_table, safe_open


def parse_meta(contents):
    # Documented at https://marcs.astro.uu.se/krz_format.html

    # Check if spherical or not. Could also do this by checking for model type:
    # MODEL TYPE = 3 (spherical)
    # MODEL TYPE = 0 (plane parallel)
    spherical_pattern = ".+SPHERICAL, RADIUS=\s(?P<radius>[\d\.E\+]+)\s+(?P<radius_unit>.+)"
    is_spherical = any(re.findall(spherical_pattern, contents))
    
    pattern = (
        "(?P<header>^.+)\n"# V[turb|TURB]\s+(?P<microturbulence>[\d\.]+) km/s.*\n"
        "T EFF=\s?(?P<teff>[\d\.]+)\s+GRAV=\s?(?P<logg>[\-\d\.]+)\s+MODEL TYPE= (?P<model_type>\d) WLSTD= (?P<standard_wavelength>\d+)"
    )
    if is_spherical:
        pattern += spherical_pattern + ".+\n"
    else:
        pattern += ".+\n"

    pattern += "\s*(?P<absorption_by_H>[01]) (?P<absorption_by_H2_plus>[01]) (?P<absorption_by_H_minus>[01]) (?P<rayleigh_scattering_by_H>[01]) (?P<absorption_by_He>[01]) (?P<absorption_by_He_plus>[01]) (?P<absorption_by_He_minus>[01]) (?P<rayleigh_scattering_by_He>[01]) (?P<absorption_by_metals_for_cool_stars_Si_Mg_Al_C_Fe>[01]) (?P<absorption_by_metals_for_intermediate_stars_N_O_Si_plus_Mg_plus_Ca_plus>[01]) (?P<absorption_by_metals_for_hot_stars_C_N_O_Ne_Mg_Si_S_Fe>[01]) (?P<thomson_scattering_by_electrons>[01]) (?P<rayleigh_scattering_by_h2>[01]).+OPACITY SWITCHES\n"

    #normalized abundances for the first 99 elements in the periodic table: N_atom/N_all where N are number densities. The sum of all abundances is 1. Positive numbers are fractions while negative numbers are log10 of fractions. The last number in line 13 is the number of depth layers in the model.
    for i, element in enumerate(periodic_table[:99], start=1):
        pattern += f"\s*(?P<log10_normalized_abundance_{element}>[\-\d\.]+) "
        if i % 10 == 0:
            pattern = pattern.rstrip() + "\n"

    pattern += "\s*(?P<n_depth>\d+)"
    
    meta = OrderedDict(re.match(pattern, contents).groupdict())
    if not is_spherical:
        meta.update(radius=0, radius_unit="cm")
    
    # Assign dtypes.
    dtype_patterns = [
        ("teff$", float),
        ("flux$", float),
        ("logg$", float),
        ("radius$", float),
        ("standard_wavelength", float),
        ("log10_normalized_abundance_\w+", float),
        ("n_depth", int),
        ("absorption_by", lambda v: bool(int(v))),
        ("scattering_by", lambda v: bool(int(v))),
    ]
    for key in meta.keys():
        for pattern, dtype in dtype_patterns:
            if re.findall(pattern, key):
                meta[key] = dtype(meta[key])
                break    

    # Correct the normalized abundance columns so they are all log10 base.
    # Numbers that are positive are not yet log10
    for element in periodic_table[:99]:
        key = f"log10_normalized_abundance_{element}"
        if meta[key] > 0:
            meta[key] = np.log10(meta[key])

    # Add overall metallicity keyword, taking Fe as representative of metallicity.
    # Remember that the normalized_abundance keywords are number densities, so to convert to the "12 scale" we need:
    #   log_10(X) = log_10(N_X/N_H) + 12
    # and for [X/H] we need:
    #   [X/H] = log_10(X) - log_10(X_Solar)

    # TODO: Assuming 7.45 Fe solar composition. Check this.
    meta["m_h"] = meta["log10_normalized_abundance_Fe"] - meta["log10_normalized_abundance_H"] + 12 - 7.45
    
    # TODO: Should we just calculate abundances for ease?
    meta["grid_keywords"] = ("teff", "logg", "m_h")
    return meta



def read_kurucz(fp_or_path, structure_start=14):

    filename, contents, content_stream = safe_open(fp_or_path)

    meta = parse_meta(contents)

    if meta["radius"] > 0:
        # Spherical models have six columns.
        usecols = (0, 1, 2, 3, 4, 5)
        column_locations = [
            (structure_start, ("tau", "T", "XNE", "numdens_other", "Density", "Depth"))
        ]      
    else:
        # Plane-parallel models have four columns
        usecols = (0, 1, 2, 3, 4)
        column_locations = [
            (structure_start, ("RHOX", "T", "XNE", "numdens_other", "Density", "Depth"))
        ]      
    
    data = {}
    for skiprows, column_names in column_locations:
        values = np.loadtxt(
            content_stream, 
            skiprows=skiprows, 
            max_rows=meta["n_depth"],
            delimiter=",",
            usecols=usecols
        )
        data.update(dict(zip(column_names, values.T)))

    column_descriptions = {
        "RHOX": ("Mass column density", "g/cm^2"),
        "tau": ("Optical depth at standard wavelength", "-"),
        "T": ("Temperature", "K"),
        "XNE": ("Number density of free electrons", "1/cm^3"),
        "numdens_other": ("Number density of all other particles (except electrons)", "1/cm^3"),
        "Density": ("Density", "g/cm^3"),
        "Depth": ("Height relative to the reference radius given in the metadata", "cm")
    }
    descriptions = { k: desc for k, (desc, unit) in column_descriptions.items() if k in data }
    units = { k: unit for k, (desc, unit) in column_descriptions.items() if k in data  }

    return Photosphere(
        data=data,
        units=units,
        descriptions=descriptions,
        meta=meta
    )

def identify_kurucz(origin, *args, **kwargs):
    return (isinstance(args[0], str) and \
            args[0].lower().endswith((".krz", ".krz.gz")))

registry.register_reader("kurucz", Photosphere, read_kurucz)
registry.register_identifier("kurucz", Photosphere, identify_kurucz)