
import gzip
from os import replace
import re
import numpy as np
from datetime import datetime
from astropy.io import registry
from collections import OrderedDict

from .photosphere import Photosphere
from .utils import periodic_table

def parse_meta(contents):

    patterns = [
        "(?P<file_id>^.+)\n",
        "\s+(?P<teff>[\d\.]+)\s+Teff \[(?P<teff_unit>.+)\].+yyyymmdd\=(?P<updated>\d+)\n",
        "\s+(?P<flux>[\d\.E\-\+]+)\s+Flux \[(?P<flux_unit>.+)\]\n",
        "\s+(?P<logg>[\d\.E\-\+]+)\s+Surface gravity \[(?P<surface_gravity_unit>.+)\]\n",
        "\s+(?P<microturbulence>[\d\.]+)\s+Microturbulence parameter \[(?P<microturbulence_unit>.+)\]\n",
        "\s+(?P<mass>[\d\.]+)\s+(No m|m)ass.+\n", #\s+Mass \[(?P<mass_unit>.+)\]\n"
        "\s+(?P<m_h>[\d\.\+\-]+)\s(?P<alpha_m>[\d\-\.\+]+) Metallicity.+\n",
        "\s+(?P<radius>[\d\.E\-\+]+).+(1 cm r|R)adius.+\n", #Radius \[(?P<radius_unit>.+)\] at Tau\(Rosseland\)=(?P<tau_rosseland_at_radius>[\d\.]+)\n"    
        "\s+(?P<luminosity>[\d\.E\-\+]+) Luminosity \[(?P<luminosity_unit>.+)\].*\n",
        "\s+(?P<convection_alpha>[\d\.E\-\+]+)\s+(?P<convection_nu>[\d\.E\-\+]+)\s+(?P<convection_y>[\d\.E\+]+)\s+(?P<convection_beta>[\d\.E\+]+).+\n",
        "\s+(?P<X>[\d\.E\-\+]+)\s+(?P<Y>[\d\.E\-\+]+)\s+(?P<Z>[\d\.E\-\+]+)\s+are X, Y and Z, +12C\/13C=(?P<isotope_ratio_12C_13C>\d+).*\n",
    ]
    meta = OrderedDict([])
    for pattern in patterns:
        meta.update(next(re.finditer(pattern, contents)).groupdict())
        
    pattern = "Logarithmic chemical number abundances, H always 12.00\n"
    for i, element in enumerate(periodic_table[:92], start=1):
        pattern += f"\s+(?P<log_abundance_{element}>[\d\.\-]+)"
        if (i % 10) == 0:
            pattern += "\n"

    pattern += (
        "\n"
        "\s+(?P<n_depth>\d+) Number of depth points"
    )
    if isinstance(contents, bytes):
        contents = contents.decode("utf-8")

    meta.update(next(re.finditer(pattern, contents)).groupdict())

    # Assign dtypes.
    dtype_patterns = [
        ("teff$", float),
        ("flux$", float),
        ("logg", float),
        ("microturbulence$", float),
        ("mass$", float),
        ("m_h", float),
        ("alpha_m", float),
        ("radius$", float),
        ("tau_rosseland_at_radius", float),
        ("luminosity$", float),
        ("convection_alpha", float),
        ("convection_nu", float),
        ("convection_y", float),
        ("convection_beta", float),
        ("isotope_ratio_12C_13C", float),
        ("X", float),
        ("Y", float),
        ("Z", float),
        ("log_abundance_\w+", float),
        ("n_depth", int),
        ("updated", lambda date: datetime.strptime(date, "%Y%m%d")),
    ]

    for key in meta.keys():
        for pattern, dtype in dtype_patterns:
            if re.match(pattern, key):
                meta[key] = dtype(meta[key])
                break
            
    # We want to keep the position of the surface gravity in the dict, but what we have is not
    # actually log10(surface gravity), it's just surface gravity.
    meta["logg"] = np.log10(meta["logg"])

    # 
    return meta


def loadtxt(filename, skiprows, max_rows, replace_as_nan="******"):
    
    can_opener = gzip.open if filename.lower().endswith(".gz") else open
    with can_opener(filename, "r") as fp:
        lines = fp.readlines()

    data = []
    for line in lines[skiprows:skiprows + max_rows]:
        if isinstance(line, bytes):
            line = line.decode("utf-8")
        if replace_as_nan is not None:
            line = line.replace(replace_as_nan, "NaN")
        data.append([each.group() for each in re.finditer("(-?\d{1,6}(\.\d{1,6})?(E[\-\+]\d{1,2})?|(NaN))", line)])
    return np.array(data, dtype=float)





def read_marcs(filename, structure_start=25):
    
    can_opener = gzip.open if filename.lower().endswith(".gz") else open
    with can_opener(filename, "r") as fp:
        contents = fp.read()

    if isinstance(contents, bytes):
        contents = contents.decode("utf-8")

    meta = parse_meta(contents)
    meta["filename"] = filename

    S, N = (structure_start, meta["n_depth"])
    column_locations = [
        [S, ("k", "lgTauR", "lgTau5", "Depth", "T", "Pe", "Pg", "Prad", "Pturb")],
        [S + N + 1, ("k", "lgTauR", "KappaRoss", "Density", "Mu", "Vconv", "Fconv/F", "RHOX")],
        [S + (N + 1)*2 + 1, ("k", "lgPgas", "H I", "H-", "H2", "H2+", "H2O", "OH", "CH", "CO", "CN", "C2")],
        [S + (N + 1)*3 + 1, ("k", "N2", "O2", "NO", "NH", "TiO", "2H2", "HCN", "C2H", "HS", "SiH", "C3H")],
        [S + (N + 1)*4 + 1, ("k", "C3", "CS", "SiC", "SiC2", "NS", "SiN", "SiO", "SO", "S2", "SiS", "Other")]
    ]

    data = OrderedDict([])
    for skiprows, column_names in column_locations:
        values = loadtxt(filename, skiprows=skiprows, max_rows=N)
        data.update(dict(zip(column_names, values.T)))

    column_descriptions = {
        "k": ("Depth point", "-"),
        "lgTauR": ("log(tau(Rosseland))", "-"),
        "lgTau5": ("log(tau(5000 Angstroms))", "-"),
        "Depth": ("Depth (depth=0 and tau(Rosseland)=1.0)", "cm"),
        "T": ("Temperature", "K"),
        "Pe": ("Electron pressure", "dyn/cm^2"),
        "Pg": ("Gas pressure", "dyn/cm^2"),
        "Prad": ("Radiation pressure", "dyn/cm^2"),
        "Pturb": ("Turbulence pressure", "dyn/cm^2"),
        "KappaRoss": ("Rosseland opacity", "cm^2/g"),
        "Density": ("Density", "g/cm^3"),
        "Mu": ("Mean molecular weight", "amu"),
        "Vconv": ("Convective velocity", "cm/s"),
        "RHOX": ("Mass column density", "g/cm^2"),
    }
    # Add descriptions for partial pressure columns
    index = list(data.keys()).index("H I")
    for partial_pressure_column_name in list(data.keys())[index:]:
        column_descriptions[partial_pressure_column_name] = (
            f"log({partial_pressure_column_name} pressure/[1 dyn/cm^2])",
            "-"
        )

    descriptions = { k: desc for k, (desc, unit) in column_descriptions.items() }
    units = { k: unit for k, (desc, unit) in column_descriptions.items() }

    data["k"] = data["k"].astype(int)

    photosphere = Photosphere(
        data=data,
        units=units,
        descriptions=descriptions,
        meta=meta
    )
    return photosphere


registry.register_reader("marcs", Photosphere, read_marcs)