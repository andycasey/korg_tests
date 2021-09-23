import re
import numpy as np
from astropy.io import registry
from typing import OrderedDict

from grok.transitions import Transitions

_strip_quotes = lambda s: s.strip("' ")
_strip_quotes_bytes = lambda s: _strip_quotes(s.decode("utf-8"))

def parse_references(lines):
    reference_header = "References:"
    for i, line in enumerate(lines):
        if line.strip() == reference_header:
            break
    else:
        raise ValueError(f"No references header '{reference_header}' found.")

    delimiter = ". "
    references = OrderedDict([])
    for line in lines[i+1:]:
        number, *reference_text = line.split(delimiter)
        number = int(number.lstrip())
        references[number] = delimiter.join(reference_text).rstrip()
    return references



def read_extract_stellar_long_output(path):
    """
    Read transitions that have been extracted from the Vienna Atomic Line Database
    version 3 (VALD3), using the "extract stellar" in long format.

    Documentation at https://www.astro.uu.se/valdwiki/select_output

    """

    with open(path, "r") as fp:
        lines = fp.readlines()
    
    header_pattern = "\s*(?P<lambda_start>[\d\.]+),\s*(?P<lambda_end>[\d\.]+),\s*(?P<n_selected>\d+),\s*(?P<n_processed>\d+),\s*(?P<v_micro>[\d\.]+), Wavelength region, lines selected, lines processed, Vmicro"

    meta = re.match(header_pattern, lines[0])
    if not meta:
        raise ValueError(f"Cannot match header with '{header_pattern}'")
    
    meta = meta.groupdict()
    for k in meta.keys():
        dtype = int if k.startswith("n_") else float
        meta[k] = dtype(meta[k])

    data = np.genfromtxt(
        path,
        skip_header=3,
        max_rows=meta["n_selected"],
        delimiter=",",
        dtype={
            "names": (
                "species",
                "lambda_air", 
                "excitation_potential_lower",
                "v_micro",
                "loggf", 
                "radiative_damping_constant",
                "stark_damping_constant",
                "vanderwaals_damping_constant",
                "lande_factor",
                "lande_depth",
                "reference"
            ),
            "formats": ("U10", float, float, float, float, float, float, float, float, float, "U100")
        },
        converters={
            0: _strip_quotes_bytes,
            10: _strip_quotes_bytes,
        }
    )

    meta = dict(references=parse_references(lines))

    return Transitions(data=data, meta=meta)


def read_extract_all_or_extract_element(path):
    """
    Read transitions that have been extracted from the Vienna Atomic Line Database
    version 3 (VALD3), using the "extract all" or "extract element" formats.

    Documentation at https://www.astro.uu.se/valdwiki/presformat_output
    """

    with open(path, "r") as fp:
        lines = fp.readlines()
    
    header_rows = 2
    keys_and_dtypes = (
        ("species", str),
        ("lambda_air", float),
        ("loggf", float),
        ("excitation_potential_lower", float),
        ("j_lower", float),
        ("excitation_potential_upper", float),
        ("j_upper", float),
        ("lande_factor_lower", float),
        ("lande_factor_upper", float),
        ("lande_factor_mean", float),
        ("radiative_damping_constant", float),
        ("stark_damping_constant", float),
        ("vanderwaals_damping_constant", float),
        ("line1", str),
        ("line2", str),
        ("reference", str),
    )
    data = OrderedDict([(key, []) for key, dtype in keys_and_dtypes])
    for i, line in enumerate(lines[header_rows:]):
        if line.strip() == "References:":
            break

        if i % 4 == 0:
            if line.count(",") < 10: 
                # Could be something like a footnote.
                break
            
            for (key, dtype), value in zip(keys_and_dtypes, line.split(",")):
                formatted_value = dtype(value)
                if key in ("species", ):
                    formatted_value = _strip_quotes(formatted_value)                
                data[key].append(formatted_value)

        elif i % 4 == 1:
            data["line1"].append(_strip_quotes(line))
        elif i % 4 == 2:
            data["line2"].append(_strip_quotes(line))
        elif i % 4 == 3:
            data["reference"].append(_strip_quotes(line))
    
    meta = dict(references=parse_references(lines[i+1:]))

    return Transitions(data=data, meta=meta)
    

def read_vald(path):
    """
    Read transitions exported from the Vienna Atomic Line Database
    version 3 (VALD3).
    
    Here we assume the format is not known, and we try our best.
    """

    methods = [
        read_extract_all_or_extract_element,
        read_extract_stellar_long_output,
    ]
    for method in methods:
        try:
            t = method(path)
        except:
            continue
        else:
            return t
    else:
        raise IOError(f"Unable to read VALD format. Tried methods: {methods}")
            
registry.register_reader("vald", Transitions, read_vald)
registry.register_reader("vald.stellar", Transitions, read_extract_stellar_long_output)