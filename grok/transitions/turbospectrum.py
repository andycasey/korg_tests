import re
import numpy as np
from astropy.io import registry

from grok.transitions import Transitions
from grok.utils import safe_open

_header_pattern = "'\s*(?P<turbospectrum_species_as_float>[\d\.]+)\s*'\s*(?P<ionisation>\d+)\s+(?P<num>\d+)"
_line_pattern = "\s*(?P<lambda>[\d\.]+)\s*(?P<E_lower>[\-\d\.]+)\s*(?P<log_gf>[\-\d\.]+)\s*(?P<gamma_vdW_or_ABO>[\-\d\.]+)\s*(?P<g_upper>[\d\.]+)\s*(?P<gamma_rad>[\d\.E\+\-]+)\s'(?P<lower_orbital_type>\w)' '(?P<upper_orbital_type>\w)'\s+(?P<equivalent_width>[\d\.]+)\s+(?P<equivalent_width_error>[\d\.]+)\s'(?P<comment>.+)'"

# Amazingly, there does not seem to be a python string formatting that does the following:
_format_log_gf = lambda log_gf: "{0:< #6.3f}".format(log_gf)[:6]
# You might think that "{0:< #6.4g}" would work, but try it for -0.002 :head_exploding:
_line_template = "{line[lambda]:10.3f} {line[E_lower]:6.3f} {formatted_log_gf:s} {line[gamma_vdW_or_ABO]:8.3f} {line[g_upper]:6.1f} {line[gamma_rad]:9.2E} '{line[lower_orbital_type]:s}' '{line[upper_orbital_type]:s}' {line[equivalent_width]:5.1f} {line[equivalent_width_error]:6.1f} '{line[comment]}'"


def read_transitions(path):
    """
    Read transitions formatted for Turbospectrum.
    
    :param path:
        The path where the transitions are written.
    """

    path, content, content_stream = safe_open(path)
    lines = content.split("\n")

    keys_as_floats = (
        "turbospectrum_species_as_float", "lambda", "E_lower", "log_gf",
        "gamma_vdW_or_ABO", "g_upper", "gamma_rad",
        "equivalent_width", "equivalent_width_error"
    )

    # https://github.com/bertrandplez/Turbospectrum2019/blob/master/Utilities/vald3line-BPz-freeformat.f#448
    # g_upper = j_upper * 2 + 1

    i, data = (0, [])
    while i < len(lines):
        if len(lines[i].strip()) == 0:
            i += 1
            continue

        common = re.match(_header_pattern, lines[i]).groupdict()
        common["ionisation"] = int(common["ionisation"])
        num = int(common.pop("num"))

        common["representation"] = lines[i + 1][1:-1].strip()

        for j, line in enumerate(lines[i+2:i+2+num], start=i+2):
            transition = re.match(_line_pattern, line).groupdict()
            
            row = { **common, **transition }

            # Format things.
            for key in keys_as_floats:
                row[key] = float(row[key])
            data.append(row)

            raise NotImplementedError("parse vdW or ABO")

        i += 2 + num
    
    return Transitions(data=data)


def write_transitions(transitions, path):
    """
    Write transitions to disk in a format that Turbospectrum accepts.
    
    :param transitions:
        The atomic and molecular transitions.
    
    :param path:
        The path to store the transitions on disk.
    """

    lines = []

    for group in transitions.group_by(("turbospectrum_species_as_float", "representation")).groups:

        group.sort("lambda")

        # Add header information.
        is_molecule = group[0]["turbospectrum_species_as_float"] > 100
        if is_molecule:
            lines.extend([
                f"'{group[0]['turbospectrum_species_as_float']: >11.6f}         ' {group[0]['ionisation']: >4.0f} {len(group): >9.0f}",
                f"'{group[0]['representation']:7s}'"
            ])
        else:
            repr = f"'{group[0]['turbospectrum_species_as_float']: >8.3f}            '"
            lines.extend([
                f"{repr} {group[0]['ionisation']: >4.0f} {len(group): >9.0f}",
                f"'{group[0]['representation']:7s}'"
            ])
        
        for line in group:
            lines.append(
                _line_template.format(
                    line=line,
                    formatted_log_gf=_format_log_gf(line["log_gf"])
                )
            )

    
    with open(path, "w") as fp:
        fp.write("\n".join(lines))

    return None



def identify_turbospectrum(origin, *args, **kwargs):
    """
    Identify a line list being in Turbospectrum format.

    There's no accepted nomenclature for line list filenames, so we will check
    the header of the file.
    """
    first_line = args[1].readline()
    if isinstance(first_line, bytes):
        first_line = first_line.decode("utf-8")

    # Reset pointer.
    args[1].seek(0)

    return (re.match(_header_pattern, first_line) is not None)

registry.register_reader("turbospectrum", Transitions, read_transitions)
registry.register_writer("turbospectrum", Transitions, write_transitions)
registry.register_identifier("turbospectrum", Transitions, identify_turbospectrum)