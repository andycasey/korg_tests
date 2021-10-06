import re
import numpy as np
from astropy.io import registry
from astropy import units as u

from grok.transitions import (Species, Transition, Transitions)
from grok.utils import safe_open

_header_pattern = "'\s*(?P<turbospectrum_species_as_float>[\d\.]+)\s*'\s*(?P<ionisation>\d+)\s+(?P<num>\d+)"
_line_pattern = "\s*(?P<lambda_air>[\d\.]+)\s*(?P<E_lower>[\-\d\.]+)\s*(?P<log_gf>[\-\d\.]+)\s*(?P<vdW>[\-\d\.]+)\s*(?P<g_upper>[\d\.]+)\s*(?P<gamma_rad>[\d\.E\+\-]+)\s'(?P<lower_orbital_type>\w)' '(?P<upper_orbital_type>\w)'\s+(?P<equivalent_width>[\d\.]+)\s+(?P<equivalent_width_error>[\d\.]+)\s'(?P<comment>.+)'"

# Amazingly, there does not seem to be a python string formatting that does the following:
_format_log_gf = lambda log_gf: "{0:< #6.3f}".format(log_gf)[:6]
# You might think that "{0:< #6.4g}" would work, but try it for -0.002 :head_exploding:
_line_template = "{line.lambda_air.value:10.3f} {line.E_lower.value:6.3f} {formatted_log_gf:s} {line.vdW_compact:8.3f} {line.g_upper:6.1f} {line.gamma_rad.value:9.2E} '{line.lower_orbital_type:s}' '{line.upper_orbital_type:s}' {line.equivalent_width:5.1f} {line.equivalent_width_error:6.1f} '{line.comment}'"


def should_skip(transition):
    """
    Returns a boolean flag whether a transition should be excluded (True) or included (False)
    from Turbospectrum.

    The logic here follows that from the VALD export script by Bertrand Plez.
    """
    if not (2 > transition.species >= 0) \
    or transition.E_lower >= (15 * u.eV):
        return True

    # TODO: autoionization lines.
    # why Z < 3 for any? that would exclude MgH, etc.
    
    return False


def update_missing_transition_data(transition):
    """
    Validate (and update) transition data. This is done by the script that converts VALD line lists
    to Turbospectrum formats. Specifically it checks the van der Waals broadening constant and
    gamma radiative. If it doesn't like the values, it changes them.

    This returns a copy of the transition, with the updated data.
    """
    t = transition.copy()
    if np.isclose(t.vdW, 0, atol=1-e5):
        t.vdW = lookup_approximate_vdW(t)
    
    t.gamma_rad = parse_gamma_rad(t)
    return t


def parse_gamma_rad(transition):
    #https://github.com/bertrandplez/Turbospectrum2019/blob/master/Utilities/vald3line-BPz-freeformat.f#L420-428

    if np.isclose(transition.gamma_rad.value, 0, atol=1e-5)
        return 1e5 * (1/u.s)
    elif transition.gamma_rad.value > 3:
        return 10**transition.gamma_rad.value * (1/u.s)


def lookup_approximate_vdW(transition):

    default_value = 2.5 # Mackle et al. 1975 A&A 38, 239

    neutral_damping = { 
        "Na": 2.0, # Holweger 1971 A&A 10, 128
        "Si": 1.3, # Holweger 1973 A&A 26, 275
        "Ca": 1.8, # O'Neill & Smith 1980 A&A 81, 100
        "Fe": 1.4  # Simmons & Blackwell 1982 A&A 112, 209 and Magain & Zhao 1996 A&A 305, 245
    }
    ionized_damping = {
        "Ca": 1.4, # from fit of H&K in the HM model to the fts intensity spectrum
        "Sr": 1.8, # from fit of Sr II 4077.724 in the HM model to the fts intensity spectrum
        "Ba": 3.0 # Holweger & Muller 1974 Solar Physics 39, 19
    }
    is_molecule = (len([Z for Z in transition.species.Zs if Z > 0]) > 1)
    if not is_molecule:
        if transition.species.charge == 0:
            try:
                return neutral_damping[transition.species.atoms[0]]
            except KeyError:
                return default_value
        elif transition.species.charge == 1:
            try:
                return ionized_damping[transition.species.atoms[0]]
            except KeyError:
                return default_value
    else:
        return default_value


def read_transitions(path):
    """
    Read transitions formatted for Turbospectrum.
    
    :param path:
        The path where the transitions are written.
    """

    path, content, content_stream = safe_open(path)
    lines = content.split("\n")

    keys_as_floats = (
        "lambda_air", "E_lower", "log_gf",
        "vdW", "g_upper", "gamma_rad",
        "equivalent_width", "equivalent_width_error"
    )    
    i, transitions = (0, [])
    while i < len(lines):
        if len(lines[i].strip()) == 0:
            i += 1
            continue

        common = re.match(_header_pattern, lines[i]).groupdict()
        num = int(common.pop("num"))

        # Note that we are interpreting the ionisation state from the species representation,
        # and ignoring the 'ionisation' in 'common'.
        common["species"] = Species(common["turbospectrum_species_as_float"])
        # Update the charge, using the representation like "Th II", since Turbospectrum does
        # not encode the ionisation in this species representation.
        common["species"].charge = Species(lines[i + 1][1:-1]).charge
        
        for j, line in enumerate(lines[i+2:i+2+num], start=i+2):
            transition = re.match(_line_pattern, line).groupdict()
            
            row = { **common, **transition }
            
            # Format things.
            for key in keys_as_floats:
                row[key] = float(row[key])

            # Add units.
            row["lambda_air"] *= u.Angstrom
            row["E_lower"] *= u.eV
            row["gamma_rad"] *= (1/u.s)

            # https://github.com/bertrandplez/Turbospectrum2019/blob/master/Utilities/vald3line-BPz-freeformat.f#448
            # g_upper = j_upper * 2 + 1
            row["j_upper"] = 0.5 * (row["g_upper"] - 1)
            transitions.append(Transition(
                lambda_vacuum=None,
                **row
            ))
        i += 2 + num
    
    return Transitions(transitions)
    

def write_transitions(transitions, path):
    """
    Write transitions to disk in a format that Turbospectrum accepts.
    
    :param transitions:
        The atomic and molecular transitions.
    
    :param path:
        The path to store the transitions on disk.
    """

    # Sort the right way.
    group_names = list(map(str, sorted(set([float(line.species.compact) for line in transitions]))))
    
    lines = []
    for i, group_name in enumerate(group_names):

        group = sorted(
            filter(lambda t: str(float(t.species.compact)) == group_name, transitions),
            key=lambda t: t.lambda_vacuum
        )

        # Add header information.
        species = group[0].species

        formula = "".join([f"{Z:0>2.0f}" for Z in species.Zs if Z > 0])
        isotope = "".join([f"{I:0>3.0f}" for Z, I in zip(species.Zs, species.isotopes) if Z > 0])

        # Left-pad and strip formula.
        formula = formula.lstrip("0")
        formula = f"{formula: >4}"

        compact = f"{formula}.{isotope}"

        lines.extend([
            f"'{compact: <20}' {species.charge + 1: >4.0f} {len(group): >9.0f}",
            f"'{str(species):7s}'"
        ])

        for line in group:

            # Turbospectrum makes some approximations and lookups if the input values are not what
            # they want.
            updated_line = update_missing_transition_data(line)
            lines.append(
                _line_template.format(
                    line=updated_line,
                    formatted_log_gf=_format_log_gf(line.log_gf)
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
    if args[1] is None:
        return False
    
    first_line = args[1].readline()
    if isinstance(first_line, bytes):
        first_line = first_line.decode("utf-8")

    # Reset pointer.
    args[1].seek(0)

    return (re.match(_header_pattern, first_line) is not None)

registry.register_reader("turbospectrum", Transitions, read_transitions)
registry.register_writer("turbospectrum", Transitions, write_transitions)
registry.register_identifier("turbospectrum", Transitions, identify_turbospectrum)