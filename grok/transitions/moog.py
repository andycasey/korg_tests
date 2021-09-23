
import numpy as np
from astropy.io import registry
from collections import OrderedDict

from grok.transitions import Transitions


def read_moog(path):
    """
    Read a MOOG-formatted list of transitions.

    :param path:
        The path of the line list.
    """

    with open(path, "r") as fp:
        lines = fp.readlines()
    
    data = OrderedDict([
        ("lambda", []),
        ("species_as_float", []),
        ("excitation_potential_lower", []),
        ("loggf", []),
        ("vanderwaals_damping_constant", []),
        ("dissociation_energy", []),
        ("comment", [])
    ])
    
    for i, line in enumerate(lines, start=1):
        s = line.split()
        if len(s) < 1:
            continue

        try:
            lambda_, species, ep, loggf = map(float, s[:4])
        except:
            if i == 1:
                continue
            else:
                raise IOError(f"Invalid contents on line {i}")
            
        else:
            damping, dissoc, comment = (np.nan, np.nan, "")
            if len(s) > 4:
                try:
                    damping, dissoc, *_ = line[40:60].strip().split()
                    damping, dissoc = (float(damping), float(dissoc))
                except: 
                    None
                comment = line[50:].strip()
                
        data["lambda"].append(lambda_)
        data["species_as_float"].append(species)
        data["excitation_potential_lower"].append(ep)
        data["loggf"].append(loggf)
        data["vanderwaals_damping_constant"].append(damping)
        data["dissociation_energy"].append(dissoc)
        data["comment"].append(comment)

    return Transitions(data=data)


def write_moog(transitions, path):
    """
    Write a list of atomic and molecular transitions to disk in a format that
    is friendly to MOOG.
    
    :param transitions:
        The atomic and molecular transitions to write.
    
    :param path:
        The path to store the transitions.
    """

    fmt = "{:10.3f}{:10.5f}{:10.3f}{:10.3f}{}{}{}{}"
    space = " "*10
    with open(path, "w") as f:
        f.write("\n")

        for line in transitions:
            C6 = space if np.ma.is_masked(line['vanderwaals_damping_constant']) or np.isnan(line['vanderwaals_damping_constant']) else "{:10.3f}".format(line['vanderwaals_damping_constant'])
            try:
                D0 = line["dissociation_energy"]
            except:
                D0 = np.nan
            D0 = space if np.ma.is_masked(D0) or np.isnan(D0) else "{:10.3}".format(D0)
            try:
                comment = '' if np.ma.is_masked(line['comment']) else line['comment']
            except:
                comment = ""
            f.write(
                fmt.format(
                    line['lambda'],
                    line['species_as_float'],
                    line['excitation_potential_lower'],
                    line['loggf'],
                    C6,
                    D0,
                    space,
                    comment
                ) + "\n"
            )

    return None

registry.register_reader("moog", Transitions, read_moog)
registry.register_writer("moog", Transitions, write_moog)