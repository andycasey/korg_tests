
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
        ("wavelength", []),
        ("species", []),
        ("excitation_potential", []),
        ("loggf", []),
        ("vanderwaals_damping", []),
        ("dissociation_energy", []),
        ("comment", [])
    ])
    
    for i, line in enumerate(lines, start=1):
        s = line.split()
        if len(s) < 1:
            continue

        try:
            wavelength, species, ep, loggf = map(float, s[:4])
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
                
        data["wavelength"].append(wavelength)
        data["species"].append(species)
        data["excitation_potential"].append(ep)
        data["loggf"].append(loggf)
        data["vanderwaals_damping"].append(damping)
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
            C6 = space if np.ma.is_masked(line['vanderwaals_damping']) or np.isnan(line['vanderwaals_damping']) else "{:10.3f}".format(line['vanderwaals_damping'])
            D0 = space if np.ma.is_masked(line['dissociation_energy']) or np.isnan(line['dissociation_energy']) else "{:10.3}".format(line['dissociation_energy'])
            comment = '' if np.ma.is_masked(line['comment']) else line['comment']
            f.write(fmt.format(line['wavelength'],line['species'],line['excitation_potential'],line['loggf'],C6,D0,space,line['comment'])+"\n")

    return None

registry.register_reader("moog", Transitions, read_moog)
registry.register_writer("moog", Transitions, write_moog)
