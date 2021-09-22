
import numpy as np
from textwrap import dedent
from astropy.io import registry

from grok.photospheres import Photosphere


def write_photosphere_for_moog(photosphere, path, format=None):
    """
    Write a photosphere to disk in a format that is known to MOOG.
    
    :param photosphere:
        The photosphere, an instance of `photospheres.Photosphere`.
    
    :param path:
        The path to store the photosphere.
    
    :param format:
        The format to write for the photosphere. This can be one of:

        - WEBMARCS
        - KURUCZ
    """

    if format is None:
        # Try to intelligently figure it out.
        if "XNE" in photosphere.dtype.names:
            return write_photosphere_for_moog(photosphere, path, format="KURUCZ")
        else:
            return write_photosphere_for_moog(photosphere, path, format="WEBMARCS")


    available_formats = {
        #"NEWMARCS": 
        # WEBMARCs wants Ne, but if we give Pe then MOOG will convert it.
        # Recall: Ne = (Pe [dyn/cm^2] / T [K]) / (k_Boltzmann [erg / K])

        "WEBMARCS": " {i:>3.0f} {i:>3.0f} {line[lgTau5]:10.3e} {i:>3.0f} {line[T]:10.3e} {line[Pe]:10.3e} {line[Pg]:10.3e}",
        #"WEB2MARC":
        "KURUCZ": " {line[RHOX]:.8e} {line[T]:10.3e}{line[P]:10.3e}{line[XNE]:10.3e}{line[ABROSS]:10.3e}",
        #"NEXTGEN":
        #"BEGN":
        #"KURTYPE":
        #"KUR-PADOVA":
        #"GENERIC": 
    }
    format = str(format).strip().upper()
    try:
        line_format = available_formats[format]
    except KeyError:
        raise ValueError(f"Format '{format}' not one of the known formats: {', '.join(list(available_formats.keys()))}")
    

    output = dedent(f"""
        {format}
            PHOTOSPHERE {photosphere.meta['teff']:.0f} / {photosphere.meta['logg']:.3f} / {photosphere.meta['m_h']:.3f} / {photosphere.meta.get('microturbulence', 0):.3f}
        NTAU       {photosphere.meta['n_depth']}
        """
    ).lstrip()

    # Add standard wavelength
    if format == "WEBMARCS":
        output += "5000.0\n"

    for i, line in enumerate(photosphere, start=1):
        output += line_format.format(i=i, line=line) + "\n"
        
    output += f"        {photosphere.meta['microturbulence']:.3f}\n"
    output += f"NATOMS        0     {photosphere.meta['m_h']:.3f}\n"
    output += "NMOL          0\n"
    # MOOG11 fails to read if you don't add an extra line
    output += "\n"

    with open(path, "w") as fp:
        fp.write(output)
    
    return None


def parse_single_spectrum(lines):
    """
    Parse the header, dispersion and depth information for a single spectrum 
    that has been output to a summary synthesis file.

    :param lines:
        A list of string lines partitioned from a MOOG summary output file.

    :returns:
        A three-length tuple containing the (1) header rows, (2) an array of
        dispersion values, and (3) an array of intensity values.
    """

    # Skip over the header information
    for i, line in enumerate(lines):
        if line.startswith("MODEL:"): break

    else:
        raise ValueError("could not find model information for spectrum")

    # Get the dispersion information.
    start, end, delta, _ = np.array(lines[i + 1].strip().split(), dtype=float)

    # If MOOG doesn't have opacity contributions at a given wavelength, it will
    # just spit out ****** entries.
    #_pre_format = lambda l: l.strip().replace("*******", " 0.0000").replace("-0.0000"," 0.0000").split()
    def _pre_format(l):
        l = l.replace("*******", " 0.0000").rstrip()
        assert len(l) % 7 == 0, len(l)
        return [l[7*i:7*(i+1)] for i in range(len(l)//7)]

    depths = np.array(
        sum([_pre_format(line) for line in lines[i + 2:]], []),
        dtype=float)

    dispersion = np.arange(start, end + delta, delta)[:depths.size]
    intensity = 1.0 - depths
    
    # Parse the headers into metadata
    meta = {
        "raw": lines[:i + 2]
    }
    return (dispersion, intensity, meta)


def parse_summary_synth_output(summary_out_path):
    """
    Parse the summary output from a MOOG synth operation.

    :param summary_out:
        The path of the summary output file produced by MOOG.
    """

    with open(summary_out_path, "r") as fp:
        summary = fp.readlines()

    # There could be multiple spectra in this output file, each separated by 
    # headers. Scan to find the headers, then extract the information.
    partition_indices = [i for i, line in enumerate(summary) \
        if line.lower().startswith("all abundances not listed below differ")] \
        + [None]

    spectra = []
    for i, start in enumerate(partition_indices[:-1]):
        end = partition_indices[i + 1]
        spectra.append(parse_single_spectrum(summary[start:end]))

    return spectra


def write_marcs_photosphere_for_moog(photosphere, path):
    return write_photosphere_for_moog(photosphere, path, format="WEBMARCS")

def write_kurucz_photosphere_for_moog(photosphere, path):
    return write_photosphere_for_moog(photosphere, path, format="KURUCZ")

registry.register_writer("moog", Photosphere, write_photosphere_for_moog)
registry.register_writer("moog.marcs", Photosphere, write_marcs_photosphere_for_moog)
registry.register_writer("moog.kurucz", Photosphere, write_kurucz_photosphere_for_moog)
