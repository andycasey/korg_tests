import gzip
import numpy as np
from textwrap import dedent
from astropy.io import registry

from grok.photospheres import Photosphere


def write_photosphere_for_turbospectrum(photosphere, path):
    """
    Write a photosphere in an ASCII format that Turbospectrum will recognise.

    :param photosphere:
        The photosphere, an instance of `photospheres.Photosphere`.
    
    :param path:
        The path to store the photosphere.    
    """

    # TODO: A hack because we are testing things without interpolation here.
    if photosphere.meta.get("read_format", None) == "marcs":
        print(f"WARNING: Copying MARCS file directly. I hope you're not interpolating atmospheres!")

        original_path = photosphere.meta["filename"]
        can_opener = gzip.open if original_path.lower().endswith(".gz") else open
        with can_opener(original_path, "r") as fp:
            content = fp.read()
        with open(path, "w") as fp:
            fp.write(content.decode("utf-8"))
        return None

    # TODO: Turbospectrum describes this format as 'KURUCZ', but it looks like an ATLAS-style format to me.
    # NOTE: Turbospectrum does not read the abundance information from the photosphere. It reads it from the control file.
    output = (
        f"'KURUCZ' {photosphere.meta['n_depth']} {photosphere.meta.get('standard_wavelength', 5000):.0f} {photosphere.meta['logg']:.2f} 0 0.\n"
    )

    # See https://github.com/bertrandplez/Turbospectrum2019/blob/master/source-v19.1/babsma.f#L727
    line_format = " {line[RHOX]:.8e} {line[T]: >8.1f} {line[P]:.3e} {line[XNE]:.3e} {line[ABROSS]:.3e} {line[ACCRAD]:.3e} {line[VTURB]:.3e} {line[FLXCNV]:.3e}\n"

    for i, line in enumerate(photosphere, start=1):
        output += line_format.format(i=i, line=line)

    with open(path, "w") as fp:
        fp.write(output)
    return None



registry.register_writer("turbospectrum", Photosphere, write_photosphere_for_turbospectrum)