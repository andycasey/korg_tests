import numpy as np
import gzip
import os
import logging
from astropy.io import fits
from io import FileIO, StringIO
from shutil import copy

def read_benchmark_star(path):
    """
    Read a benchmark star spectrum that was downloaded from https://www.blancocuaresma.com/s/benchmarkstars
    
    :param path:
        The path to the FITS file.
    
    :returns:
        A three length tuple containing the wavelength (in Angstroms), flux, and flux uncertainty.
    """

    with fits.open(path) as image:
        flux = image[0].data
        sigma = image[1].data
        wl = image[0].header["CRVAL1"] \
           + np.arange(flux.size) * image[0].header["CDELT1"]

        if image[0].header["CUNIT1"] == "NM":
            wl *= 10.0
        else:
            raise NotImplementedError
        
    return (wl, flux, sigma)



def safe_open(fp_or_path):

    if isinstance(fp_or_path, (FileIO, gzip.GzipFile)):
        path = os.path.realpath(fp_or_path.name)
        contents = fp_or_path.read()
        if isinstance(contents, bytes):
            contents = contents.decode("utf-8")

    else:
        can_opener = gzip.open if fp_or_path.lower().endswith(".gz") else open
        with can_opener(fp_or_path, "r") as fp:
            contents = fp.read()

        if isinstance(contents, bytes):
            contents = contents.decode("utf-8")

        path = fp_or_path

    content_stream = StringIO(contents)

    return (path, contents, content_stream)


def copy_or_write(input_item, destination_path, **kwargs):
    """
    Copy or write an input object (e.g., photosphere, transitions) to the given
    destination path.

    :param input_item:
        The relevant input object. This could be a photosphere, a set of transitions,
        or a path to one of those items.
    
    :param destination_path:
        The path where the input object should be written.
    """

    os.makedirs(os.path.dirname(destination_path), exist_ok=True)

    if isinstance(input_item, (str, bytes)) and os.path.exists(input_item):
        copy(input_item, destination_path)

    else:
        input_item.write(destination_path, **kwargs)
    


class CustomFormatter(logging.Formatter):

    grey = "\x1b[37;20m"
    white = "\x1b[38;20m"
    yellow = "\x1b[33;20m"
    red = "\x1b[31;20m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"
    format = "%(asctime)s [%(levelname)s] %(message)s (%(filename)s:%(lineno)d)"

    FORMATS = {
        logging.DEBUG: grey + format + reset,
        logging.INFO: white + format + reset,
        logging.WARNING: yellow + format + reset,
        logging.ERROR: red + format + reset,
        logging.CRITICAL: bold_red + format + reset
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)


def get_logger(level=logging.DEBUG):
    """
    Get a logger for this code base.
    """
    logger = logging.getLogger("grok")
    logger.setLevel(level)

    ch = logging.StreamHandler()
    ch.setLevel(level)

    ch.setFormatter(CustomFormatter())
    logger.addHandler(ch)
    return logger


periodic_table = """H                                                  He
                    Li Be                               B  C  N  O  F  Ne
                    Na Mg                               Al Si P  S  Cl Ar
                    K  Ca Sc Ti V  Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr
                    Rb Sr Y  Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I  Xe
                    Cs Ba Lu Hf Ta W  Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn
                    Fr Ra Lr Rf Db Sg Bh Hs Mt Ds Rg Cn UUt"""

lanthanoids    =   "La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb"
actinoids      =   "Ac Th Pa U  Np Pu Am Cm Bk Cf Es Fm Md No"

periodic_table = periodic_table.replace(" Ba ", " Ba " + lanthanoids + " ") \
    .replace(" Ra ", " Ra " + actinoids + " ").split()

    