import gzip
import os
from io import FileIO, StringIO
from shutil import copy


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

    if isinstance(input_item, (str, bytes)):    
        copy(input_item, destination_path)

    else:
        input_item.write(destination_path, **kwargs)
        


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