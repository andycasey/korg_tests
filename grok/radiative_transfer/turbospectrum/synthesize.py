import os
import numpy as np
import subprocess
from collections import OrderedDict
from pkg_resources import resource_stream

from grok.radiative_transfer.utils import get_default_lambdas
from grok.utils import copy_or_write

def turbospectrum_bsyn(
        photosphere,
        transitions,
        lambdas=None,
        abundances=None,
        isotopes=None,
        opacities=None,
        verbose=False,
        dir=None,
        **kwargs
    ):
    """
    Synthesize a stellar spectrum using Turbospectrum's `bsyn` routine.
    
    :param photosphere:
        The model photosphere.
    
    :param transitions:
        The atomic and molecular line transitions.
        
    :param lambdas: [optional]
        A three length tuple containing the start wavelength, end wavelength, and the
        wavelength step size. If `None` is given then this will default to nearly the 
        range of the transition list, with a step size of 0.01 Angstroms.

    # TODO: Docs for abundances, isotopes.

    :param opacities: [optional]
        A path where the calculated opacities are stored. If None is given, then
        these will be calculated by Turbospectrum.

    :param verbose: [optional]
        Provide verbose outputs.

    :param dir: [optional]
        The directory to execute Turbospectrum from.
    """    

    if lambdas is not None:
        lambda_min, lambda_max, lambda_delta = lambdas
    else:
        lambda_min, lambda_max, lambda_delta = get_default_lambdas(transitions)

    _path = lambda basename: os.path.join(dir or "", basename)


    transition_basename_format = "transitions_{i}"
    modelinput_basename = "model"
    modelopac_basename = "opac"
    result_basename = "result"

    if opacities is not None:
        copy_or_write(
            opacities,
            _path(modelopac_basename)
        )        

    # Turbospectrum allows for transitions to be written into multiple different files.
    # TODO: Need to find out if that's ever a legitamite use case. Otherwise let's merge
    #       all to one file.
    T = 1
    copy_or_write(
        transitions,
        _path(transition_basename_format.format(i=0)), 
        format=kwargs.get("transition_format", "turbospectrum")
    )

    # Write photosphere.
    copy_or_write(
        photosphere,
        _path(modelinput_basename),
        format=kwargs.get("photosphere_format", "turbospectrum")
    )

    kwds = dict(
        lambda_min=lambda_min,
        lambda_max=lambda_max,
        lambda_delta=lambda_delta,
        marcs_file_flag_str=".false.",
        metallicity=0,          # TODO: parse from abundances
        alpha_fe=0,             # TODO: parse from abundances
        helium_abundance=0,     # TODO: parse from abundances
        r_process_abundance=0,  # TODO: parse from abundances
        s_process_abundance=0,  # TODO: parse from abundances
        modelinput_basename=modelinput_basename,
        modelopac_basename=modelopac_basename,
        result_basename=result_basename,
        microturbulence=2.0, # TODO
        N_transition_paths=T,
        transition_paths_formatted="\n".join(
            [transition_basename_format.format(i=i) for i in range(T)]
        )
    )

    if abundances is not None or isotopes is not None:
        raise NotImplementedError        
    


    # Symbolically link the turbospectrum data folder.
    # TODO: UGH, this is a pain in the ass.
    os.symlink(
        "/uufs/chpc.utah.edu/common/home/u6020307/Turbospectrum2019/DATA",
        _path("DATA")
    )

    if opacities is None:
        # Execute Turbospectrum's babsma_lu

        # Write the babsma control file.
        with resource_stream(__name__, "babsma.template") as fp:
            babsma_template = fp.read()
            if isinstance(babsma_template, bytes):
                babsma_template = babsma_template.decode("utf-8")

        babsma_contents = babsma_template.format(**kwds)
        input_path = _path("babsma.par")
        with open(input_path, "w") as fp:
            fp.write(babsma_contents)

        # Execute Turbospectrum's babsma
        process = subprocess.run(
            ["babsma_lu"],
            cwd=dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            input=babsma_contents,
            encoding="ascii"
        )
        if process.returncode != 0:
            raise RuntimeError(process.stderr)
            
    # Write the bsyn_lu control file.
    with resource_stream(__name__, "bsyn.template") as fp:
        bsyn_template = fp.read()
        if isinstance(bsyn_template, bytes):
            bsyn_template = bsyn_template.decode("utf-8")
            
    contents = bsyn_template.format(**kwds)
    input_path = _path("bsyn.par")
    with open(input_path, "w") as fp:
        fp.write(contents)

    # Execute Turbospectrum's bsyn_lu
    process = subprocess.run(
        ["bsyn_lu"],
        cwd=dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        input=contents,
        encoding="ascii"
    )
    if process.returncode != 0:
        raise RuntimeError(process.stderr)

    # Load the spectrum from the result file.
    data = np.loadtxt(_path(result_basename))

    wavelength, rectified_flux, flux = data.T
    result = OrderedDict([
        ("wavelength", wavelength),
        ("wavelength_unit", "Angstrom"),
        ("rectified_flux", rectified_flux),
        ("flux", flux),
        ("flux_unit", "1e-17 erg / (Angstrom cm2 s)"),
    ])
    
    meta = dict(
        dir=dir,
        # TODO: lots of stuff..
    )

    return (result, meta)