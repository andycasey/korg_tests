import os
import subprocess
from pkg_resources import resource_stream

from grok.radiative_transfer.utils import get_default_lambdas
from grok.utils import copy_or_write

def turbospectrum_bsyn(
        photosphere,
        transitions,
        lambdas=None,
        abundances=None,
        isotopes=None,
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

    with resource_stream(__name__, "bsyn.template") as fp:
        template = fp.read()
        if isinstance(template, bytes):
            template = template.decode("utf-8")
        
    transition_basename_format = "transitions_{i}"
    modelinput_basename = "model"
    modelopac_basename = "opac"
    result_basename = "result"
    
    marcs_file_flag = "marcs" in photosphere.meta["read_format"].lower()

    # Turbospectrum allows for transitions to be written into multiple different files.
    if not isinstance(transitions, (tuple, list)):
        transitions = (transitions, )
    
    T = len(transitions)
    transition_format = kwargs.get("transition_format", "turbospectrum")
    for i, each in enumerate(transitions):
        copy_or_write(
            each,
            _path(transition_basename_format.format(i=i)), 
            format=transition_format
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
        marcs_file_flag=marcs_file_flag,
        metallicity=0,          # TODO: parse from abundances
        alpha_fe=0,             # TODO: parse from abundances
        r_process_abundance=0,  # TODO: parse from abundances
        s_process_abundance=0,  # TODO: parse from abundances
        modelinput_basename=modelinput_basename,
        modelopac_basename=modelopac_basename,
        result_basename=result_basename,
        N_transition_paths=T,
        transition_paths_formatted="\n".join(
            [transition_basename_format.format(i=i) for i in range(T)]
        )
    )

    if abundances is not None or isotopes is not None:
        raise NotImplementedError        
    
    # Write the control file.
    contents = template.format(**kwds)
    input_path = _path("bsyn.par")
    with open(input_path, "w") as fp:
        fp.write(contents)

    # Execute Turbospectrum.
    p = subprocess.run(
        ["bsyn_lu"],
        cwd=dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        input=contents,
        encoding="ascii"
    )

    foo = p.returncode, p.stdout, p.stderr

    raise a
