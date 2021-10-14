import os
import numpy as np
import subprocess
from collections import OrderedDict
from pkg_resources import resource_stream
from time import time

from grok.radiative_transfer.utils import get_default_lambdas
from grok.utils import copy_or_write



from grok.transitions.utils import (air_to_vacuum, vacuum_to_air)
from grok.utils import copy_or_write

def synthesize(
        photosphere,
        transitions,
        lambdas=None,
        dir=None,
        **kwargs
    ):
    """
    Execute synthesis with Korg.
    """

    transitions_path_basename = "transitions"
    photosphere_path_basename = "atmosphere"


    if lambdas is not None:
        lambda_air_min, lambda_air_max, lambda_air_delta = lambdas
    else:
        lambda_air_min, lambda_air_max, lambda_air_delta = get_default_lambdas(transitions)

    _path = lambda basename: os.path.join(dir or "", basename)

    copy_or_write(
        transitions,
        _path(transitions_path_basename),
        format=kwargs.get("transitions_format", "vald")
    )

    copy_or_write(
        photosphere,
        _path(photosphere_path_basename),
        format=kwargs.get("photosphere_format", "korg")
    )
    
    kwds = dict(
        # Korg works in vacuum wavelengths.
        lambda_vacuum_min=air_to_vacuum(lambda_air_min),
        lambda_vacuum_max=air_to_vacuum(lambda_air_max),
        atmosphere_path=photosphere_path_basename,
        linelist_path=transitions_path_basename,
        metallicity=photosphere.meta["m_h"]
    )

    # Write the template.
    with resource_stream(__name__, "template.jl") as fp:
        template = fp.read()
        if isinstance(template, bytes):
            template = template.decode("utf-8")

    contents = template.format(**kwds)
    input_path = _path("grok.jl")
    with open(input_path, "w") as fp:
        fp.write(contents)

    # Execute it in a subprocess like we do for Turbospectrum and MOOG.
    # (I don't know how Julia works, but this could introduce some weird time penalty on Julia. Should check on that.)
    t_init = time()
    process = subprocess.run(
        ["julia", "grok.jl"],
        cwd=dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        #input=babsma_contents,
        encoding="ascii"
    )
    t_korg = time() - t_init
    if process.returncode != 0:
        raise RuntimeError(process.stderr)    

    raise a
