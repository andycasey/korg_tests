import os
import numpy as np
import subprocess
import re
from collections import OrderedDict
from pkg_resources import resource_stream
from time import time
from astropy import units as u

from grok.radiative_transfer.utils import get_default_lambdas
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

    julia_script_basename = "grok.jl"
    photosphere_path_basename = "atmosphere"

    if lambdas is not None:
        lambda_air_min, lambda_air_max, lambda_air_delta = lambdas
    else:
        lambda_air_min, lambda_air_max, lambda_air_delta = get_default_lambdas(transitions)

    _path = lambda basename: os.path.join(dir or "", basename)

    # Write the photosphere first.
    copy_or_write(
        photosphere,
        _path(photosphere_path_basename),
        format=kwargs.get("photosphere_format", "kurucz")
    )    
    
    kwds = dict(
        # Korg works in vacuum wavelengths.
        # TODO: Don't assume units for lambdas.
        lambda_vacuum_min=air_to_vacuum(lambda_air_min * u.Angstrom).to("Angstrom").value,
        lambda_vacuum_max=air_to_vacuum(lambda_air_max * u.Angstrom).to("Angstrom").value,
        atmosphere_path=photosphere_path_basename,
        metallicity=photosphere.meta["m_h"],
        # I'm just giving a different metallicity for the initial run so that people don't 
        # incorrectly think the result from the second run is actually being cached.
        fake_metallicity=photosphere.meta["m_h"] - 0.25,
        korg_read_transitions_format=kwargs.get("korg_read_transitions_format", "vald")
    )

    # I wish I knew some Julia... eeek!
    if isinstance(transitions, (list, tuple)) and len(transitions) == 2:
        template_path = "template_two_linelists.jl"
        for i, each in enumerate(transitions):
            basename = f"transitions_{i:.0f}"
            copy_or_write(
                each,
                _path(basename),
                format=kwargs.get("transitions_format", "vald.stellar")
            )
            kwds[f"linelist_path_{i}"] = basename

    else:
        template_path = "template.jl"
        transitions_path_basename = "transitions"
        copy_or_write(
            transitions,
            _path(transitions_path_basename),
            format=kwargs.get("transitions_format", "vald.stellar")
        )
        kwds["linelist_path"] = transitions_path_basename

    # Write the fake line list template.
    with resource_stream(__name__, "weak_linelist.template") as fp:
        template = fp.read()
        if isinstance(template, bytes):
            template = template.decode("utf-8")

    with open(_path("linelist.fake"), "w") as fp:
        fp.write(template)

    # Write the template.
    with resource_stream(__name__, template_path) as fp:
        template = fp.read()
        if isinstance(template, bytes):
            template = template.decode("utf-8")

    contents = template.format(**kwds)
    input_path = _path(julia_script_basename)
    with open(input_path, "w") as fp:
        fp.write(contents)

    # Execute it in a subprocess like we do for Turbospectrum and MOOG.
    # (I don't know how Julia works, but this could introduce some weird time penalty on Julia. Should check on that.)
    t_init = time()
    process = subprocess.run([
            "julia", 
            julia_script_basename
        ],
        cwd=dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding="utf-8"
    )
    t_subprocess = time() - t_init
    if process.returncode != 0:
        raise RuntimeError(process.stderr)    
    
    # Write the stdout and stderr.
    with open(_path("stdout"), "w") as fp:
        fp.write(process.stdout)
    with open(_path("stderr"), "w") as fp:
        fp.write(process.stderr)

    # Parse for the synthesis time.
    pattern = "Synthesizing..\n\s*(?P<seconds>[\d\.]+) seconds"
    search = re.finditer(pattern, process.stdout)
    t_first = float(next(search).groupdict()["seconds"])
    t_second = float(next(search).groupdict()["seconds"])
    
    flux = np.loadtxt(_path("spectrum.out"))
    wavelength = np.arange(lambda_air_min, lambda_air_max + lambda_air_delta, lambda_air_delta)[:len(flux)]

    continuum = np.loadtxt(_path("continuum.out"))


    result = OrderedDict([
        ("wavelength", wavelength),
        ("wavelength_unit", "Angstrom"),
        ("flux", flux),
        ("flux_unit", "erg / (Angstrom cm2 s"),
        ("continuum", continuum),
        ("rectified_flux", flux / continuum)
    ])

    meta = dict(
        dir=dir,
        wallclock_time=t_second,
        t_first=t_first,
        t_subprocess=t_subprocess,
        stdout=process.stdout,
        stderr=process.stderr
    )
    # Parse the version used.
    match = re.search("Korg v(?P<korg_version_major>\d+)\.(?P<korg_version_minor>\d+)\.(?P<korg_version_micro>\d+)", process.stdout)
    meta.update(match.groupdict())

    return (result, meta)