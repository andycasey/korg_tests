
import os
import numpy as np
import subprocess
from collections import OrderedDict
from pkg_resources import resource_stream
from tempfile import mkdtemp
from time import time

from grok.transitions import Transitions
from grok.radiative_transfer.moog.io import parse_summary_synth_output
from grok.radiative_transfer.utils import get_default_lambdas
from grok.utils import copy_or_write

def moog_synthesize(
        photosphere,
        transitions,
        strong_transitions=None,
        lambdas=None,
        abundances=None,
        isotopes=None,
        terminal="x11",
        atmosphere_flag=0,
        molecules_flag=1,
        trudamp_flag=0,
        lines_flag=0,
        flux_int_flag=0,
        damping_flag=1,
        units_flag=0,
        scat_flag=1,
        opacit=0,
        opacity_contribution=2.0,
        n_chunks=1,
        verbose=False,
        dir=None,
        **kwargs
    ):
    
    if lambdas is not None:
        lambda_min, lambda_max, lambda_delta = lambdas
    else:
        lambda_min, lambda_max, lambda_delta = get_default_lambdas(transitions)

    N = 1 # number of syntheses to do
    
    _path = lambda basename: os.path.join(dir or "", basename)

    # Write photosphere and transitions.
    model_in, lines_in = (_path("model.in"), _path("lines.in"))
    copy_or_write(
        photosphere,
        model_in,
        format=kwargs.get("photosphere_format", "moog")
    )

    if isinstance(transitions, Transitions):
        # Cull transitions outside of the linelist, and sort. Otherwise MOOG dies.    
        use_transitions = Transitions(sorted(
            filter(
                lambda t: (lambda_max + opacity_contribution) >= t.lambda_air.value >= (lambda_min - opacity_contribution),
                transitions
            ),
            key=lambda t: t.lambda_air
        ))
    else:
        # You're living dangerously!
        use_transitions = transitions
    
    print(f"Writing {len(use_transitions)} to {lines_in}")

    copy_or_write(
        use_transitions,
        lines_in,
        format=kwargs.get("transitions_format", "moog")
    )
    
    with resource_stream(__name__, "moog_synth.template") as fp:
        template = fp.read()
    
        if isinstance(template, bytes):
            template = template.decode("utf-8")

    kwds = dict(
        terminal=terminal,
        atmosphere_flag=atmosphere_flag,
        molecules_flag=molecules_flag,
        trudamp_flag=trudamp_flag,
        lines_flag=lines_flag,
        damping_flag=damping_flag,
        flux_int_flag=flux_int_flag,
        units_flag=units_flag,
        scat_flag=scat_flag,
        opacit=opacit,
        opacity_contribution=opacity_contribution,
        lambda_min=lambda_min,
        lambda_max=lambda_max,
        lambda_delta=lambda_delta
    )
    if verbose:
        kwds.update(dict(
            atmosphere_flag=2,
            molecules_flag=2,
            lines_flag=3
        ))

    if strong_transitions is not None:
        stronglines_in = _path("stronglines.in")
        copy_or_write(
            strong_transitions,
            stronglines_in,
            format=kwargs.get("transitions_format", "moog"),
            # MOOG likes a blank header for a transition list, but not for a strong transition list.
            # because of course that makes sense
            include_header=False
        )
        kwds.update(
            stronglines_in=os.path.basename(stronglines_in),
            strong_flag=1,
        )
    else:
        kwds.update(
            stronglines_in="",
            strong_flag=0
        )
        
    # Abundances.
    if abundances is None:
        kwds["abundances_formatted"] = "0 1"
    else:
        raise NotImplementedError

    # Isotopes.
    if isotopes is None:
        kwds["isotopes_formatted"] = f"0 {N:.0f}"

    # I/O files:
    kwds.update(
        dict(
            standard_out="synth.std.out",
            summary_out="synth.sum.out",
            model_in=os.path.basename(model_in),
            lines_in=os.path.basename(lines_in)
        )
    )

    if n_chunks > 1:
        chunk_size = (lambda_max - lambda_min) / n_chunks

        wallclock_time = 0
        spectrum = OrderedDict([
            ("wavelength", []),
            ("wavelength_unit", "Angstrom"),
            ("rectified_flux", [])
        ])
        for chunk in range(n_chunks):
            lambda_min = lambda_min + chunk_size * chunk
            lambda_max = lambda_min + chunk_size * (chunk + 1)
            lines_in = _path(f"lines.in{chunk}")

            kwds.update(
                lines_in=os.path.basename(lines_in),
                lambda_min=lambda_min,
                lambda_max=lambda_max,
                standard_out=f"synth.std.out.{chunk}",
                summary_out=f"synth.sum.out.{chunk}"
            )

            chunk_transitions = Transitions(sorted(
                filter(
                    lambda t: (lambda_max + opacity_contribution) >= t.lambda_air.value >= (lambda_min - opacity_contribution),
                    use_transitions
                ),
                key=lambda t: t.lambda_air
            )) 

            copy_or_write(
                chunk_transitions,
                lines_in,
                format=kwargs.get("transitions_format", "moog")
            )                      

            # Write the control file.
            contents = template.format(**kwds)
            control_path = _path("batch.par")
            with open(control_path, "w") as fp:
                fp.write(contents)
            with open(_path(f"batch.par{chunk}"), "w") as fp:
                fp.write(contents)
            
            print(f"Executing chunk {chunk} for MOOGSILENT")

            t_init = time()
            process = subprocess.run(
                ["MOOGSILENT"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=os.path.dirname(control_path),
                input=os.path.basename(control_path) + "\n"*100,
                encoding="ascii"
            )
            wallclock_time += (time() - t_init)
            if process.returncode != 0:
                raise RuntimeError(process.stderr)    

            # Copy outputs.    
            output = parse_summary_synth_output(_path(kwds["summary_out"]))

            wavelength, rectified_flux, meta = output[0]
            spectrum["wavelength"].extend(wavelength)
            spectrum["rectified_flux"].extend(rectified_flux)
            meta["dir"] = dir

        meta["wallclock_time"] = wallclock_time
        return (spectrum, meta)

    else:
    
        # Write the control file.
        contents = template.format(**kwds)
        control_path = _path("batch.par")
        with open(control_path, "w") as fp:
            fp.write(contents)

        # Execute MOOG(SILENT).
        print(f"Executing MOOGSILENT in {_path('')}")

        t_init = time()
        process = subprocess.run(
            ["MOOGSILENT"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            cwd=os.path.dirname(control_path),
            input=os.path.basename(control_path) + "\n"*100,
            encoding="ascii"
        )
        t_moogsilent = time() - t_init
        if process.returncode != 0:
            raise RuntimeError(process.stderr)

        # Read the output.
        output = parse_summary_synth_output(_path(kwds["summary_out"]))
        wavelength, rectified_flux, meta = output[0]
        
        spectrum = OrderedDict([
            ("wavelength", wavelength),
            ("wavelength_unit", "Angstrom"),
            ("rectified_flux", rectified_flux),
        ])
        
        meta["dir"] = dir
        meta["wallclock_time"] = t_moogsilent
        
        return (spectrum, meta)
    
