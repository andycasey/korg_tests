
import os
from tempfile import mkdtemp
from pkg_resources import resource_stream

from grok.radiative_transfer.moog.io import parse_summary_synth_output
from grok.radiative_transfer.moog.utils import moogsilent
from grok.radiative_transfer.utils import get_default_lambdas
from grok.utils import copy_or_write

def moog_synthesize(
        photosphere,
        transitions,
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

    # Cull transitions outside of the linelist, and sort. Otherwise MOOG dies.
    mask = \
            (transitions["wavelength"] >= (lambda_min - opacity_contribution)) \
        *   (transitions["wavelength"] <= (lambda_max + opacity_contribution))
    use_transitions = transitions[mask]
    use_transitions.sort("wavelength")
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

    # Write the control file.
    contents = template.format(**kwds)
    with open(_path("batch.par"), "w") as fp:
        fp.write(contents)

    # Execute MOOG(SILENT).
    moogsilent(_path("batch.par"))

    # Read the output.
    spectra = parse_summary_synth_output(_path(kwds["summary_out"]))
    
    wl, flux, meta = spectra[0]

    meta["dir"] = dir
    return (wl, flux, meta)
    
