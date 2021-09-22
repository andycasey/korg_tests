
import os
from tempfile import mkdtemp
from pkg_resources import resource_stream

from grok.radiative_transfer.moog.io import parse_summary_synth_output
from grok.radiative_transfer.moog.utils import moogsilent

def moog_synthesize(
        photosphere,
        transitions,
        wavelengths=None,
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
        mkdtemp_kwargs=None,
        **kwargs
    ):
    
    if wavelengths is not None:
        wavelength_start, wavelength_end, wavelength_delta = wavelengths
    else:
        wavelength_delta = 0.01
        wavelength_start = transitions["wavelength"].min() - 2
        wavelength_end = transitions["wavelength"].max() + 2

    N = 1 # number of syntheses to do

    # Create a working directory.
    mkdtemp_kwargs = mkdtemp_kwargs or dict()
    dir = mkdtemp(**mkdtemp_kwargs)
    
    _path = lambda basename: os.path.join(dir, basename)

    # Write photosphere and transitions.
    model_in, lines_in = (_path("model.in"), _path("lines.in"))
    photosphere.write(model_in, format=kwargs.get("photosphere_format", "moog"))

    # Cull transitions outside of the linelist, and sort. Otherwise MOOG dies.
    mask = \
            (transitions["wavelength"] >= (wavelength_start - opacity_contribution)) \
        *   (transitions["wavelength"] <= (wavelength_end + opacity_contribution))
    use_transitions = transitions[mask]
    use_transitions.sort("wavelength")
    use_transitions.write(lines_in, format=kwargs.get("transitions_format", "moog"))
    
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
        wavelength_start=wavelength_start,
        wavelength_end=wavelength_end,
        wavelength_delta=wavelength_delta
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

    # TODO: what the fuck moog
    meta["dir"] = dir
    return (wl, flux, meta)
    
