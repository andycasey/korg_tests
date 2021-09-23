from tempfile import mkdtemp
from grok.radiative_transfer.moog.synthesize import moog_synthesize
from grok.radiative_transfer.turbospectrum.synthesize import turbospectrum_bsyn

def synthesize(
        photosphere,
        transitions,
        method,
        lambdas=None,
        abundances=None,
        options=None,
        verbose=False,
        mkdtemp_kwargs=None,
        **kwargs
    ):
    """
    Synthesize a stellar spectrum.

    :param photosphere:
        The photosphere to use.
    
    :param transitions:
        The list of transitions to include in synthesis.
    
    :param method:
        The radiative transfer code to use for synthesis.

    :param lambdas: [optional]
        A three length tuple containing the start wavelength, end wavelength, and the
        wavelength step size. If `None` is given then this will default to nearly the 
        range of the transition list, with a step size of 0.01 Angstroms.

    :param options: [optional]
        A dictionary of options to provide that are specific to individual radiative
        transfer codes.

    :param verbose: [optional]
        Provide verbose outputs. The exact verbosity given depends on the radiative 
        transfer code used.

    :param mkdtemp_kwargs: [optional]
        Keyword arguments to supply to `tempfile.mkdtemp` to create a temporary directory
        before executing the synthesis.
    """

    available_methods = {
        "moog": moog_synthesize,
        "turbospectrum": turbospectrum_bsyn
    }
    try:
        _synthesize = available_methods[str(method).lower().strip()]
    except KeyError:
        raise ValueError(f"Method '{method}' unknown. Available methods: {', '.join(available_methods.keys())}")

    options = options or dict()
    dir = mkdtemp(**(mkdtemp_kwargs or dict()))

    result = _synthesize(
        photosphere,
        transitions,
        lambdas=lambdas,
        abundances=abundances,
        verbose=verbose,
        dir=dir,
        **options
    )

    # TODO: Clean up temporary directory?
    return result
