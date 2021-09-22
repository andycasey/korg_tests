from grok.radiative_transfer.moog.synthesize import moog_synthesize

def synthesize(
        photosphere,
        transitions,
        method,
        wavelengths=None,
        abundances=None,
        options=None,
        verbose=False,
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

    :param wavelengths: [optional]
        A three length tuple containing the start wavelength, end wavelength, and the
        wavelength step size. If `None` is given then this will default to the range
        of the transition list.

    :param options: [optional]
        A dictionary of options to provide that are specific to individual radiative
        transfer codes.

    :param verbose: [optional]
        Provide verbose outputs. The exact verbosity given depends on the radiative 
        transfer code used.
    """

    available_methods = {
        "moog": moog_synthesize,
    }
    try:
        _synthesize = available_methods[str(method).lower().strip()]
    except KeyError:
        raise ValueError(f"Method '{method}' unknown. Available methods: {', '.join(available_methods.keys())}")
    
    kwds = (options or dict())
    kwds.update(**kwargs)

    return _synthesize(
        photosphere,
        transitions,
        wavelengths=wavelengths,
        abundances=abundances,
        verbose=verbose,
        **kwds
    )
