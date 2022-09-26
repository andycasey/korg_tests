
from tempfile import tempdir
from pysme.linelist.linelist import LineList
from pysme.linelist.vald import ValdFile

from pysme.sme import SME_Structure
from pysme.atmosphere.atmosphere import Atmosphere
from pysme.abund import Abund
from pysme.synthesize import synthesize_spectrum
from time import time
from collections import OrderedDict
from typing import Union, Tuple, List
import numpy as np

from grok.transitions import Transitions, Transition
from grok.photospheres import Photosphere
from grok.synthesis.utils import get_default_lambdas

def repr_species(species):
    atoms = "".join(species.atoms)
    translate = {
        "HFe": "FeH",
        "HO": "OH", 
        "CC": "C2"
    }
    atoms = translate.get(atoms, atoms)
    return f"{atoms} {species.charge + 1}"

def sme_synthesize(
    photosphere: Union[str, Photosphere],
    transitions: Union[str, List[str], Tuple[str]],
    abundances,
    lambdas: Tuple[float, float, float],
    input_transitions_format="vald",
    **kwargs
):
    
    # TODO: SME expects air wavelengths, but lambda_min etc is vacuum
    lambda_min, lambda_max, lambda_delta = lambdas
    N = int((lambda_max - lambda_min) / lambda_delta)
    wavelength = np.linspace(lambda_min, lambda_max, 1 + N)

    sme = SME_Structure()
    sme.atmo = Atmosphere(
        teff=photosphere.meta["teff"],
        logg=photosphere.meta["logg"],
        method="embedded",
        geom="SPH" if photosphere.is_spherical_geometry else "PP",
        radius=photosphere.meta["radius"],
        height=np.zeros(len(photosphere)),
        vturb=photosphere.meta["microturbulence"],
        lonh=photosphere.meta["convection_alpha"],
        wlstd=5000,
        depth="RHOX",
        interp="TAU",
        rhox=photosphere["RHOX"],
        tau=10**photosphere["lgTau5"],
        temp=photosphere["T"],
        rho=photosphere["Density"],
        xne=photosphere["XNE"],
        xna=photosphere["numdens_other"],
        # abund
        abund=np.array([v for k, v in photosphere.meta.items() if k.startswith("log_abundance_")]),
        abund_format="H=12", # TODO: Check that this is followed, otherwise convert abundances to SME format.
        # I downloaded the MARCS models from the SME website and confirmed that every photosphere had these opflags
        # for the standard composition (see check_sme_photospheres.py)
        opflag=np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0]), 
    )
    sme.teff = photosphere.meta["teff"]
    sme.logg = photosphere.meta["logg"]
    sme.monh = photosphere.meta["m_h"]
    #sme.abund = #a # TODO: What's actually used?
    sme.vmic = photosphere.meta["microturbulence"]
    sme.vmac = 0
    sme.vsini = 0


    sme.abund = Abund(photosphere.meta["m_h"], "asplund2009")
    
    if isinstance(transitions, (Transition, Transitions)):
        linelist = to_linelist(transitions)
    else:
        if isinstance(transitions, str):
            transitions = [transitions]
        if input_transitions_format == "vald":
            linelist = ValdFile(transitions[0])
            for each in transitions[1:]:
                linelist.append(ValdFile(each))
        else:
            _transitions = []
            for each in transitions:
                _transitions.extend(Transitions.read(each))

            linelist = to_linelist(Transitions(_transitions))

    sme.linelist = linelist

    # SME expects air wavelengths.
    sme.wave = wavelength

    sme.normalize_by_continuum = True
    
    # TODO: Where to start timing from...?
    t_init = time()
    sme = synthesize_spectrum(sme)
    t_synthesis = time() - t_init

    spectrum = OrderedDict([
        ("wavelength", sme.wave[0]),
        ("wavelength_unit", "Angstrom"),
        ("rectified_flux", sme.synth[0])
    ])
    timing = dict(
        t_synthesis=t_synthesis, 
    )
    meta = dict(sme_structure=sme.to_dict())
    return (spectrum, timing, meta)


def to_linelist(transitions):
    linedata = []
    kwargs = {}
    for k in ( "term_lower", "term_upper"):#, "j_lo", "e_upp", "j_up", "lande_lower", "lande_upper", "lande"):
        kwargs[k] = []
    for t in transitions:
        linedata.append([
            repr_species(t.species),
            t.lambda_vacuum.value,
            t.log_gf,
            t.E_lower.value,
            t.gamma_rad.value,
            t.gamma_stark.value,
            t.vdW,
            t.j_lower,
            t.E_upper.value if t.E_upper is not None else 0,
            t.j_upper,
            t.lande_factor_lower,
            t.lande_factor_upper,
            t.lande_factor_mean,
            t.species.charge + 1
        ])
        #kwargs["j_lo"].append(t.j_lower)
        #kwargs["e_upp"].append(t.E_upper.value)
        #kwargs["j_up"].append(t.j_upper)
        #kwargs["lande_lower"].append(t.lande_factor_lower)
        #kwargs["lande_upper"].append(t.lande_factor_upper)
        #kwargs["lande"].append(t.lande_factor_mean)
        #kwargs["ionization"].append(t.species.charge + 1)
        try:
            other, lower, upper = t.comment.split(":")
            lower = " ".join(lower.split(" ")[:-1])
        except:
            kwargs["term_lower"].append(" ")
            kwargs["term_upper"].append(" ")
        else:
            kwargs["term_lower"].append(lower)
            kwargs["term_upper"].append(upper)

    from pandas import DataFrame
    linedata = DataFrame(
        linedata,
        columns=["species", "wlcent", "gflog", "excit", "gamrad", "gamqst", "gamvw",
            "j_lo",
            "e_upp",
            "j_up",
            "lande_lower",
            "lande_upper",
            "lande",
            "ionization"
            ]
    )
    linelist = LineList(linedata, **kwargs)
    linelist.sort("wlcent")
    return linelist
