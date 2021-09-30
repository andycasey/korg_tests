
import numpy as np
from astropy.constants import (e, m_e, c)

from grok.transitions.species import Species


class Line(object):

    """
    Represent an individual line transition.
        """

    def __init__(
            self,
            lambda_vacuum,
            log_gf,
            species,
            E_lower,
            gamma_rad=None,
            gamma_stark=None,
            vdW=None,
            # All the extra possible stuff.
            lambda_air=None,
            E_upper=None,
            j_lower=None,
            j_upper=None,
            lande_factor_lower=None,
            lande_factor_upper=None,
            lande_factor_mean=None,
            # TODO: What is the translation between lande_factor_lower/upper/mean and lande_factor_depth?
            lande_factor=None,
            lande_depth=None,
            reference=None,
            comment=None,
            equivalent_width=None,
            equivalent_width_error=None,
            **kwargs
        ):
                
        if lambda_vacuum is None:
            # Calculate it from lambda_air.
            if lambda_air is None:
                raise ValueError("Wavelength (lambda) must be given as lambda_vacuum or lambda_air.")
            raise NotImplementedError
        
        if not isinstance(species, Species):
            species = Species(species)

        if gamma_stark is None or vdW is None:
            gamma_stark_approx, vdW_approx = approximate_gammas(lambda_vacuum, species, E_lower)
            gamma_stark = gamma_stark or gamma_stark_approx
            vdW = vdW or vdW_approx

        if vdW < 0:
            # If the van der Waals constant is negative, we assume it is log(\Gamma_vdW)
            vdW = 10**vdW
        elif vdW > 1:
            # If the van der Waals constant is > 1 we assume that it's packed ABO
            raise NotImplementedError("check to see how this is packed in MOOG vs others")
            

        gamma_rad = gamma_rad \
                    or approximate_radiative_gamma(lambda_vacuum, log_gf)
                
        # Store everything.
        self.lambda_vacuum = lambda_vacuum
        self.species = species
        self.log_gf = log_gf
        self.E_lower = E_lower

        self.gamma_stark = gamma_stark
        self.gamma_rad = gamma_rad
        self.vdW = vdW
        self.log_gf = log_gf
        
        self.E_upper = E_upper
        self.j_lower = j_lower
        self.j_upper = j_upper
        self.lande_factor_lower = lande_factor_lower
        self.lande_factor_upper = lande_factor_upper
        self.lande_factor_mean = lande_factor_mean
        # TODO: What is the translation between lande_factor_lower/upper/mean and lande_factor_depth?
        self.lande_factor = lande_factor
        self.lande_depth = lande_depth
        self.reference = reference
        self.comment = comment
        self.equivalent_width = equivalent_width
        self.equivalent_width_error = equivalent_width_error

        return None


def approximate_radiative_gamma(lambda_vacuum, log_gf):
    # TODO: assumes lambda_vacuum is in cm, or has a unit.
    return 8 * np.pi**2 * e.cgs**2 / (m.cgs * c.cgs * lambda_vacuum**2) * 10**log_gf


def approximate_gammas(lambda_vacuum, species, ionization_energies=None):
    raise a