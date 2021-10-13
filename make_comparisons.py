
"""
Compare Solar spectrum with Castelli/Kurucz (ATLAS9) model photospheres.
"""

import os
from re import A
import logging as log
import numpy as np
import pickle
import matplotlib.pyplot as plt
from distutils.dir_util import copy_tree
            
from astropy.convolution import (Gaussian1DKernel, convolve)
from glob import glob
from tqdm import tqdm
from grok import (Transitions, Photosphere, PhotosphereInterpolator)
from grok.radiative_transfer.turbospectrum.transitions import should_keep
from grok.radiative_transfer import synthesize
from grok.utils import read_benchmark_star

# Combinations of `star`, `transitions`, and `radiative transfer code`.

OVERWRITE = False

stars = [
    # Sun.
    dict(
        description="Sun",
        observation_path="data/ATLAS.Sun_47000.fits",
        model_kwargs=dict(  
            # NOTE: Not using interpolation.
            #photosphere_path="data/photospheres/castelli-kurucz-2004/ap00t5750g45k2odfnew.dat",
            photosphere_grid_wildmask="data/photospheres/castelli-kurucz-2004/*.dat",
            #photosphere_read_format="atlas9",
            photosphere_point=dict(
                teff=5777,
                logg=4.4,
                m_h=0,
                alpha_m=0
            )
        )
    ),   
    # Arcturus
    dict(
        description="Arcturus",
        observation_path="data/ATLAS.Arcturus_47000.fits",
        model_kwargs=dict(
            # NOTE: Not using interpolation.
            #photosphere_path = "data/photospheres/marcs_mod/s4250_g+1.5_m1.0_t02_st_z-0.50_a+0.20_c+0.00_n+0.00_o+0.20_r+0.00_s+0.00.mod.gz",
            #photosphere_path = "data/photospheres/castelli-kurucz-2004/am05t4250g15k2odfnew.dat",
            photosphere_grid_wildmask="data/photospheres/marcs_mod/s*_m1.0_t02_st_*a+0.20*",
            #photosphere_grid_wildmask="data/photospheres/castelli-kurucz-2004/*.dat",
            #photosphere_read_format="atlas9",
            photosphere_point=dict(
                teff=4286,
                logg=1.64,
                m_h=-0.53, 
                alpha_m=0,
                neighbours=8 # for interpolation.
            )
        )
    ),
    # HD49933
    dict(
        description="HD49933",
        observation_path="data/HARPS.Archive_HD49933_47000.fits",
        model_kwargs=dict(
            #photosphere_path="data/photospheres/castelli-kurucz-2004/am05t6500g40k2odfnew.dat",
            photosphere_grid_wildmask="data/photospheres/castelli-kurucz-2004/*.dat",
            #photosphere_read_format="atlas9",
            photosphere_point=dict(
                teff=6635,
                logg=4.20,
                m_h=-0.46,
                alpha_m=0
            )
        )
    ),
    # HD122563
    dict(
        description="HD122563",
        observation_path="data/ESPaDOnS_HD122563_47000.fits",
        model_kwargs=dict(
            # TODO: Switch to MARCS.
            photosphere_grid_wildmask="data/photospheres/marcs_mod/s*_m1.0_t02_st_*a+0.40*",
            #photosphere_path="data/photospheres/marcs_mod/s4500_g+1.5_m1.0_t02_st_z-2.50_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod.gz",
            #photosphere_grid_wildmask="data/photospheres/castelli-kurucz-2004/*.dat",
            #photosphere_read_format="atlas9",
            photosphere_point=dict(
                teff=4587,
                logg=1.61,
                m_h=-2.74, 
                alpha_m=+0.4,
            )
        )
    )

]


all_transition_kwds = [
    dict(
        paths=(
            "data/transitions/all-3660-3680.vald", 
        ),
        strong_path=None,
        lambdas=(3660, 3680, 0.01),
    ),
    dict(
        paths=(
            "data/transitions/all-3930-3950.vald", 
        ),
        # Ca II lines.
        strong_path="data/transitions/strong-3930-3950.vald",
        lambdas=(3930, 3950, 0.01),
    ),    
    dict(
        paths=(
            "data/transitions/all-5160-5176.vald", 
            "data/transitions/all-5176-5190.vald"
        ),
        # Mg II lines.
        strong_path="data/transitions/strong-5160-5190.vald",
        lambdas=(5160, 5190, 0.01),
    ),
    dict(
        paths=(
            "data/transitions/all-6540-6559.vald",
            "data/transitions/all-6559-6578.vald"
        ),
        strong_path=None,
        lambdas=(6540, 6578, 0.01),
    ),
]

methods = {
    "moog": dict(),
    "turbospectrum": dict(
        skip_irrelevant_transitions=True,
        update_missing_data=True
    ),
}

output_prefix = lambda star_description, lambda_min, lambda_max, method: f"{star_description}_{lambda_min:.0f}_{lambda_max:.0f}_{method}_interpolated"


failed = []
for star_kwds in stars:
    star_description = star_kwds["description"]

    if not OVERWRITE:
        any_missing = False
        for j, transition_kwds in enumerate(all_transition_kwds):
            lambda_min, lambda_max, lambda_step = lambdas = transition_kwds["lambdas"]

            for method in methods.keys():                
                output_path = f"executions/{output_prefix(star_description, lambda_min, lambda_max, method)}.pkl"
                
                if os.path.exists(output_path):
                    continue
                else:
                    any_missing = True
                    break
                
            if any_missing:
                break
        
        if not any_missing:
            continue


    star_description = star_kwds["description"]

    if "photosphere_path" in star_kwds["model_kwargs"]:
        photosphere = Photosphere.read(star_kwds["model_kwargs"]["photosphere_path"])

    else:
        # Load the photosphere.
        photospheres = []
        for i, path in enumerate(tqdm(glob(star_kwds["model_kwargs"]["photosphere_grid_wildmask"]), desc="Loading photospheres")):
            photospheres.append(
                Photosphere.read(
                    path, 
                    format=star_kwds["model_kwargs"].get("photosphere_read_format", None)
                )
            )

        photosphere_interpolator = PhotosphereInterpolator(
            photospheres, 
            interpolate_log_quantities=("P", "XNE", "Pg", "Pe", ),
            #filter_values=dict(
            #    # When interpolating on VELSND, only keep nearby photospheres that have this value.
            #    VELSND=lambda values: np.all(values > 0, axis=1),
            #    FLXCNV=lambda values: np.any(values > 0, axis=1)
            #)
        )

        photosphere = photosphere_interpolator(
            **star_kwds["model_kwargs"]["photosphere_point"],
        )
    
    for j, transition_kwds in enumerate(all_transition_kwds):

        all_transitions = []
        for path in transition_kwds["paths"]:
            transitions = Transitions.read(path, format="vald")

            # Restrict to the same set of transitions that Turbospectrum would use.
            transitions = Transitions([t for t in transitions if should_keep(t)])

            all_transitions.extend(transitions)

        strong_path = transition_kwds.get("strong_path", None)
        if strong_path is None:
            strong_transitions = None
        else:
            strong_transitions = Transitions.read(strong_path, format="vald")
        

        lambda_min, lambda_max, lambda_step = lambdas = transition_kwds["lambdas"]

        for method, options in methods.items():
            
            if star_description == "HD49933" and method == "moog":
                options["n_chunks"] = 2
            basename = output_prefix(star_description, lambda_min, lambda_max, method)
            
            output_path = f"executions/{basename}.pkl"
            output_dir = f"executions/{basename}"
            
            if os.path.exists(output_path) and not OVERWRITE:
                print(f"Skipping {star_description} between {lambda_min} and {lambda_max} as {output_path} exists")
                continue

            if method == "moog":
                # MOOG simply cannot handle a huge list of transitions.
                # Instead we will restrict the TiO transitions to those that are very strong,
                # up until a maximum number of transitions.
                max_transitions = 10_000
                restrict_on = lambda t: sorted(t.species.atoms) == ["O", "Ti"]

                teff = star_kwds["model_kwargs"]["photosphere_point"]["teff"]
                def strength(t, teff):
                    return t.log_gf - t.E_lower.value * (5040/teff) + np.log(t.lambda_air.value)

                non_restricted, S = 0, []
                for t in all_transitions:
                    if restrict_on(t):
                        S.append(strength(t, teff))
                    else:
                        non_restricted += 1

                if non_restricted > max_transitions:
                    raise RuntimeError("Need to restrict on something else.")
                        
                min_S = np.sort(S)[-max_transitions + non_restricted]

                transitions = [t for t in all_transitions if not restrict_on(t) or strength(t, teff) >= min_S]

                # Separate transitions and strong transitions for MOOG
                if strong_transitions is not None:
                    transitions = [t for t in transitions if t not in strong_transitions]
                
                #transitions = Transitions([t for t in all_transitions if not t.is_molecule])
                transitions = Transitions(transitions)
                options.update(strong_transitions=strong_transitions)

            else:
                transitions = Transitions(all_transitions)

            print(f"Executing {method} for {star_description} between {lambda_min} and {lambda_max} with {len(transitions)} transitions")

            try:
                    
                spectrum, meta = synthesize(
                    photosphere,
                    transitions,
                    method=method,
                    lambdas=lambdas,
                    options=options,
                    # I ran an example with Turbospectrum that took 5 seconds.
                    # MOOG took 15 hours before I stopped it.
                    # With verbose mode, the text log file was 242 GB hahahahaha
                    # ....but still no spectra
                    #verbose=True
                )

            except:
                log.exception("Exception occurred")
                failed.append((star_description, lambdas, method))
                raise
                continue

            print(f"That took {meta['wallclock_time']} seconds")

            with open(output_path, "wb") as fp:
                pickle.dump((spectrum, meta), fp)
            
            copy_tree(meta["dir"], output_dir)
            
            print(f"Wrote to {output_path}")
            print(f"Copied execution contents to {output_dir}")
