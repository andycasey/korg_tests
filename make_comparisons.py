
"""
Compare synthesise methods for some benchmark stars.
"""

import os
from re import A
import logging as log
import numpy as np
import pickle
import matplotlib.pyplot as plt
import yaml
from collections import OrderedDict
from distutils.dir_util import copy_tree
            
from astropy.convolution import (Gaussian1DKernel, convolve)
from glob import glob
from tqdm import tqdm
from grok import (Transitions, Photosphere, PhotosphereInterpolator)
from grok.photospheres.interpolator import NewPhotosphereInterpolator
from grok.radiative_transfer.turbospectrum.transitions import should_keep
from grok.radiative_transfer import synthesize
from grok.utils import read_benchmark_star


OVERWRITE = False

# Combinations of `star`, `transitions`, and `methods`.
with open("comparisons.yml", "r") as fp:
    config = yaml.load(fp, Loader=yaml.FullLoader)

output_prefix = lambda star_description, lambda_min, lambda_max, method_description: f"{star_description}_{lambda_min:.0f}_{lambda_max:.0f}_{method_description}"

window = 2

for method_description, options in config["methods"].items():
    
    method, *desc = method_description.split("_")
    if not method.startswith("moog"): 
        print("SKIPPING")
        continue

    for star_description, star in config["stars"].items():

        photosphere = None

        for transition_kwds in config["transitions"]:

            # Copy options instance
            options_instance = {}
            options_instance.update(options)

            lambda_min, lambda_max, lambda_step = lambdas = transition_kwds["lambdas"]

            print(f"Checking {method_description} for {star_description} and {lambda_min} to {lambda_max}")

            basename = output_prefix(star_description, lambda_min, lambda_max, method_description)

            output_path = f"executions/{basename}.pkl"
            output_dir = f"executions/{basename}"
            if os.path.exists(output_path) and not OVERWRITE:
                print(f"Skipping {method_description} on {star_description} from {lambda_min} to {lambda_max}")
                continue

            if photosphere is None:
                print(f"Creating photosphere..")
                photospheres = []
                for i, path in enumerate(tqdm(glob(star["model_kwargs"]["photosphere_grid_wildmask"]), desc="Loading photospheres")):
                    photospheres.append(
                        Photosphere.read(
                            path, 
                            format=star["model_kwargs"].get("photosphere_read_format", None)
                        )
                    )

                #photosphere_interpolator = PhotosphereInterpolator(
                #    photospheres, 
                #    interpolate_log_quantities=("P", "XNE", "Pg", "Pe", ),
                #)
                photosphere_interpolator = NewPhotosphereInterpolator(photospheres)
                photosphere = photosphere_interpolator(**star["model_kwargs"]["photosphere_point"])

            else:
                print(f"Using photosphere calculated for {star_description} from the last region")
            

            # Some bespoke handling options for transitions.
            if method == "moog":
                # MOOG simply cannot handle a huge list of transitions.
                # Instead we will restrict the TiO transitions to those that are very strong,
                # up until a maximum number of transitions.
                
                # Using max_transitions = 25_000 requires us to do ~50 chunks for the bluest region AND to
                # set opacity_contribution = 1.0, so it just cannot be practically done because MOOG fails.
                max_transitions_and_chunks = {
                    3660: (10_000, 2),
                    3930: (10_000, 2),
                    5160: (10_000, 2), # (25k, 10) fails
                    6540: (15_000, 2), # (25k, 20) fails 
                    15000: (100_000, 20) 
                }

                max_transitions, n_chunks = max_transitions_and_chunks[lambda_min]
                
                restrict_on = lambda t: sorted(t.species.atoms) == ["O", "Ti"]

                teff = star["model_kwargs"]["photosphere_point"]["teff"]
                def strength(t, teff):
                    return t.log_gf - t.E_lower.value * (5040/teff) + np.log(t.lambda_air.value)

                _transitions = []
                for path in transition_kwds["paths"]:
                    print(f"Reading in transitions from {path}")
                    _transitions.extend(Transitions.read(path))

                print(f"There were {len(_transitions)} before removing strange lines and things outside the range")
                _transitions = [t for t in _transitions \
                        if should_keep(t) and (lambda_max + window) >= t.lambda_air.value >= (lambda_min - window)
                ]

                if len(_transitions) > max_transitions:
                    print(f"There were {len(_transitions)} before restricting on strength")

                    non_restricted, S = 0, []
                    for t in _transitions:
                        if restrict_on(t):
                            S.append(strength(t, teff))
                        else:
                            non_restricted += 1

                    if non_restricted > max_transitions:
                        raise RuntimeError("Need to restrict on something else.")
                            
                    min_S = np.sort(S)[-max_transitions + non_restricted]

                    transitions = [t for t in _transitions if not restrict_on(t) or strength(t, teff) >= min_S]
                    print(f"Now there are {len(transitions)} after restricting on strength")
                    
                else:
                    print(f"Not restricting on strength")
                    transitions = _transitions

                # Separate transitions and strong transitions for MOOG
                strong_path = transition_kwds.get("strong_path", None)
                if strong_path is None:
                    strong_transitions = None
                    N_strong_transitions = 0
                else:
                    strong_transitions_format = "vald" if strong_path.endswith(".vald") else None
                    strong_transitions = Transitions.read(strong_path, format=strong_transitions_format)
                    N_strong_transitions = len(strong_transitions)
                    before = len(transitions)
                    transitions = [t for t in transitions if t not in strong_transitions]
                    after = len(transitions)
                    assert (after - before) <= N_strong_transitions
                    print(f"Now there are {len(transitions)} after adding {N_strong_transitions} strong transitions")
            
                transitions = Transitions(transitions)
                options_instance.update(
                    strong_transitions=strong_transitions,
                    n_chunks=n_chunks
                )
                
                transitions_desc = f"using {len(transitions)} transitions and {N_strong_transitions} strong transitions"

            elif method == "korg":# and lambda_min != 15_000:
                # If we're running Korg and it's not the (15,000 - 15,500) region then the
                # line list is already in VALD format, so we will supply this directly to Korg.
                # If it's not that region, we need to load things.
                transitions = transition_kwds["paths"]
                transitions_desc = f"using these line lists: {', '.join(transitions)}"

                if len(transitions) == 1:
                    # This is annoying, but Korg doesn't like it if we pass it a list of paths.
                    transitions = transitions[0]

            elif method == "turbospectrum" and lambda_min == 15_000:
                # Supply line list directly to TS, since it is already in TS format.
                transitions = transition_kwds["paths"]
                transitions_desc = f"using these line lists: {', '.join(transitions)}"
            
            else:
                _transitions = []
                for path in transition_kwds["paths"]:
                    print(f"Reading in transitions from {path}")
                    _transitions.extend(Transitions.read(path))
                
                # Restrict to things within the wavelength region, which Turbospectrum would use.
                # (Otherwise Turbospectrum will die, and there's no point sending things TS would not use to MOOG)
                transitions = Transitions([
                    t for t in _transitions \
                        if should_keep(t) and (lambda_max + window) >= t.lambda_air.value >= (lambda_min - window)
                ])
                transitions_desc = f"using {len(transitions)} transitions"
            
            # Check if we are supplying hydrogen lines to this region, and whether we should be.
            #if not transition_kwds["has_hydrogen_lines"] and "hydrogen_lines" in options_instance:
            #    print("Turning off hydrogen lines because there aren't any in this region")
            #    options_instance["hydrogen_lines"] = False

            print(f"Executing {method} ({method_description}) for {star_description} between {lambda_min} and {lambda_max} {transitions_desc} with {options_instance}")

            try:
                spectrum, meta = synthesize(
                    photosphere,
                    transitions,
                    method=method,
                    lambdas=lambdas,
                    options=options_instance,
                )

            except:
                log.exception("Exception occurred")
                raise

            print(f"  That took {meta['wallclock_time']} seconds")

            with open(output_path, "wb") as fp:
                pickle.dump((spectrum, meta), fp)
            
            copy_tree(meta["dir"], output_dir)
            
            print(f"  Wrote to {output_path}")
            print(f"  Copied execution contents to {output_dir}")

raise RuntimeError("Done")