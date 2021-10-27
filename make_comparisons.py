
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
from grok.radiative_transfer.turbospectrum.transitions import should_keep
from grok.radiative_transfer import synthesize
from grok.utils import read_benchmark_star


OVERWRITE = False

# Combinations of `star`, `transitions`, and `methods`.
with open("comparisons.yml", "r") as fp:
    config = yaml.load(fp, Loader=yaml.FullLoader)

output_prefix = lambda star_description, lambda_min, lambda_max, method: f"{star_description}_{lambda_min:.0f}_{lambda_max:.0f}_{method}-alpha"

window = 2

for method, options in config["methods"].items():

    for star_description, star in config["stars"].items():

        photosphere = None

        for transition_kwds in config["transitions"]:
            
            lambda_min, lambda_max, lambda_step = lambdas = transition_kwds["lambdas"]

            print(f"Checking {method} for {star_description} and {lambda_min} to {lambda_max}")

            basename = output_prefix(star_description, lambda_min, lambda_max, method)

            output_path = f"executions/{basename}.pkl"
            output_dir = f"executions/{basename}"
            if os.path.exists(output_path) and not OVERWRITE:
                print(f"Skipping {method} on {star_description} from {lambda_min} to {lambda_max}")
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

                photosphere_interpolator = PhotosphereInterpolator(
                    photospheres, 
                    interpolate_log_quantities=("P", "XNE", "Pg", "Pe", ),
                )

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
                    strong_transitions = Transitions.read(strong_path, format="vald")
                    N_strong_transitions = len(strong_transitions)
                    before = len(transitions)
                    transitions = [t for t in transitions if t not in strong_transitions]
                    after = len(transitions)
                    assert (after - before) <= N_strong_transitions
                    print(f"Now there are {len(transitions)} after adding {N_strong_transitions} strong transitions")
            
                transitions = Transitions(transitions)
                options.update(
                    strong_transitions=strong_transitions,
                    n_chunks=n_chunks
                )
                
                transitions_desc = f"using {len(transitions)} transitions and {N_strong_transitions} strong transitions"

            elif method == "korg" and lambda_min != 15_000:
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

            print(f"Executing {method} for {star_description} between {lambda_min} and {lambda_max} {transitions_desc}")

            try:
                spectrum, meta = synthesize(
                    photosphere,
                    transitions,
                    method=method,
                    lambdas=lambdas,
                    options=options,
                )

            except:
                log.exception("Exception occurred")
                raise
                continue

            print(f"  That took {meta['wallclock_time']} seconds")

            with open(output_path, "wb") as fp:
                pickle.dump((spectrum, meta), fp)
            
            copy_tree(meta["dir"], output_dir)
            
            print(f"  Wrote to {output_path}")
            print(f"  Copied execution contents to {output_dir}")

raise RuntimeError("Done")


failed = []
for star_kwds in config["stars"]:
    star_description = star_kwds["description"]

    if not OVERWRITE:
        any_missing = False
        for j, transition_kwds in enumerate(config["transitions"]):
            lambda_min, lambda_max, lambda_step = lambdas = transition_kwds["lambdas"]

            for method in config["methods"].keys():                
                output_path = f"executions/{output_prefix(star_description, lambda_min, lambda_max, method)}.pkl"

                
                print(f"Checking {output_path}")
                if os.path.exists(output_path):
                    print("  present")
                    continue
                else:
                    print("  missing")
                    any_missing = True
                    break
                
            if any_missing:
                break
        
        if not any_missing:
            print(f"Skipping {star_description} because we have all files")
            continue


    star_description = star_kwds["description"]

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
    )

    photosphere = photosphere_interpolator(
        **star_kwds["model_kwargs"]["photosphere_point"],
    )
    
    for j, transition_kwds in enumerate(config["transitions"]):

        N_total_possible_transitions = 0
        all_transitions = []
        for path in transition_kwds["paths"]:
            transitions = Transitions.read(path)

            print(f"Read in {len(transitions)} from {path}")
            N_total_possible_transitions += len(transitions)

            # Restrict to the same set of transitions that Turbospectrum would use.
            transitions = Transitions([t for t in transitions if should_keep(t)])

            # Restrict to things +/- 1 Angstrom in the lambdas 
            lambda_min, lambda_max, lambda_delta = transition_kwds["lambdas"]
            transitions = Transitions([t for t in transitions if (lambda_max + 1) >= t.lambda_air.value >= (lambda_min - 1)])

            all_transitions.extend(transitions)

        print(f"There were a total of {N_total_possible_transitions} transitions we could use")
        print(f"We are using up to {len(all_transitions)} of them ({N_total_possible_transitions - len(all_transitions)} removed)")

        strong_path = transition_kwds.get("strong_path", None)
        if strong_path is None:
            strong_transitions = None
        else:
            strong_transitions = Transitions.read(strong_path, format="vald")
    
        lambda_min, lambda_max, lambda_step = lambdas = transition_kwds["lambdas"]

        for method, options in config["methods"].items():
            
            if star_description == "HD49933" and method == "moog":
                options["n_chunks"] = 2

            if transition_kwds["lambdas"][0] == 15_000 and method == "moog":
                # 93505 transitions. MOOG will die
                options["n_chunks"] = 10

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
                N_transitions = len(transitions)
                
            elif method == "korg" and transition_kwds["lambdas"][0] != 15_000:
                # If it's the 15_000 region we need to load the TS line list and write to Vald
                # for all others, just supply the vald line list
                transitions = transition_kwds["paths"]
                if len(transitions) == 1:
                    transitions = transitions[0]
                N_transitions = -1
                
            elif method == "turbospectrum" and transition_kwds["lambdas"][0] == 15_000:
                transitions = transition_kwds["paths"]
                N_transitions = -2

            else:
                transitions = Transitions(all_transitions)
                N_transitions = len(transitions)


