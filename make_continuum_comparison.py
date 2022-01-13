
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
from grok.transitions import Transition
from grok.transitions.utils import (air_to_vacuum, vacuum_to_air)
from grok.utils import read_benchmark_star
from astropy import units as u

OVERWRITE = False

suffix = "-v0.7.0"
# Combinations of `star`, `transitions`, and `methods`.
with open("comparisons.yml", "r") as fp:
    config = yaml.load(fp, Loader=yaml.FullLoader)

lambdas = (lambda_min, lambda_max, lambda_delta) = (2_000, 10_000, 1.0)

output_prefix = lambda star_description, lambda_min, lambda_max, method_description: f"continuum_{star_description}_{lambda_min}_{lambda_max}_{method_description}"

for method_description, options in config["methods"].items():
    
    method, *desc = method_description.split("_")
    if len(desc) > 0:
        continue
    
    for star_description, star in config["stars"].items():

        basename = output_prefix(star_description, lambda_min, lambda_max, method_description)
        output_dir = f"executions/{basename}{suffix}"
        output_path = f"executions/{basename}{suffix}.out"
        if os.path.exists(output_path) and not OVERWRITE:
            print(f"Skipping {method_description} {star_description} because {output_path} exists")
            continue

        print(f"Creating photosphere..")
        photospheres = []
        for i, path in enumerate(tqdm(glob(star["model_kwargs"]["photosphere_grid_wildmask"]), desc="Loading photospheres")):
            photospheres.append(
                Photosphere.read(
                    path, 
                    format=star["model_kwargs"].get("photosphere_read_format", None)
                )
            )

        photosphere_interpolator = NewPhotosphereInterpolator(photospheres)

        photosphere = photosphere_interpolator(**star["model_kwargs"]["photosphere_point"])

        print(f"Executing {method} ({method_description}) for {star_description} between {lambda_min} and {lambda_max}")

        transitions = Transitions.read("data/transitions/fake-line.vald")
        
        if method_description == "moog":
                
            # If we are running MOOG, then we need a transition at least every 10 A to avoid MOOG dying.
            # (We can't get around this by just supplying a single strong line, or something like that.)
            data = dict([(k.lstrip("_"), v) for k, v in transitions[0].__dict__.items()])
            data["lambda_vacuum"] = None

            t = []
            step = 10
            for lambda_air in np.arange(lambda_min, lambda_max + step, step):
                data["lambda_air"] = lambda_air * u.Angstrom
                t.append(Transition(**data))
        
            transitions = Transitions(t)

        try:
            spectrum, meta = synthesize(
                photosphere,
                transitions,
                method=method,
                lambdas=lambdas,
                options=None
            )

        except:
            log.exception("Exception occurred")
            raise
            continue

        print(f"  That took {meta['wallclock_time']} seconds")

        if method_description == "korg":
            wavelength_vacuum = np.arange(meta["lambda_vacuum_min"], meta["lambda_vacuum_max"], meta["lambda_vacuum_delta"])
            continuum = spectrum["continuum"]

            wavelength_air = vacuum_to_air(wavelength_vacuum * u.Angstrom).to("Angstrom").value
            
            assert wavelength_air.shape == continuum.shape
            spectrum = np.vstack([wavelength_air, continuum])
        elif method_description == "turbospectrum":
            wavelength = np.arange(lambda_min, lambda_max + lambda_delta, lambda_delta)
            continuum = np.loadtxt(os.path.join(meta["dir"], "result"), usecols=(2,))
            spectrum = np.vstack([wavelength, continuum])
        elif method_description == "moog":
            wavelength = meta["continuum_lambda_air"]
            continuum = meta["continuum_flux"]
            spectrum = np.vstack([wavelength, continuum])

        else:
            raise NotImplementedError(f"Method {method_description} not implemented")

        np.savetxt(output_path, spectrum)
        
        copy_tree(meta["dir"], output_dir)
        
        print(f"  Wrote to {output_path}")
        print(f"  Copied execution contents to {output_dir}")



# Make some plots.
methods = ("korg", "turbospectrum", "moog")
S = 4 # stars
fig, axes = plt.subplots(4, 1)

plot_kwargs = {
    "korg": dict(c="tab:red"),
    "turbospectrum": dict(c="k"),
    "moog": dict(c="tab:blue")
}
scales = {
    "turbospectrum": 10**8,
    "moog": np.pi
}
for i, ax in enumerate(axes):

    star_description, star = list(config["stars"].items())[i]

    for j, method in enumerate(methods):
        
        scale = scales.get(method, 1)
        basename = output_prefix(star_description, lambda_min, lambda_max, methods[j])
        output_path = f"executions/{basename}{suffix}.out"

        spectrum = np.loadtxt(output_path)
        print(method, star_description, np.log10(np.median(spectrum[1])), scale)

        ax.plot(spectrum[0], scale * spectrum[1], label=method.title(), **plot_kwargs[method])

    ax.set_xlim(lambda_min, lambda_max)

    if i == 3:
        ax.legend(loc="lower right", frameon=False)
    
    if ax.is_last_row():
        ax.set_xlabel(r"$\lambda$")
    else:
        ax.set_xticklabels([])
    ax.set_ylabel(star_description)
    #ax.semilogy()

fig.tight_layout()
fig.savefig(f"continuum{suffix}.png", dpi=300)
