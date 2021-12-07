import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pickle
import yaml
from collections import OrderedDict
from matplotlib import gridspec
from grok.utils import read_benchmark_star
from matplotlib.patches import ConnectionPatch
from scipy.ndimage import gaussian_filter
from grok.transitions.utils import (air_to_vacuum, vacuum_to_air)


korg_version = (0, 5, 0) # TODO: Check this from the meta of each file that it matches
korg_major_version, korg_minor_version, korg_micro_version = korg_version

# Specify some paths. The output_prefix should match that from `make_comparisons.py`
output_prefix = lambda star_description, lambda_min, lambda_max, method_description: f"{star_description}_{lambda_min:.0f}_{lambda_max:.0f}_{method_description}"

figure_path = lambda star_description, **_: f"compare-{star_description}-v{korg_major_version}.{korg_minor_version}.{korg_micro_version}.png"
wallclock_figure_path = f"wallclock-v{korg_major_version}.{korg_minor_version}.{korg_micro_version}.png"

# Load in the same file that we used in `make_comparisons.py`
with open("comparisons.yml", "r") as fp:
    config = yaml.load(fp, Loader=yaml.FullLoader)

# Some plotting keywords
plot_kwargs = {
    "moog": dict(c="#CCCCCC", zorder=1, ms=0),
    "turbospectrum": dict(c="#666666", zorder=2, ms=0),
    "korg": dict(c="tab:red", zorder=5, ls=":", ms=0),
    "observation": dict(c="k", zorder=10, ms=0),
    "korg_turbospectrum_diff": dict(c="tab:red", alpha=0.5, ms=0)
}
ylim = (0, 1.1)
diff_ylim = (-0.10, 0.10)
log_abs_diff = False

def convolve_flux(wavelengths, fluxes, R):

    sigma = np.mean(wavelengths / R / 2)
    return gaussian_filter(fluxes, sigma)


total_times = { method: [] for method in config["methods"].keys() }

T = len(config["transitions"])

for star_description, star in config["stars"].items():
    
    obs_wl, obs_flux, obs_sigma = read_benchmark_star(star["observation_path"])
    bad = (obs_flux <= 0)
    obs_flux[bad] = np.interp(obs_wl[bad], obs_flux[~bad], obs_flux[~bad])
    #obs_flux[obs_flux <= 0] = np.nan
    fig = plt.figure(constrained_layout=True, figsize=(30, 3))

    lambda_ptp = np.array([np.ptp(ea["lambdas"][:2]) for ea in config["transitions"]])
    gs = gridspec.GridSpec(
        ncols=T,
        nrows=2,
        figure=fig,
        width_ratios=lambda_ptp / np.min(lambda_ptp),
        height_ratios=[1, 5]
    )

    diff_axes = [fig.add_subplot(gs[0, i]) for i in range(T)]
    axes = [fig.add_subplot(gs[1, i]) for i in range(T)]
    
    for diff_ax in diff_axes:
        diff_ax.set_xticks([])
        if not diff_ax.get_subplotspec().is_first_col():
            diff_ax.set_yticks([])

    for t, (ax, diff_ax, transition_kwds) in enumerate(zip(axes, diff_axes, config["transitions"])):

        lambda_min, lambda_max, _ = transition_kwds["lambdas"]

        si, ei = obs_wl.searchsorted([lambda_min, lambda_max])

        ax.plot(
            obs_wl[si:ei],
            obs_flux[si:ei],
            label=star_description,
            **plot_kwargs["observation"],
        )

        synthesized_spectra = {}
        for method_description, _ in config["methods"].items(): 
            method, *desc = method_description.split("_")
            output_path = f"executions/{output_prefix(star_description, lambda_min, lambda_max, method_description)}.pkl"


            print(f"Loading {output_path}")
            if not os.path.exists(output_path):
                continue
            
            with open(output_path, "rb") as fp:
                spectrum, meta = pickle.load(fp)

            # Check version in meta.
            if method == "korg":
                # TODO: naming here inconsistent: korg_major_version and korg_version_major
                actual_korg_major_version = int(meta.get("korg_version_major", -1))
                actual_korg_minor_version = int(meta.get("korg_version_minor", -1))
                actual_korg_micro_version = int(meta.get("korg_version_micro", -1))
                
                if actual_korg_major_version != korg_major_version \
                or actual_korg_minor_version != korg_minor_version \
                or actual_korg_micro_version != korg_micro_version:
                    print(
                        f"WARNING: We are making a figure for Korg v{korg_major_version}.{korg_minor_version}.{korg_micro_version}, "
                        f"but the metadata in {output_path} suggests this comparison was executed with Korg v{actual_korg_major_version}.{actual_korg_minor_version}.{actual_korg_micro_version}. "
                        f"I won't make a fuss, but you should ensure that this is what you want."
                    )

            wl, flux = (np.array(spectrum["wavelength"]), np.array(spectrum["rectified_flux"]))
            wl = wl[:len(flux)]
            flux = flux[:len(wl)]
            if len(wl) == 0:
                print(f"NO SPECTRUM IN {output_path}")
                continue


            total_times[method_description].append(meta["wallclock_time"])

            if desc and desc[0] == "no-hydrogen-lines":
                print(f"Not plotting {method} ({method_description}) in the spectrum figure, since we only did it for comparison.")
                continue

            #convolved_flux = convolve_flux(wl, flux, star.get("R", 25_000))
            synthesized_spectra[method] = [wl, flux]


            ax.plot(
                wl,
                flux,
                label=method,
                **plot_kwargs[method],
            )

            #interpolated_flux = np.interp(
            #    obs_wl[si:ei],
            #    wl,
            #    flux
            #)
            #diff_ax.plot(
            #    obs_wl[si:ei],
            #    interpolated_flux - obs_flux[si:ei],
            #    **plot_kwargs[method]
            #)

        korg_wl_air, korg_flux = synthesized_spectra["korg"]
        ts_wl_air, ts_flux = synthesized_spectra["turbospectrum"]

        korg_flux_air = np.interp(
            ts_wl_air,
            korg_wl_air,
            korg_flux
        )

        diff_flux = ts_flux - korg_flux_air
        if log_abs_diff:
            diff_flux = np.abs(diff_flux)

        diff_ax.plot(
            ts_wl_air,
            diff_flux,
            **plot_kwargs["korg_turbospectrum_diff"],
        )

        diff_ax.set_xlim(lambda_min, lambda_max)
        diff_ax.set_ylim(diff_ylim)

        if log_abs_diff:
            diff_ax.set_yscale("log")

        ax.set_xlim(lambda_min, lambda_max)
        ax.set_ylim(ylim)

    fig.text(0.5, 0.03, 'Wavelength', ha='center')
    axes[-1].legend(frameon=False, loc="lower right")

    plt.draw()

    # Deal with spines, labels, and common limits.
    for i, (ax, diff_ax, transition_kwds) in enumerate(zip(axes, diff_axes, config["transitions"])):
        lambda_min, lambda_max, _ = transition_kwds["lambdas"]

        if ax.get_subplotspec().is_first_col():
            ax.set_ylabel(r"Rectified flux")
            hide_spines = ("right", )
        else:
            ax.set_yticks([])
            if ax.get_subplotspec().is_last_col():
                hide_spines = ("left", )                
            else:
                hide_spines = ("left", "right")
            
        for spine in hide_spines:
            ax.spines[spine].set_visible(False)
            diff_ax.spines[spine].set_visible(False)
                
        if not ax.get_subplotspec().is_last_col():
            next_ax = axes[i + 1]
            next_lambda_min = config["transitions"][i + 1]["lambdas"][0]

            for value in ylim:
                ax.add_artist(
                    ConnectionPatch(
                        xyA=(lambda_max, value), 
                        xyB=(next_lambda_min, value), 
                        coordsA="data", coordsB="data", 
                        axesA=ax, axesB=next_ax,
                        color="k",
                        lw=plt.rcParams["axes.linewidth"]
                    )
                )

            d = .015  # how big to make the diagonal lines in axes coordinates
            # arguments to pass to plot, just so we don't keep repeating them
            xc = 0.125
            xo = 0.005
            kwargs = dict(ms=0, transform=ax.transAxes, color='k', lw=plt.rcParams["axes.linewidth"], clip_on=False)
            ax.plot((1 + xc + xo - d, 1 + xc + xo + d), (-d, +d), **kwargs)  # top-right diagonal
            ax.plot((1 + xc + xo - d, 1 + xc + xo + d), (1 - d, 1 + d), **kwargs)        # bottom right
            ax.plot((1 + xc - xo - d, 1 + xc - xo + d), (-d, +d), **kwargs)  # top-right diagonal
            ax.plot((1 + xc - xo - d, 1 + xc - xo + d), (1 - d, 1 + d), **kwargs)        # bottom right
            
            next_diff_ax = diff_axes[i + 1]
            for value in diff_ylim:
                diff_ax.add_artist(
                    ConnectionPatch(
                        xyA=(lambda_max, value), 
                        xyB=(next_lambda_min, value), 
                        coordsA="data", coordsB="data", 
                        axesA=diff_ax, axesB=next_diff_ax,
                        color="k",
                        lw=plt.rcParams["axes.linewidth"]
                    )
                )   

            """
            ax_ratio = np.ptp(ax.get_xlim())/np.ptp(ax.get_ylim())
            ax_diff_ratio = 1/(np.ptp(diff_ax.get_xlim())/np.ptp(diff_ax.get_ylim()))
            d = 0.015 # ax_diff_ratio / ax_ratio
            
              # how big to make the diagonal lines in axes coordinates
            # arguments to pass to plot, just so we don't keep repeating them
            xc = 0.125
            xo = 0.005
            kwargs = dict(transform=diff_ax.transAxes, color='k', lw=plt.rcParams["axes.linewidth"], clip_on=False)
            diff_ax.plot((1 + xc + xo - d, 1 + xc + xo + d), (-d, +d), **kwargs)  # top-right diagonal
            diff_ax.plot((1 + xc + xo - d, 1 + xc + xo + d), (1 - d, 1 + d), **kwargs)        # bottom right
            diff_ax.plot((1 + xc - xo - d, 1 + xc - xo + d), (-d, +d), **kwargs)  # top-right diagonal
            diff_ax.plot((1 + xc - xo - d, 1 + xc - xo + d), (1 - d, 1 + d), **kwargs)        # bottom right
            """

        ax.set_xlim(lambda_min, lambda_max)
        ax.set_ylim(ylim)
        ax.set_xlabel(" ")

    path = figure_path(star_description)
    fig.savefig(path, dpi=300)

    print(f"Created figure {path}")



plt.style.use("latex_figure.mplstyle")


korg_version = (0, 4, 1)

# Group by method, as we turn hydrogen lines on/off for Korg and TS.
latex_labels = {
    "moog": r"$\textrm{MOOG}$",
    "korg": r"$\textrm{{Korg}}$",#~({0}.{1}.{2})}}$".format(*korg_version),
    "turbospectrum": r"$\textrm{Turbospectrum}$",#~(19.1.3)}$",
}

# Some fake data.
np.random.seed(0)

N = 4 * 5 # stars * regions
total_times["moog_no-hydrogen-lines"] = total_times.pop("moog")
total_times["moog"] = []


t = { k: np.sum(v) for k, v in total_times.items()}

order = ("korg", "moog", "turbospectrum")

if order is None:
    method_times = OrderedDict(sorted([(k, v) for k, v in t.items() if "_" not in k], key=lambda item: item[1] if item[1] > 0 else np.inf))
else:
    method_times = OrderedDict([(k, t[k]) for k in order])


fig, ax = plt.subplots(figsize=(6, 3))

yticks = list(range(len(method_times)))
ax.barh(
    y=yticks,
    width=method_times.values(),
    color=[("tab:blue" if method == "korg" else "#cccccc") for method in method_times.keys()],
    lw=0,
)
# Do the same for the no hydrogen lines.
ax.barh(
    y=yticks,
    width=[np.sum(total_times.get(f"{method}_no-hydrogen-lines", method_times[method])) for method in method_times.keys()],
    height=0.6,
    color="k",#[("tab:red" if method == "korg" else "#666666") for method in method_times.keys()],
    alpha=0.5,
    zorder=10,
)

ax.barh(y=[0], width=[0], color="#cccccc", label=r"$\textrm{With~hydrogen~lines}$")
ax.barh(y=[0], width=[0], color="#666666", label=r"$\textrm{Without~hydrogen~lines}$")

# Draw arrow for MOOG.
ax.arrow(
    np.sum(total_times["moog_no-hydrogen-lines"]),
    list(method_times.keys()).index("moog"),
    dx=50,
    dy=0,
    length_includes_head=True,
    head_starts_at_zero=False,
    width=0.01,
    head_width=0.1,
    head_length=20,
    color="#666666"
)

ax.legend(frameon=False)

ax.set_yticks(yticks)
ax.set_yticklabels([latex_labels.get(k, k) for k in method_times.keys()])
ax.yaxis.set_tick_params(width=0)
ax.set_xlabel(r"$\textrm{Total~execution~time~for~all~syntheses}$ $\textrm{[seconds]}$")
fig.tight_layout()

fig.savefig(f"time-elapsed-v{korg_major_version}.{korg_minor_version}.{korg_micro_version}.png", dpi=300)


#print(f"Created {wallclock_figure_path}")

order = ("korg", "moog", "turbospectrum")

S = len(config["stars"])
T = len(config["transitions"])

fig, axes = plt.subplots(S, T, figsize=(12, 10))

for s, (ax_row, (star_description, star)) in enumerate(zip(axes, config["stars"].items())):
    for t, (ax, transition_kwds) in enumerate(zip(ax_row, config["transitions"])):

        lambda_min, lambda_max, _ = transition_kwds["lambdas"]

        total_times = {}
        for method_description, _ in config["methods"].items(): 
            method, *desc = method_description.split("_")
            output_path = f"executions/{output_prefix(star_description, lambda_min, lambda_max, method_description)}.pkl"

            if not os.path.exists(output_path):
                total_times[method_description] = 0
                continue

            with open(output_path, "rb") as fp:
                spectrum, meta = pickle.load(fp)
            
            total_times[method_description] = meta["wallclock_time"]

        total_times["moog_no-hydrogen-lines"] = total_times.pop("moog")
        total_times["moog"] = 0

        times = { k: v for k, v in total_times.items()}
        method_times = OrderedDict([(k, times[k]) for k in order])


        yticks = list(range(len(method_times)))
        ax.barh(
            y=yticks,
            width=method_times.values(),
            color=[("tab:blue" if method == "korg" else "#cccccc") for method in method_times.keys()],
            lw=0,
        )
        # Do the same for the no hydrogen lines.
        ax.barh(
            y=yticks,
            width=[np.sum(total_times.get(f"{method}_no-hydrogen-lines", method_times[method])) for method in method_times.keys()],
            height=0.6,
            color="k",#[("tab:red" if method == "korg" else "#666666") for method in method_times.keys()],
            alpha=0.5,
            zorder=10,
        )

        ax.barh(y=[0], width=[0], color="#cccccc", label=r"$\textrm{With~hydrogen~lines}$")
        ax.barh(y=[0], width=[0], color="#666666", label=r"$\textrm{Without~hydrogen~lines}$")

        # Draw arrow for MOOG.
        """
        ax.arrow(
            np.sum(total_times["moog_no-hydrogen-lines"]),
            list(method_times.keys()).index("moog"),
            dx=np.sum(total_times["moog_no-hydrogen-lines"]) * 0.1,
            dy=0,
            length_includes_head=True,
            head_starts_at_zero=False,
            width=0.01,
            head_width=0.1,
            head_length=20,
            color="#666666"
        )
        """

        #ax.legend(frameon=False)

        if ax.is_first_row():
            ax.set_title(f"${lambda_min:.0f}\\textrm{{~to~}}{lambda_max:.0f}\\textrm{{\AA}}$", fontsize=12)
        
        if ax.is_last_col():
            twin_ax = ax.twinx()
            twin_ax.set_yticks([np.mean(ax.get_ylim())])
            twin_ax.yaxis.set_tick_params(width=0)
            twin_ax.set_ylim(ax.get_ylim())
            twin_ax.set_yticklabels([f"$\\textrm{{{star_description}}}$"], rotation=0)
        ax.set_yticks(yticks)
        if ax.is_first_col():
            ax.set_yticklabels([latex_labels.get(k, k) for k in method_times.keys()])
        else:
            ax.set_yticklabels([])
        ax.yaxis.set_tick_params(width=0)

        if ax.is_last_row() and t == 2:
            ax.set_xlabel(r"$\textrm{Execution~time}$ $\textrm{[seconds]}$")

        if not ax.is_last_row():
            ax.set_xticklabels([])


for ax_column in axes.T:
    xlim_max = max([ax.get_xlim()[1] for ax in ax_column])
    for ax in ax_column:
        ax.set_xlim(0, xlim_max)


fig.tight_layout()
fig.savefig(f"time-elapsed-breakdown-v{korg_major_version}.{korg_minor_version}.{korg_micro_version}.png", dpi=300)
