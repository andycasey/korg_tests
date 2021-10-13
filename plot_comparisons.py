import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pickle
from matplotlib import gridspec
from grok.utils import read_benchmark_star
from matplotlib.patches import ConnectionPatch
from scipy.ndimage import gaussian_filter

output_prefix = lambda star_description, lambda_min, lambda_max, method: f"{star_description}_{lambda_min:.0f}_{lambda_max:.0f}_{method}_interpolated"
figure_path = lambda star_description, **_: f"compare-{star_description}-interpolated.png"

stars = [
    # Arcturus
    dict(
        description="Arcturus",
        observation_path="data/ATLAS.Arcturus_47000.fits",
        model_kwargs=dict(
            # NOTE: Not using interpolation.
            #photosphere_path = "data/photospheres/marcs_mod/s4250_g+1.5_m1.0_t02_st_z-0.50_a+0.20_c+0.00_n+0.00_o+0.20_r+0.00_s+0.00.mod.gz",
            photosphere_path = "data/photospheres/castelli-kurucz-2004/am05t4250g15k2odfnew.dat",
            #photosphere_grid_wildmask="data/photospheres/marcs_mod/s*_m1.0_t02_st_*a+0.20*",
            photosphere_grid_wildmask="data/photospheres/castelli-kurucz-2004/*.dat",
            photosphere_read_format="atlas9",
            photosphere_point=dict(
                teff=4286,
                logg=1.64,
                m_h=-0.53, 
                alpha_m=0
            )
        )
    ),
    # Sun.
    dict(
        description="Sun",
        observation_path="data/ATLAS.Sun_47000.fits",
        model_kwargs=dict(  
            # NOTE: Not using interpolation.
            photosphere_path="data/photospheres/castelli-kurucz-2004/ap00t5750g45k2odfnew.dat",
            
            photosphere_grid_wildmask="data/photospheres/castelli-kurucz-2004/*.dat",
            photosphere_read_format="atlas9",
            photosphere_point=dict(
                teff=5777,
                logg=4.4,
                m_h=0,
                alpha_m=0
            )
        )
    ),   
    # HD49933
    dict(
        description="HD49933",
        observation_path="data/HARPS.Archive_HD49933_47000.fits",
        model_kwargs=dict(
            photosphere_path="data/photospheres/castelli-kurucz-2004/am05t6500g40k2odfnew.dat",
            photosphere_grid_wildmask="data/photospheres/castelli-kurucz-2004/*.dat",
            photosphere_read_format="atlas9",
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
            #photosphere_grid_wildmask="data/photospheres/marcs_mod/s*_m1.0_t02_st_*a+0.20*",
            photosphere_path="data/photospheres/marcs_mod/s4500_g+1.5_m1.0_t02_st_z-2.50_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod.gz",
            photosphere_grid_wildmask="data/photospheres/castelli-kurucz-2004/*.dat",
            photosphere_read_format="atlas9",
            photosphere_point=dict(
                teff=4587,
                logg=1.61,
                m_h=-2.74, 
                alpha_m=+0.4
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
    "korg": dict(),
}

from scipy.stats.distributions import norm as normal

def convolve_flux(wavelengths, fluxes, R):

    sigma = np.mean(wavelengths / R / 2)
    return gaussian_filter(fluxes, sigma)



total_times = { method: [] for method in methods.keys() }


plot_kwargs = {
    "moog": dict(c="tab:red"),
    "turbospectrum": dict(c="tab:blue", ls=":"),
    "korg": dict(c="#666666", zorder=5),
    "observation": dict(c="k", zorder=10)
}


T = len(all_transition_kwds)
ylim = (0, 1.1)
diff_ylim = (-0.10, 0.10)



for star in stars:
    
    star_description = star["description"]

    obs_wl, obs_flux, obs_sigma = read_benchmark_star(star["observation_path"])
    bad = (obs_flux <= 0)
    obs_flux[bad] = np.interp(obs_wl[bad], obs_flux[~bad], obs_flux[~bad])
    #obs_flux[obs_flux <= 0] = np.nan
    fig = plt.figure(constrained_layout=True, figsize=(30, 3))

    lambda_ptp = np.array([np.ptp(ea["lambdas"][:2]) for ea in all_transition_kwds])
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

    for t, (ax, diff_ax, transition_kwds) in enumerate(zip(axes, diff_axes, all_transition_kwds)):

        lambda_min, lambda_max, _ = transition_kwds["lambdas"]

        si, ei = obs_wl.searchsorted([lambda_min, lambda_max])

        ax.plot(
            obs_wl[si:ei],
            obs_flux[si:ei],
            label=star_description,
            **plot_kwargs["observation"],
        )

        for method, _ in methods.items(): 
            output_path = f"executions/{output_prefix(star_description, lambda_min, lambda_max, method)}.pkl"
            if method == "korg":
                output_path = output_path[:-4]
                if output_path.endswith("_interpolated"):
                    output_path = output_path[:-13] # to remove the _interpolated

            print(f"LOading {output_path}")
            if not os.path.exists(output_path):
                continue
            
            if method == "korg":
                flux = np.loadtxt(output_path)
                flux /= np.max(flux)
                wl = np.arange(lambda_min, lambda_max + 0.01, 0.01)[:len(flux)]
                meta = dict(wallclock_time=0)

                convolved_flux = convolve_flux(wl, flux, star.get("R", 15000))
                convolved_flux = gaussian_filter(flux, 2)
                #    return gaussian_filter(fluxes, sigma)


            else:
                with open(output_path, "rb") as fp:
                    spectrum, meta = pickle.load(fp)

                wl, flux = (np.array(spectrum["wavelength"]), np.array(spectrum["rectified_flux"]))
                if len(wl) == 0:
                    print(f"NO SPECTRUM IN {output_path}")
                    continue

                convolved_flux = convolve_flux(wl, flux, star.get("R", 25_000))

            total_times[method].append(meta["wallclock_time"])

            ax.plot(
                wl,
                convolved_flux,
                label=method,
                **plot_kwargs[method],
            )

            interpolated_flux = np.interp(
                obs_wl[si:ei],
                wl,
                convolved_flux
            )
            diff_ax.plot(
                obs_wl[si:ei],
                interpolated_flux - obs_flux[si:ei],
                **plot_kwargs[method]
            )

        diff_ax.set_xlim(lambda_min, lambda_max)
        diff_ax.set_ylim(diff_ylim)
        ax.set_xlim(lambda_min, lambda_max)
        ax.set_ylim(ylim)

    fig.text(0.5, 0.03, 'Wavelength', ha='center')
    axes[-1].legend(frameon=False, loc="lower right")

    plt.draw()

    # Deal with spines, labels, and common limits.
    for i, (ax, diff_ax, transition_kwds) in enumerate(zip(axes, diff_axes, all_transition_kwds)):
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
            next_lambda_min = all_transition_kwds[i + 1]["lambdas"][0]

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
            kwargs = dict(transform=ax.transAxes, color='k', lw=plt.rcParams["axes.linewidth"], clip_on=False)
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



fig, ax = plt.subplots()

yticks = list(range(len(total_times)))

ax.barh(
    y=yticks,
    width=[np.sum(v) for v in total_times.values()],
)
ax.set_yticks(yticks)
ax.set_yticklabels(list(total_times.keys()))
ax.yaxis.set_tick_params(width=0)
ax.set_xlabel("Total execution time / s")
fig.tight_layout()
fig.savefig("a.png")