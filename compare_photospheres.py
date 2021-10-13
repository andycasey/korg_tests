import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from glob import glob
from grok import Photosphere, PhotosphereInterpolator

try:
    photospheres
except NameError:
        
    photospheres = []
    for path in tqdm(glob("data/photospheres/castelli-kurucz-2004/*.dat")):
        photospheres.append(Photosphere.read(path))

pi = PhotosphereInterpolator(
    photospheres, 
    interpolate_log_quantities=("P", "XNE", "Pg", "Pe", ),
)
point = dict(teff=6635, logg=4.2, m_h=-0.46, alpha_m=0.4)
neighbours = 8
photosphere, common_opacity_scale, quantities, indices, interpolated_column_names \
    = pi(**point, full_output=True, neighbours=neighbours)

columns = ("T", "P", "XNE", "ABROSS")

K = int(np.ceil(np.sqrt(len(columns))))

fig, axes = plt.subplots(K, K)
for ax, column_name in zip(axes.flat, columns):
    ax.plot(
        photosphere["RHOX"],
        photosphere[column_name],
        c="k",
        zorder=10
    )

    #for index in pi.nearest_neighbour_indices([point[k] for k in ("teff", "logg", "m_h", "alpha_m")], neighbours):
    for j, index in enumerate(indices):
        ax.plot(
            photospheres[index]["RHOX"],
            photospheres[index][column_name],
            c="tab:blue",
            ls=":"
        )
        ax.plot(
            photospheres[index]["RHOX"],
            quantities[j, interpolated_column_names.index(column_name)],
            c="tab:red"
        )

    ax.set_ylabel(column_name)
    ax.set_xlabel("RHOX")
fig.tight_layout()
fig.savefig("tmp.png")



raise a


marcs = Photosphere.read("data/photospheres/marcs_mod/s4250_g+1.5_m1.0_t02_st_z-0.50_a+0.20_c+0.00_n+0.00_o+0.20_r+0.00_s+0.00.mod.gz")
castk = Photosphere.read("data/photospheres/castelli-kurucz-2004/am05t4250g15k2odfnew.dat")


fig, axes = plt.subplots(2, 2, figsize=(12, 12))
axes = axes.flat

axes[0].plot(
    marcs["RHOX"],
    marcs["Pg"],
    label="marcs"
)
axes[0].plot(
    castk["RHOX"],
    castk["P"],
    label="CK",
    ls=":"
)
axes[0].set_xlabel("RHOX")
axes[0].set_ylabel("P")
axes[0].loglog()
axes[0].legend()


axes[1].plot(
    marcs["RHOX"],
    marcs["T"],
    label="marcs"
)
axes[1].plot(
    castk["RHOX"],
    castk["T"],
    label="CK",
    ls=":"
)
axes[1].set_xlabel("RHOX")
axes[1].set_ylabel("T")

axes[1].loglog()

axes[2].plot(
    marcs["RHOX"],
    marcs["KappaRoss"],
    label="marcs"
)
axes[2].plot(
    castk["RHOX"],
    castk["ABROSS"],
    label="CK",
    ls=":"
)
axes[2].set_xlabel("RHOX")
axes[2].set_ylabel("ABROSS (ck) / KappaRoss (marcs)")

axes[2].loglog()
axes[2].legend()

axes[3].plot(
    marcs["RHOX"],
    marcs["XNE"],
    label="marcs"
)
axes[3].plot(
    castk["RHOX"],
    castk["XNE"],
    label="CK",
    ls=":"
)
axes[3].set_xlabel("RHOX")
axes[3].set_ylabel("XNE")

axes[3].loglog()
axes[3].legend()
fig.tight_layout()
fig.savefig("comp1.png")