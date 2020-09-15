"""
script for generating plot of delta relic density as a function of c and n
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from utils import remove_nans, OMEGA_CDM_H2, find_contour
from scipy.interpolate import interp2d, interp1d


BASE_PATH = "/home/logan/Research/DarkSun/cpp/rundata/"
DATA_FILES = [BASE_PATH + "c_vs_n1.csv", BASE_PATH + "c_vs_n2.csv"]


def generate_plot(idx):
    # Load the benchmark data and extract needed quantities
    data = pd.read_csv(DATA_FILES[idx])

    data.sort_values(by=["C", "N"], inplace=True)
    data.drop(
        [
            "ADEL",
            "DNEFF_CMB",
            "DNEFF_BBN",
            "ETA_SI_PER_MASS",
            "DEL_SI_PER_MASS",
            "XI_FO",
            "TSM_FO",
            "LEC1",
            "LEC2",
            "MU_DEL",
            "MU_ETA",
            "XI_INF",
            "XI_CMB",
            "XI_BBN",
            "RD_ETA",
        ],
        axis=1,
        inplace=True,
    )

    ns = np.unique(np.array(data["N"]))
    cs = np.unique(np.array(data["C"]))

    n_array = np.array(data["N"]).reshape((len(cs), len(ns)))
    c_array = np.array(data["C"]).reshape((len(cs), len(ns)))

    mpl.use("TkAgg")
    plt.figure(dpi=100)

    rdd_array = np.array(data["RD_DEL"]).reshape((len(cs), len(ns)))
    remove_nans(rdd_array)

    plt.contour(
        n_array, c_array, rdd_array, levels=[OMEGA_CDM_H2], colors="k",
    )
    plt.contour(
        n_array,
        c_array,
        rdd_array,
        levels=[OMEGA_CDM_H2 / 3],
        colors="firebrick",
    )
    plt.contour(
        n_array,
        c_array,
        rdd_array,
        levels=[OMEGA_CDM_H2 * 3],
        colors="steelblue",
    )

    plt.contourf(
        n_array,
        c_array,
        rdd_array,
        levels=[OMEGA_CDM_H2 / 3.0, OMEGA_CDM_H2, 3.0 * OMEGA_CDM_H2],
        colors=["firebrick", "k", "steelblue"],
        alpha=0.5,
    )

    plt.ylabel(r"$c$", fontdict={"size": 16})
    plt.xlabel(r"$N$", fontdict={"size": 16})

    NS, LOGCS = find_contour(ns, np.log10(cs), rdd_array, OMEGA_CDM_H2)
    logc_interp = interp1d(NS, LOGCS)
    print(10.0 ** logc_interp(6))
    print(10.0 ** logc_interp(7))
    print(10.0 ** logc_interp(8))

    plt.yscale("log")
    plt.xscale("log")
    plt.xlim([np.min(ns), 15])

    xticks = np.arange(5, 16)
    plt.xticks(xticks, [str(t) for t in xticks])

    lines = [
        Line2D([0], [0], color="k", linewidth=1, linestyle="-"),
        Line2D([0], [0], color="steelblue", linewidth=1, linestyle="-"),
        Line2D([0], [0], color="firebrick", linewidth=1, linestyle="-"),
    ]
    labels = [
        r"$\Omega_{\Delta} h^2 = \Omega_{\mathrm{DM}} h^2$",
        r"$\Omega_{\Delta} h^2 = 3\Omega_{\mathrm{DM}} h^2$",
        r"$\Omega_{\Delta} h^2 = \frac{1}{3}\Omega_{\mathrm{DM}} h^2$",
    ]
    plt.legend(lines, labels, frameon=False, fontsize=12)

    plt.savefig(
        "/home/logan/Research/DarkSun/cpp/analysis/figures/c_vs_n"
        + str(idx + 1)
        + ".pdf"
    )

    if idx == 1:
        plt.show()


if __name__ == "__main__":
    for idx in range(len(DATA_FILES)):
        generate_plot(idx)

