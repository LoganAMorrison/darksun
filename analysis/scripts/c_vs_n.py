"""
script for generating plot of delta relic density as a function of c and n
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from utils import remove_nans, OMEGA_CDM_H2
from scipy.interpolate import interp2d, interp1d
from skimage.measure import find_contours


def find_cont(xs, ys, zs):
    contour = find_contours(zs, OMEGA_CDM_H2)[0]

    jjs, iis = contour.T

    new_ys = np.interp(jjs, np.arange(len(ys)), ys)
    new_xs = np.interp(iis, np.arange(len(xs)), xs)

    return new_xs, new_ys


if __name__ == "__main__":
    # Load the benchmark data and extract needed quantities
    data6 = pd.read_csv(
        "/home/logan/Projects/DarkSun/cpp/rundata/c_vs_n_lam=1e-4_lec1=0.1_lec2=1.0.csv"
    )

    data6.sort_values(by=["C", "N"], inplace=True)
    data6.drop(
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
    print(data6)

    ns = np.unique(np.array(data6["N"]))
    cs = np.unique(np.array(data6["C"]))

    n_array = np.array(data6["N"]).reshape((len(cs), len(ns)))
    c_array = np.array(data6["C"]).reshape((len(cs), len(ns)))

    mpl.use("TkAgg")
    plt.figure(dpi=100)

    rdd_array = np.array(data6["RD_DEL"]).reshape((len(cs), len(ns)))
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
        colors=["firebrick", "steelblue"],
        alpha=0.5,
    )

    plt.ylabel(r"$c$", fontdict={"size": 16})
    plt.xlabel(r"$N$", fontdict={"size": 16})

    NS, LOGCS = find_cont(ns, np.log10(cs), rdd_array)
    logc_interp = interp1d(NS, LOGCS)
    print(10.0 ** logc_interp(6))
    print(10.0 ** logc_interp(7))
    print(10.0 ** logc_interp(8))

    plt.yscale("log")
    plt.xscale("log")
    plt.xlim([np.min(ns), 12])

    plt.xticks(
        [5, 6, 7, 8, 9, 10, 11, 12],
        ["5", "6", "7", "8", "9", "10", "12", "13"],
    )
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

    plt.savefig("/home/logan/Projects/DarkSun/cpp/analysis/figures/c_vs_n.pdf")
