#!/usr/bin/python

"""
Script for generating the benchmark plots for the Dark SU(N) model.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scipy.ndimage import gaussian_filter
from scipy.interpolate import interp2d, interp1d
from scipy.interpolate import RectBivariateSpline
from matplotlib.lines import Line2D
from scipy.optimize import curve_fit

from utils import (
    remove_nans,
    OMEGA_CDM_H2,
    SI_BOUND,
    find_contour,
    NEFF_BBN_BOUND,
    NEFF_CMB_BOUND,
)

FILE_NAMES = [
    "bm_lec1=1_lec2=0_xi_inf=1e-2",
    "bm_lec1=0.1_lec2=1_xi_inf=1e-2",
    "bm_lec1=0.1_lec2=1_xi_inf=5e-2",
    "bm_lec1=1e-3_lec2=1_xi_inf=1e-2",
]

DATA_FILES = [
    os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..",
        "..",
        "rundata",
        fn + ".csv",
    )
    for fn in FILE_NAMES
]

Y_LOW_BOUNDS = [1e-4, 1e-4, 1e-7, 1e-4]
# TITLES = [
#    r"$\lambda_1=1, \ \lambda_2=0, \ \xi_{\infty}=10^{-2}$",
#    r"$\lambda_1=0.1, \ \lambda_2=1, \ \xi_{\infty}=10^{-2}$",
#    r"$\lambda_1=0.1, \ \lambda_2=1, \ \xi_{\infty}=5\times10^{-2}$",
#    r"$\lambda_1=10^{-3}, \ \lambda_2=1, \ \xi_{\infty}=10^{-2}$",
# ]
TITLES = ["Minimal", "Decorrelated", "Hot", "Low SI"]


def generate_bm_plot(idx):
    # Load the benchmark data and extract needed quantities
    data = pd.read_csv(DATA_FILES[idx])

    data.sort_values(by=["LAM", "N"], inplace=True)

    ns = np.unique(np.array(data["N"]))
    lams = np.unique(np.array(data["LAM"]))

    n_array = np.array(data["N"]).reshape((len(lams), len(ns)))
    lam_array = np.array(data["LAM"]).reshape((len(lams), len(ns)))

    # Relic densities
    rde_array = np.array(data["RD_ETA"]).reshape((len(lams), len(ns)))
    rdd_array = np.array(data["RD_DEL"]).reshape((len(lams), len(ns)))
    remove_nans(rde_array)
    remove_nans(rdd_array)
    rdt_array = rde_array + rdd_array  # total

    rd_cont = np.log10(find_contour(ns, lams, rde_array, OMEGA_CDM_H2))
    print(TITLES[idx])
    print(curve_fit(lambda x, m, b: m * x + b, rd_cont[0], rd_cont[1]))

    # Self interactions constraints
    sie_array = np.array(data["ETA_SI_PER_MASS"]).reshape((len(lams), len(ns)))
    sid_array = np.array(data["DEL_SI_PER_MASS"]).reshape((len(lams), len(ns)))
    remove_nans(sie_array)
    remove_nans(sid_array)

    # BBN + CMB constraints
    bbn_array = np.array(data["DNEFF_BBN"]).reshape((len(lams), len(ns)))
    cmb_array = np.array(data["DNEFF_CMB"]).reshape((len(lams), len(ns)))
    remove_nans(bbn_array)
    remove_nans(cmb_array)

    # Get the SI contours
    del_si_cont = find_contour(ns, lams, sid_array, SI_BOUND)
    eta_si_cont = find_contour(ns, lams, sie_array, SI_BOUND)
    si_del = interp1d(del_si_cont[0], del_si_cont[1])
    si_eta = interp1d(eta_si_cont[0], eta_si_cont[1])

    del_si_cont = find_contour(
        ns, lams, sid_array, level=SI_BOUND / 10.0 ** 1.5
    )
    eta_si_cont = find_contour(
        ns, lams, sie_array, level=SI_BOUND / 10.0 ** 1.5
    )
    si2_del = interp1d(del_si_cont[0], del_si_cont[1])
    si2_eta = interp1d(eta_si_cont[0], eta_si_cont[1])

    plt.figure(dpi=100)
    #    plt.contour(
    #        n_array,
    #        lam_array,
    #        rde_array,
    #        levels=[OMEGA_CDM_H2 / 3, OMEGA_CDM_H2, 3 * OMEGA_CDM_H2],
    #        colors=["mediumorchid"],
    #        linewidths=[1, 2, 1],
    #        linestyles=["--", "-", "-."],
    #    )
    # Delta Contours
    plt.contour(
        n_array,
        lam_array,
        rdd_array,
        levels=[OMEGA_CDM_H2 / 3, OMEGA_CDM_H2, OMEGA_CDM_H2 * 3],
        colors=["firebrick"],
        linewidths=[1, 2, 1],
        linestyles=["--", "-", "-."],
    )
    plt.contour(
        n_array,
        lam_array,
        rdt_array,
        levels=[OMEGA_CDM_H2 / 3, OMEGA_CDM_H2, 3 * OMEGA_CDM_H2],
        colors=["k"],
        linewidths=[1, 2, 1],
        linestyles=["--", "-", "-."],
    )
    plt.contour(
        n_array,
        lam_array,
        rdt_array,
        levels=[OMEGA_CDM_H2 / 3, 3 * OMEGA_CDM_H2],
        colors=["k"],
        linewidths=1,
        linestyles=["--", "-."],
    )

    # BBN + CMB contours
    plt.contourf(
        n_array,
        lam_array,
        bbn_array,
        levels=[NEFF_BBN_BOUND[0] + NEFF_BBN_BOUND[1], 1e20],
        colors=["goldenrod"],
    )
    plt.contourf(
        n_array,
        lam_array,
        cmb_array,
        levels=[NEFF_CMB_BOUND[0] + NEFF_CMB_BOUND[1], 1e20],
        colors=["teal"],
    )

    # Delta self interaction constraint
    ns = np.linspace(5, 6, 100)
    sis = si_del(ns)
    si2s = si2_del(ns)
    plt.fill_between(ns, sis, color="steelblue", alpha=0.6)
    plt.plot(ns, si2s, color="Peru", alpha=0.8, ls="--", lw=3)

    ns = np.linspace(6, 8, 100)
    sis = si_del(ns)
    si2s = si2_del(ns)
    # plt.fill_between(ns, sis, color="steelblue", alpha=0.2)
    plt.plot(ns, si2s, color="Peru", alpha=0.5, ls="--", lw=3)

    ns = np.linspace(6, 7.5, 100)
    y1, n1 = si_del(6), 6
    y2, n2 = Y_LOW_BOUNDS[idx], 7.5
    slope = (np.log10(y2) - np.log10(y1)) / (np.log10(n2) - np.log10(n1))
    ylow = 10 ** np.array(
        [np.log10(y1) + slope * (np.log10(n) - np.log10(n1)) for n in ns]
    )
    plt.fill_between(ns, sis, ylow, color="steelblue", alpha=0.2)
    plt.fill_between(ns, ylow, color="steelblue", alpha=0.6)

    ns = np.linspace(7.5, 8.0, 100)
    plt.fill_between(ns, sis, color="steelblue", alpha=0.2)

    # Eta self interaction constraint
    ns = np.linspace(8, 30, 100)
    sis = si_eta(ns)
    si2s = si2_eta(ns)
    plt.fill_between(ns, sis, color="steelblue", alpha=0.6)
    plt.plot(ns, si2s, color="Peru", alpha=0.8, ls="--", lw=3)

    ns = np.linspace(7, 8, 100)
    sis = si_eta(ns)
    si2s = si2_eta(ns)
    # plt.fill_between(ns, sis, color="steelblue", alpha=0.2)
    plt.plot(ns, si2s, color="Peru", alpha=0.5, ls="--", lw=3)

    plt.ylabel(r"$\Lambda \ (\mathrm{GeV})$", fontsize=16)
    plt.xlabel(r"$N$", fontsize=16)

    # Construct a custom legend
    lines = [
        Line2D([0], [0], color="k", linewidth=1, linestyle="-"),
        Line2D([0], [0], color="k", linewidth=1, linestyle="-."),
        Line2D([0], [0], color="k", linewidth=1, linestyle="--"),
        # Line2D([0], [0], color="mediumorchid", linewidth=1, linestyle="-"),
        Line2D([0], [0], color="firebrick", linewidth=1, linestyle="-"),
        Line2D(
            [0], [0], color="steelblue", linewidth=3, alpha=0.5, linestyle="-"
        ),
        Line2D([0], [0], color="Peru", linewidth=3, alpha=0.5, linestyle="-"),
    ]
    labels = [
        r"$\Omega_{\bar{\eta}'+\Delta} h^2 = \Omega_{\mathrm{DM}} h^2$",
        r"$\Omega_{\bar{\eta}'+\Delta} h^2 = 3\Omega_{\mathrm{DM}} h^2$",
        r"$\Omega_{\bar{\eta}'+\Delta} h^2 = \frac{1}{3}\Omega_{\mathrm{DM}} h^2$",
        # r"$\bar{\eta}'$",
        r"$\Delta$",
        r"$\sigma_{\mathrm{S.I.}} / m > 0.1 \mathrm{cm}^2/\mathrm{g}$",
        r"$\sigma_{\mathrm{S.I.}} / m > 10^{-2.5} \mathrm{cm}^2/\mathrm{g}$",
    ]
    plt.legend(lines, labels, frameon=False, fontsize=12)

    plt.yscale("log")
    plt.xscale("log")
    plt.xlim([np.min(ns), 30])
    plt.ylim([Y_LOW_BOUNDS[idx], 10])

    plt.xticks(
        [5, 6, 7, 8, 9, 10, 20, 30],
        ["5", "6", "7", "8", "9", "10", "20", "30"],
    )

    plt.title(TITLES[idx], fontsize=16)

    plt.tight_layout()
    plt.savefig(
        "/home/logan/Research/DarkSun/cpp/analysis/figures/"
        + FILE_NAMES[idx]
        + ".pdf"
    )


if __name__ == "__main__":
    for idx in range(len(DATA_FILES)):
        generate_bm_plot(idx)
