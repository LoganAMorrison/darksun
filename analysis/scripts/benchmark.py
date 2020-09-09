#!/usr/bin/python

"""
Script for generating the benchmark plots for the Dark SU(N) model.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scipy.ndimage import gaussian_filter
from scipy.interpolate import interp2d, interp1d
from utils import remove_nans
from scipy.interpolate import RectBivariateSpline
from matplotlib.lines import Line2D
from skimage.measure import find_contours

OMEGA_CDM_H2 = 0.1198
SI_BOUND = 457.281  # 0.1 cm^2 / g


def self_interaction_eta(n, lam, lec1):
    return 62.012553360599640 * lec1 ** 2 / (lam ** 2 * n ** 5)


def self_interaction_del(lam):
    return 4.0 * np.pi ** 3 / lam ** 2


def find_cont(xs, ys, zs, level=SI_BOUND):
    contour = find_contours(zs, level)[0]

    jjs, iis = contour.T

    new_ys = np.interp(jjs, np.arange(len(ys)), ys)
    new_xs = np.interp(iis, np.arange(len(xs)), xs)

    return new_xs, new_ys


if __name__ == "__main__":

    # Load the benchmark data and extract needed quantities
    # data = pd.read_csv("../../rundata/bm_lec1=0.1_lec2=1.0_xi_inf=1e-3.csv")
    # data = pd.read_csv("../../rundata/bm_lec1=1.0_lec2=0.0.csv")
    data = pd.read_csv("../../rundata/bm_lec1=0.1_lec2=1.0_xi_inf=5e-2.csv")

    data.sort_values(by=["LAM", "N"], inplace=True)

    ns = np.unique(np.array(data["N"]))
    lams = np.unique(np.array(data["LAM"]))

    NS = np.linspace(np.min(ns), 20, 500)
    LAMS = np.logspace(np.log10(np.min(lams)), np.log10(np.max(lams)), 500)

    n_array = np.array(data["N"]).reshape((len(lams), len(ns)))
    lam_array = np.array(data["LAM"]).reshape((len(lams), len(ns)))

    # Relic densities
    rde_array = np.array(data["RD_ETA"]).reshape((len(lams), len(ns)))
    rdd_array = np.array(data["RD_DEL"]).reshape((len(lams), len(ns)))
    remove_nans(rde_array)
    remove_nans(rdd_array)
    rde_interp = RectBivariateSpline(ns, lams, rde_array.T, s=0)
    rdd_interp = RectBivariateSpline(ns, lams, rdd_array.T, s=10.5)
    # rdt_interp = RectBivariateSpline(ns, lams, rde_array + rdd_array)
    # rde_array = gaussian_filter(rde_array, sigma=(10, 0.5))
    # rdd_array = gaussian_filter(rdd_array, sigma=(2, 0.5))
    rdt_array = rde_array + rdd_array  # total

    # Self interactions constraints
    sie_array = np.array(data["ETA_SI_PER_MASS"]).reshape((len(lams), len(ns)))
    sid_array = np.array(data["DEL_SI_PER_MASS"]).reshape((len(lams), len(ns)))
    remove_nans(sie_array)
    remove_nans(sid_array)
    sia_array = np.where(rde_array > rdd_array, sie_array, sid_array)
    print(rde_interp(1, 1))

    # Get the SI contours
    del_si_cont = find_cont(ns, lams, sid_array)
    eta_si_cont = find_cont(ns, lams, sie_array)
    si_del = interp1d(del_si_cont[0], del_si_cont[1])
    si_eta = interp1d(eta_si_cont[0], eta_si_cont[1])

    del_si_cont = find_cont(ns, lams, sid_array, level=SI_BOUND / 10.0 ** 1.5)
    eta_si_cont = find_cont(ns, lams, sie_array, level=SI_BOUND / 10.0 ** 1.5)
    si2_del = interp1d(del_si_cont[0], del_si_cont[1])
    si2_eta = interp1d(eta_si_cont[0], eta_si_cont[1])

    plt.figure(dpi=100)
    plt.contour(
        n_array,
        lam_array,
        rde_array,
        levels=[OMEGA_CDM_H2 / 3, OMEGA_CDM_H2, 3 * OMEGA_CDM_H2],
        colors=["mediumorchid"],
        linewidths=[1, 2, 1],
        linestyles=["--", "-", "-."],
    )
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

    # Delta self interaction constraint
    ns = np.linspace(5, 6, 100)
    sis = si_del(ns)
    si2s = si2_del(ns)
    plt.fill_between(ns, sis, color="steelblue", alpha=0.6)
    plt.fill_between(ns, si2s, sis, color="Peru", alpha=0.6)

    ns = np.linspace(6, 8, 100)
    sis = si_del(ns)
    si2s = si2_del(ns)
    plt.fill_between(ns, sis, color="steelblue", alpha=0.2)
    plt.fill_between(ns, si2s, sis, color="Peru", alpha=0.2)

    # Eta self interaction constraint
    ns = np.linspace(8, 30, 100)
    sis = si_eta(ns)
    si2s = si2_eta(ns)
    plt.fill_between(ns, sis, color="steelblue", alpha=0.6)
    plt.fill_between(ns, si2s, sis, color="Peru", alpha=0.6)

    ns = np.linspace(7, 8, 100)
    sis = si_eta(ns)
    si2s = si2_eta(ns)
    plt.fill_between(ns, sis, color="steelblue", alpha=0.2)
    plt.fill_between(ns, si2s, sis, color="Peru", alpha=0.2)

    plt.ylabel(r"$\Lambda \ (\mathrm{GeV})$", fontsize=16)
    plt.xlabel(r"$N$", fontsize=16)

    # Construct a custom legend
    plt.title(r"$\lambda_1=0.1, \lambda_2 = 1$", fontsize=16)
    lines = [
        Line2D([0], [0], color="k", linewidth=1, linestyle="-"),
        Line2D([0], [0], color="k", linewidth=1, linestyle="-."),
        Line2D([0], [0], color="k", linewidth=1, linestyle="--"),
        Line2D([0], [0], color="mediumorchid", linewidth=1, linestyle="-"),
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
        r"$\bar{\eta}'$",
        r"$\Delta$",
        r"$\sigma_{\mathrm{S.I.}} / m > 0.1 \mathrm{cm}^2/\mathrm{g}$",
        r"$\sigma_{\mathrm{S.I.}} / m > 10^{-2.5} \mathrm{cm}^2/\mathrm{g}$",
    ]
    plt.legend(lines, labels, frameon=False, fontsize=12)

    plt.yscale("log")
    plt.xscale("log")
    plt.xlim([np.min(ns), 30])
    plt.ylim([1e-7, 10])

    plt.xticks(
        [5, 6, 7, 8, 9, 10, 20, 30],
        ["5", "6", "7", "8", "9", "10", "20", "30"],
    )

    plt.tight_layout()
    plt.savefig(
        "/home/logan/Research/DarkSun/cpp/analysis/figures/lec1=0.1_lec2=1_xi_inf=5e-2.pdf"
    )

