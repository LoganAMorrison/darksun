#!/usr/bin/python

"""
Script for generating the plots of xifo vs. N
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scipy.ndimage import gaussian_filter
from scipy.interpolate import interp2d, interp1d
from scipy.interpolate import RectBivariateSpline
from matplotlib.lines import Line2D
from skimage.measure import find_contours


from utils import remove_nans, StandardModel, DarkSun

OMEGA_CDM_H2 = 0.1198
SI_BOUND = 457.281  # 0.1 cm^2 / g

sm = StandardModel()


def find_cont(xs, ys, zs, level=SI_BOUND):
    contour = find_contours(zs, level)[0]

    jjs, iis = contour.T

    new_ys = np.interp(jjs, np.arange(len(ys)), ys)
    new_xs = np.interp(iis, np.arange(len(xs)), xs)

    return new_xs, new_ys


if __name__ == "__main__":

    # Load the benchmark data and extract needed quantities
    data = pd.read_csv("../../rundata/bm_lec1=0.1_lec2=1.0.csv")
    # data = pd.read_csv("../../rundata/bm_lec1=1.0_lec2=0.0.csv")
    # data = pd.read_csv("../../rundata/bm_lec1=0.1_lec2=1.0_xi_inf=5e-2.csv")
    # data = pd.read_csv("../../rundata/bm_lec1=0.1_lec2=1.0_xi_inf=1e-1.csv")

    data.sort_values(by=["LAM", "N"], inplace=True)

    c = data.iloc[0]["C"]
    lec1 = data.iloc[0]["LEC1"]
    lec2 = data.iloc[0]["LEC2"]
    mu_eta = data.iloc[0]["MU_ETA"]
    mu_del = data.iloc[0]["MU_DEL"]
    xi_inf = data.iloc[0]["XI_INF"]

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

    # Dneff
    xifo_array = np.array(data["XI_FO"]).reshape((len(lams), len(ns)))
    tsmfo_array = np.array(data["TSM_FO"]).reshape((len(lams), len(ns)))
    remove_nans(xifo_array)
    remove_nans(tsmfo_array)

    eta_rd_cont = find_cont(ns, lams, rde_array, level=OMEGA_CDM_H2)

    xifo_interp = RectBivariateSpline(ns, np.log10(lams), xifo_array)
    tsmfo_interp = RectBivariateSpline(ns, np.log10(lams), xifo_array)

    models = [
        DarkSun(n, lam, c, lec1, lec2, mu_eta, mu_del, xi_inf)
        for n, lam in zip(eta_rd_cont[0], eta_rd_cont[1])
    ]

    # Evaluate on contour
    xifo = [
        xifo_interp(x, y)[0][0]
        for x, y in zip(eta_rd_cont[0], np.log10(eta_rd_cont[1]))
    ]
    tsmfo = [
        tsmfo_interp(x, y)[0][0]
        for x, y in zip(eta_rd_cont[0], np.log10(eta_rd_cont[1]))
    ]
    energy_ratio = [
        mod.dark_geff(xi * tsm) * xi ** 4 / sm.geff(tsm)
        for xi, tsm, mod in zip(xifo, tsmfo, models)
    ]
    entropy_ratio = [
        mod.dark_heff(xi * tsm) * xi ** 3 / sm.heff(tsm)
        for xi, tsm, mod in zip(xifo, tsmfo, models)
    ]

    plt.plot(eta_rd_cont[0], energy_ratio, label="gd/gsm")
    plt.plot(eta_rd_cont[0], entropy_ratio, label="sd/ssm")
    plt.ylabel(r"$s_d/s_{\mathrm{sm}}$", fontsize=16)
    plt.xlabel(r"$N$", fontsize=16)
    plt.yscale("log")
    plt.xlim([5, 30])
    plt.legend()
    # plt.ylim([1e-4, 1e-1])
    plt.savefig(
        "/home/logan/Research/DarkSun/cpp/analysis/figures/energy_ratio_xi_inf=1e-2.pdf"
    )

