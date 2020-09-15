#!/usr/bin/python

"""
Script for generating the plots of xifo vs. N
"""
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scipy.ndimage import gaussian_filter
from scipy.interpolate import interp2d, interp1d
from scipy.interpolate import RectBivariateSpline
from matplotlib.lines import Line2D
from skimage.measure import find_contours


from utils import (
    remove_nans,
    StandardModel,
    DarkSun,
    OMEGA_CDM_H2,
    SI_BOUND,
    find_contour,
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

LABELS = [
    r"$\lambda_1=1,\lambda_2=0,\xi_{\infty}=10^{-2}$",
    r"$\lambda_1=0.1,\lambda_2=1,\xi_{\infty}=10^{-2}$",
    r"$\lambda_1=0.1,\lambda_2=1,\xi_{\infty}=5\times10^{-2}$",
    r"$\lambda_1=10^{-3},\lambda_2=1,\xi_{\infty}=10^{-2}$",
]

LINE_STYLES = [":", "--", "-.", "-"]

sm = StandardModel()


if __name__ == "__main__":

    plt.figure(dpi=150)

    for i, file in enumerate(DATA_FILES):
        # Load the benchmark data and extract needed quantities
        data = pd.read_csv(file)

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

        # Extract xi_fo and tsm_fo
        xifo_array = np.array(data["XI_FO"]).reshape((len(lams), len(ns)))
        tsmfo_array = np.array(data["TSM_FO"]).reshape((len(lams), len(ns)))
        remove_nans(xifo_array)
        remove_nans(tsmfo_array)

        eta_rd_cont = find_contour(ns, lams, rde_array, level=OMEGA_CDM_H2)

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
        dgeffs = [
            mod.dark_geff(xi * tsm)
            for xi, tsm, mod in zip(xifo, tsmfo, models)
        ]
        entropy_ratio = [
            mod.dark_geff(xi * tsm) * xi ** 4 / sm.geff(tsm)
            for xi, tsm, mod in zip(xifo, tsmfo, models)
        ]

        plt.plot(
            eta_rd_cont[0], entropy_ratio, label=LABELS[i], ls=LINE_STYLES[i]
        )

    plt.legend()
    plt.ylabel(r"$\rho_d/\rho_{\mathrm{sm}}$", fontsize=16)
    plt.xlabel(r"$N$", fontsize=16)
    plt.yscale("log")
    plt.xlim([5, 30])
    # plt.ylim([1e-4, 1e-1])
    plt.savefig(
        "/home/logan/Research/DarkSun/cpp/analysis/figures/energy_ratios.pdf"
    )

