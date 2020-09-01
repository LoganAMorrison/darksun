"""
script for generating plot of delta relic density as a function of c and n
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from utils import remove_nans, OMEGA_CDM_H2
from scipy.interpolate import interp2d

if __name__ == "__main__":
    # Load the benchmark data and extract needed quantities
    data1 = pd.read_csv("../../rundata/c_vs_n_lam=10.0.csv")
    data2 = pd.read_csv("../../rundata/c_vs_n_lam=1e-1.csv")
    data3 = pd.read_csv("../../rundata/c_vs_n_lam=1e-3.csv")
    data4 = pd.read_csv("../../rundata/c_vs_n_lam=1e-5.csv")
    data5 = pd.read_csv("../../rundata/c_vs_n_lam=1e-7.csv")

    datas = [data1, data2, data3, data4, data5]

    ns = np.unique(np.array(data1["N"]))
    cs = np.unique(np.array(data1["C"]))

    n_array = np.array(data1["N"]).reshape((len(cs), len(ns)))
    c_array = np.array(data1["C"]).reshape((len(cs), len(ns)))

    plt.figure(dpi=100)

    rdd_array = np.array(data1["RD_DEL"]).reshape((len(cs), len(ns)))
    remove_nans(rdd_array)
    rdd = interp2d(ns, np.log10(cs), rdd_array)

    nss = np.linspace(np.min(ns), 12, 200)
    css = np.logspace(np.log10(np.min(cs)), np.log10(np.max(cs)), 200)

    plt.contour(
        ns, css, rdd(ns, np.log10(css)), levels=[OMEGA_CDM_H2], colors="k",
    )
    plt.contour(
        ns, css, rdd(ns, np.log10(css)), levels=[OMEGA_CDM_H2 / 3], colors="firebrick",
    )
    plt.contour(
        nss,
        css,
        rdd(nss, np.log10(css)),
        levels=[OMEGA_CDM_H2 * 3],
        colors="steelblue",
    )

    plt.contourf(
        ns,
        css,
        rdd(ns, np.log10(css)),
        levels=[OMEGA_CDM_H2 / 3.0, OMEGA_CDM_H2, 3.0 * OMEGA_CDM_H2],
        colors=["firebrick", "steelblue"],
        alpha=0.5,
    )

    plt.ylabel(r"$c$", fontdict={"size": 16})
    plt.xlabel(r"$N$", fontdict={"size": 16})

    plt.yscale("log")
    # plt.xscale("log")
    plt.xlim([np.min(ns), 12])
    plt.show()
