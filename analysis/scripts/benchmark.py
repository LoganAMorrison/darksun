#!/usr/bin/python

"""
Script for generating the benchmark plots for the Dark SU(N) model.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OMEGA_CDM_H2 = 0.1198
SI_BOUND = 457.281  # 0.1 cm^2 / g

if __name__ == "__main__":

    # Load the benchmark data and extract needed quantities
    data = pd.read_csv("../../rundata/bm3.csv")

    ns = np.unique(np.array(data["N"]))
    lams = np.unique(np.array(data["LAM"]))

    n_array = np.array(data["N"]).reshape((len(lams), len(ns)))
    lam_array = np.array(data["LAM"]).reshape((len(lams), len(ns)))
    # Relic densities
    rde_array = np.array(data["RD_ETA"]).reshape((len(lams), len(ns)))
    rdd_array = np.array(data["RD_DEL"]).reshape((len(lams), len(ns)))
    rdt_array = rde_array + rdd_array  # total

    sie_array = np.array(data["ETA_SI_PER_MASS"]).reshape((len(lams), len(ns)))
    sid_array = np.array(data["DEL_SI_PER_MASS"]).reshape((len(lams), len(ns)))
    sia_array = (
        sid_array * rdd_array ** 2 + sie_array * rde_array ** 2
    ) / OMEGA_CDM_H2 ** 2

    print(n_array.shape)

    plt.figure(dpi=100)
    plt.contour(
        n_array, lam_array, rde_array, levels=[0.12], colors=["mediumorchid"]
    )
    plt.contour(
        n_array, lam_array, rdd_array, levels=[0.12], colors=["firebrick"]
    )
    plt.contourf(
        n_array, lam_array, sie_array, levels=[SI_BOUND, 1e30], alpha=0.5
    )
    plt.contourf(
        n_array, lam_array, sid_array, levels=[SI_BOUND, 1e30], alpha=0.5
    )

    plt.yscale("log")
    plt.xscale("log")
    plt.show()

