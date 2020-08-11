"""
Script for generating plot of 2eta -> 4eta zero-temperature cross-section.
"""

import numpy as np
import matplotlib.pyplot as plt


if __name__ == "__main__":
    logzs = np.genfromtxt("../../rundata/cs_data/log10_zs.dat")
    logcs44 = np.genfromtxt("../../rundata/cs_data/log10_cs44.dat")
    logcs66 = np.genfromtxt("../../rundata/cs_data/log10_cs66.dat")
    logcs46 = np.genfromtxt("../../rundata/cs_data/log10_cs46.dat")

    plt.figure(dpi=200)
    plt.plot(logzs, logcs44, label=r"$4\mathrm{pt-only}$", lw=1)
    plt.plot(logzs, logcs66, label=r"$6\mathrm{pt-only}$", lw=1, ls="--")
    plt.plot(logzs, logcs46, label=r"$\mathrm{interference}$", lw=1, ls="-.")

    plt.xlabel(r"$\mathrm{log}_{10}(z=\mathrm{CME}/m_{\eta'})$", fontsize=16)
    plt.ylabel(r"$\mathrm{log}_{10}(\bar{\sigma})$", fontsize=16)

    plt.ylim([-10, 10])
    plt.xlim([np.log10(4 + 1e-5), 1.5])

    plt.legend()
    plt.savefig("../figures/eta_cs.pdf")
