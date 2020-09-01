"""
Script for generating plot of 2eta -> 4eta thermally-averaged cross-section.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


if __name__ == "__main__":
    data = pd.read_csv("../../rundata/tc_data.csv")

    xs = np.array(data["X"])
    tc24s = np.array(data["TC24"])
    tc42s = np.array(data["TC42"])

    logxs = np.log(xs)

    plt.figure(dpi=200)
    plt.plot(xs, tc24s, label=r"$2\to4$", lw=1.5)
    plt.plot(xs, tc42s, label=r"$4\to2$", lw=1.5, ls="--")

    plt.yscale("log")
    plt.xscale("log")

    plt.xlabel(r"$x=m_{\bar{\eta}'}/T$", fontsize=16)

    plt.xlim([np.min(xs), np.max(xs)])
    plt.ylim([1e-30, 1e20])

    plt.text(
        10,
        1e-13,
        r"$\langle\sigma v\rangle_{2\bar{\eta}'\to4\bar{\eta}'}$",
        fontsize=16,
    )
    plt.text(
        10,
        1e9,
        r"$\langle\sigma v\rangle_{4\bar{\eta}'\to2\bar{\eta}'}$",
        fontsize=16,
    )

    plt.tight_layout()
    # plt.show()
    plt.savefig("../figures/eta_tcs.pdf")
