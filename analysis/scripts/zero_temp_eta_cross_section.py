"""
Script for generating plot of 2eta -> 4eta zero-temperature cross-section.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.offsetbox as offsetbox
from scipy.interpolate import UnivariateSpline


class EtaCrossSection:
    def __init__(self):
        self._load_cs_data()

    def _load_cs_data(self):
        # Load the data for contributions to the zero temp CS
        logzs = np.genfromtxt("../../rundata/cs_data/log10_zs.dat")
        logcs44 = np.genfromtxt("../../rundata/cs_data/log10_cs44.dat")
        logcs66 = np.genfromtxt("../../rundata/cs_data/log10_cs66.dat")
        logcs46 = np.genfromtxt("../../rundata/cs_data/log10_cs46.dat")

        # Create interpolating function
        self._scaled_cs_eta_44_interp = UnivariateSpline(
            logzs, logcs44, k=1, s=0
        )
        self._scaled_cs_eta_66_interp = UnivariateSpline(
            logzs, logcs66, k=1, s=0
        )
        self._scaled_cs_eta_46_interp = UnivariateSpline(
            logzs, logcs46, k=1, s=0
        )

    def scaled_cs_eta_44(self, scaled_cme):
        """
        Compute the scaled cross-section for 2eta->4eta using only 4pt
        interactions.

        Parameters
        ----------
        scaled_cme: float
            Center-of-mass energy scaled by the eta' mass: z = cme / meta

        Notes
        -----
        By 'scaled', we mean that all model parameters have been removed.
        Specifically a factor of [256 pi^4 lec1^2 meta^7 / 9 lam^8 n^2]^2. The
        cross-section is computed by interpolating between data generated using
        the RAMBO phase-space generator. This function includes only the 4pt
        eta interactions, i.e. the (d_mu eta)^4 term in the chiral Lagrangian.
        For energies beyond the interpolation range, we use a fit of the form
        log10(cs) = m * log10(z) + b with m = 14 and b='eta_cs_intercept44'.
        """
        logz = np.log10(scaled_cme)
        logzmin = 0.60206107706280998056
        logzmax = 2.0
        intercept = -11.116318726988425
        res = 0.0
        if logzmin <= logz <= logzmax:
            res = 10.0 ** self._scaled_cs_eta_44_interp(logz)
        elif logzmax < logz:
            res = scaled_cme ** 14 * 10.0 ** intercept
        return res

    def scaled_cs_eta_66(self, scaled_cme):
        """
        Compute the scaled cross-section for 2eta->4eta using only 6pt
        interactions

        Parameters
        ----------
        scaled_cme: float
            Center-of-mass energy scaled by the eta' mass: z = cme / meta

        Notes
        -----
        By 'scaled', we mean that all model parameters have been removed.
        Specifically a factor of [256 pi^4 lec2 meta^7 / 15 lam^8 n^2]^2. The
        cross-section is computed by interpolating between data generated using
        the RAMBO phase-space generator. This function includes only the 6pt eta
        interactions, i.e. the (d_mu eta)^6 term in the chiral Lagrangian. For
        energies beyond the interpolation range, we use a fit of the form
        log10(cs) = m * log10(z) + b with m = 14 and b='eta_cs_intercept66'.
        """
        logz = np.log10(scaled_cme)
        logzmin = 0.60206107706280998056
        logzmax = 2.0
        intercept = -12.038358477167012
        res = 0.0
        if logzmin <= logz <= logzmax:
            res = 10.0 ** self._scaled_cs_eta_66_interp(logz)
        elif logzmax < logz:
            res = scaled_cme ** 14 * 10.0 ** intercept
        return res

    def scaled_cs_eta_46(self, scaled_cme):
        """
        Compute the scaled cross-section for 2eta->4eta using only the
        interference between the 4pt and 6pt amplitudes.

        Parameters
        ----------
        scaled_cme: float
            Center-of-mass energy scaled by the eta' mass: z = cme / meta

        Notes
        -----
        By 'scaled', we mean that all model parameters have been removed.
        Specifically a factor of
            [256 pi^4 meta^7 / lam^8 n^2]^2 * [lec1^2 / 9] * [lec2 / 15]
        The cross-section is computed by interpolating between data generated
        using the RAMBO phase-space generator. This function includes only the
        4pt eta amplitude times the 6pt eta amplitude, i.e. the interference
        term: A4 * A6. For energies beyond the interpolation range, we use a
        fit of the form log10(cs) = m * log10(z) + b with m = 14 and
        b='eta_cs_intercept46'.
        """
        logz = np.log10(scaled_cme)
        logzmin = 0.60206107706280998056
        logzmax = 2.0
        intercept = -11.57800399152332
        res = 0.0
        if logzmin <= logz <= logzmax:
            res = 10.0 ** self._scaled_cs_eta_46_interp(logz)
        elif logzmax < logz:
            res = scaled_cme ** 14 * 10.0 ** intercept
        return res


class DarkSun(EtaCrossSection):
    def __init__(self, n=10, lam=0.1, lec1=0.1, lec2=1.0, mu_eta=1.0):
        super().__init__()
        self.n = n
        self.lam = lam
        self.lec1 = lec1
        self.lec2 = lec2
        self.mu_eta = mu_eta

    def meta(self):
        return self.mu_eta * self.lam / np.sqrt(self.n)

    def cross_section_2eta_4eta(self, cme):
        """
        Compute the full cross section for 2eta -> 4eta.

        Parameters
        ----------
        cme: float
            Center-of-mass energy in GeV.
        """
        scaled_cme = cme / self.meta()
        if scaled_cme < 4.0:
            return 0.0

        com = (
            6.9093374296577904e7
            * self.mu_eta ** 14
            / (self.lam ** 2 * self.n ** 11)
        )

        c44 = com * self.lec1 ** 4 / 9.0
        c66 = com * self.lec2 ** 2 / 25.0
        c46 = -2.0 * com * self.lec1 ** 2 * self.lec2 / 15.0

        return (
            c44 * self.scaled_cs_eta_44(scaled_cme)
            + c66 * self.scaled_cs_eta_66(scaled_cme)
            + c46 * self.scaled_cs_eta_46(scaled_cme)
        )


if __name__ == "__main__":
    custom_preamble = {
        "text.usetex": True,
        "text.latex.preamble": [
            r"\usepackage{amsmath}",  # for the align enivironment
        ],
    }
    plt.rcParams.update(custom_preamble)

    model = DarkSun(lec1=0.1)

    cmes = np.logspace(np.log10((4 + 1e-5)), 2.0, 150) * model.meta()

    model = DarkSun(lec1=1.0, lec2=0.0)
    css44 = np.array([model.cross_section_2eta_4eta(cme) for cme in cmes])
    model = DarkSun(lec1=0.0, lec2=1.0)
    css66 = np.array([model.cross_section_2eta_4eta(cme) for cme in cmes])
    model = DarkSun(lec1=1.0, lec2=1.0)
    cssmix = np.array([model.cross_section_2eta_4eta(cme) for cme in cmes])

    fig = plt.figure(dpi=200)
    plt.plot(cmes, css44, label=r"$\lambda_1 = 1, \lambda_2=0$", ls="-")
    plt.plot(cmes, css66, label=r"$\lambda_1 = 0, \lambda_2=1$", ls="--")
    plt.plot(cmes, cssmix, label=r"$\lambda_1 = 1, \lambda_2=1$", ls="-.")

    plt.yscale("log")
    plt.xscale("log")

    plt.ylabel(
        r"$\sigma(2\bar{\eta}'\to4\bar{\eta}') \ (\mathrm{GeV}^{-2})$",
        fontsize=16,
    )
    plt.xlabel(r"$\sqrt{s} \ (\mathrm{GeV})$", fontsize=16)

    plt.ylim([1e-11, 1e15])
    plt.xlim([np.min(cmes), np.max(cmes)])

    text1 = r"\begin{align*}"
    text2 = r"N &= 10\\"
    text3 = r"\Lambda &= 0.1\mathrm{GeV}"
    text4 = r"\end{align*}"
    text = text1 + text2 + text3 + text4

    plt.text(1.19, 1e-0, text, fontsize=16)
    plt.legend(fontsize=16, frameon=False, loc=4)

    plt.savefig("../figures/eta_cs.pdf")
