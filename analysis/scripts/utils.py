"""
utils.py

Contains a set of functions for helping in analyzing the data from the dark
SU(N) simulations.
"""

import numpy as np
from scipy.special import kv
from scipy.interpolate import UnivariateSpline
from sm_data import sm_geff_data, sm_heff_data, sm_sqrt_gstar_data, sm_temps
from skimage.measure import find_contours as _find_contours


OMEGA_CDM_H2 = 0.1198
SI_BOUND = 457.281  # 0.1 cm^2 / g
NEFF_CMB_BOUND = (2.92, 0.36)
NEFF_BBN_BOUND = (2.85, 0.28)


def find_contour(xs, ys, zs, level):
    """
    Find the controur.
    """
    contour = _find_contours(zs, level)[0]

    jjs, iis = contour.T

    new_ys = np.interp(jjs, np.arange(len(ys)), ys)
    new_xs = np.interp(iis, np.arange(len(xs)), xs)

    return new_xs, new_ys


def remove_nans(arr):
    """
    Removes the nans from an array by replacing them with the average of the
    surrounding values.

    Parameters
    ----------
    arr: np.array
        Input array.

    """
    N, M = arr.shape
    done = False

    while not done:
        for n in range(N):
            for m in range(M):
                if np.isnan(arr[n, m]):
                    avg = 0.0
                    cnt = 0
                    if n - 1 >= 0:
                        if not np.isnan(arr[n - 1, m]):
                            avg += arr[n - 1, m]
                            cnt += 1
                    if n + 1 <= N - 1:
                        if not np.isnan(arr[n + 1, m]):
                            avg += arr[n + 1, m]
                            cnt += 1
                    if m - 1 >= 0:
                        if not np.isnan(arr[n, m - 1]):
                            avg += arr[n, m - 1]
                            cnt += 1
                    if m + 1 <= M - 1:
                        if not np.isnan(arr[n, m + 1]):
                            avg += arr[n, m + 1]
                            cnt += 1
                    if cnt > 0:
                        arr[n, m] = avg / cnt
        done = np.sum(np.isnan(arr)) == 0


class DarkSun:
    def __init__(self, n, lam, c, lec1, lec2, mu_eta, mu_del, xi_inf):
        self.n = n
        self.lam = lam
        self.c = c
        self.lec1 = lec1
        self.lec2 = lec2
        self.mu_eta = mu_eta
        self.mu_del = mu_del
        self.xi_inf = xi_inf

    def m_eta(self):
        """
        Compute the mass of the eta.
        """
        return self.lam * self.mu_eta / np.sqrt(self.n)

    def m_del(self):
        """
        Compute the mass of the delta.
        """
        return self.lam * self.mu_del * self.n

    def g_del(self):
        """
        Compute the d.o.f. or the delta.
        """
        return self.n + 1

    def dark_heff(self, td):
        """
        Compute the effective d.o.f. in entropy of the dark sector.

        Parameters
        ----------
        td: float
            The dark sector temperature.
        """
        xe = self.m_eta() / td
        xd = self.m_del() / td

        ge = 1.0
        gd = self.g_del()

        pre = 45.0 / (4.0 * np.pi ** 4)
        pree = pre * ge * xe ** 3
        pred = pre * gd * xd ** 3

        bess_sum_e = sum(
            1.0 / (1.0 + k) * kv(3, (1.0 + k) * xe) for k in range(5)
        )
        bess_sum_d = kv(3, xd)

        return pree * bess_sum_e + pred * bess_sum_d

    def dark_geff(self, td):
        """
        Compute the effective d.o.f. in energy of the dark sector.

        Parameters
        ----------
        td: float
            The dark sector temperature.
        """
        xe = self.m_eta() / td
        xd = self.m_del() / td

        ge = 1.0
        gd = self.g_del()

        pre = 30.0 / (2.0 * np.pi ** 4)
        pree = pre * ge * xe ** 2
        pred = pre * gd * xd ** 2

        bess_sum_e = sum(
            1.0
            / (1.0 + k) ** 2
            * (
                (1.0 + k) * xe * kv(1, (1.0 + k) * xe)
                + 3.0 * kv(2, (1.0 + k) * xe)
            )
            for k in range(5)
        )
        bess_sum_d = xd * kv(1, xd) + 3.0 * kv(2, xd)

        return pree * bess_sum_e + pred * bess_sum_d

    def cross_section_2eta_2eta(self):
        """
        Compute the self interaction cross section of the eta.
        """
        return (
            62.012553360599640
            * self.mu_eta ** 6
            * self.lec1 ** 2
            / (self.lam ** 2 * self.n ** 5)
        )

    def cross_section_2del_2del(self):
        """
        Compute the self interaction cross section of the delta.
        """
        return 4.0 * np.pi ** 3 / self.lam ** 2


class StandardModel:
    _LOG_TEMP_MIN = -4.5
    _LOG_TEMP_MAX = 4.0
    _LOG_TEMP_STP = 0.025

    def __init__(self):
        self._heff_interp = UnivariateSpline(
            sm_temps, sm_heff_data, s=0, ext=3
        )
        self._geff_interp = UnivariateSpline(
            sm_temps, sm_heff_data, s=0, ext=3
        )
        self._sqrt_gstar_interp = UnivariateSpline(
            sm_temps, sm_sqrt_gstar_data, s=0, ext=3
        )

    def heff(self, tsm):
        """
        Returns the effective number of d.o.f. in entropy of the SM.
        """
        ltsm = np.log10(tsm)
        return self._heff_interp(ltsm)

    def geff(self, tsm):
        """
        Returns the effective number of d.o.f. in energy of the SM.
        """
        ltsm = np.log10(tsm)
        return self._geff_interp(ltsm)

    def sqrt_gstar(self, tsm):
        """
        Returns the square-root of gstar of the SM.
        """
        ltsm = np.log10(tsm)
        return self._sqrt_gstar_interp(ltsm)

