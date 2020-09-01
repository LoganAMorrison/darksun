"""
utils.py

Contains a set of functions for helping in analyzing the data from the dark
SU(N) simulations.
"""

import numpy as np


OMEGA_CDM_H2 = 0.1198
SI_BOUND = 457.281  # 0.1 cm^2 / g


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
