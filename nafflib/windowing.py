import numpy as np
import math


def hann(N, order=1):
    """
    Generates Hann's window, centered at the middle of the dataset.
    Hann's window is used to smooth data and reduce spectral leakage.

    Parameters
    ----------
    N : numpy.ndarray or similar
        Index array of the turn. Each element represents a specific turn for which
        the window value is calculated.
    Nt : int, optional
        Total number of turns. If not provided, it is determined as the maximum
        value in N. Defaults to None.
    p : int, optional
        The power to which the Hann function is raised. Defaults to 1.

    Returns
    -------
    window_signal : numpy.ndarray
        The Hann window values for each turn index in N. The window is centered
        around the middle of the dataset.

    Notes
    -----
    The formula used for the Hann window is:
    (2^p)*math.factorial(p)^2 / (math.factorial(2*p)) * (1 + cos(2*pi*(N-center)/Nt))^p
    """
    # Initialisation
    # ---------------
    Nt = np.max(N) + 1
    p = order
    # ---------------

    if np.mod(Nt, 2) == 0:
        center = Nt // 2 - 1
    else:
        center = Nt // 2
    return (
        (2**p)
        * math.factorial(p) ** 2
        / (math.factorial(2 * p))
        * (1 + np.cos(2 * np.pi * (N - center) / Nt)) ** (p)
    )
