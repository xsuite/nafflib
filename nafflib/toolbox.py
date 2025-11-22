import numpy as np

from .optimize import point_dfft,expArray
from .windowing import hann
from .naff import fundamental_frequency, N_ARANGE



# ---------------------------------------
# Taken from 10.5170/CERN-1994-002, eq. 5.1.2
def henon_map(x, px, Q, num_turns):
    """
    Simulates the Henon map for 2D phase space.

    Parameters
    ----------
    x : float
        Initial position of the particle.
    px : float
        Initial momentum of the particle.
    Q : float
        Tune of the map.
    num_turns : int
        Number of turns to simulate.

    Returns
    -------
    x,px : tuple of ndarray
        Arrays representing the position and momentum of the particle at each turn.
    """
    z_vec = np.nan * np.ones(num_turns) + 1j * np.nan * np.ones(num_turns)
    z_vec[0] = x - 1j * px
    for ii in range(num_turns - 1):
        _z = z_vec[ii]

        z_vec[ii + 1] = np.exp(2 * np.pi * 1j * Q) * (
            _z - 1j / 4 * (_z + np.conjugate(_z)) ** 2
        )
    return np.real(z_vec), -np.imag(z_vec)


# ---------------------------------------


# ---------------------------------------
# Taken from 10.5170/CERN-1994-002, eq. 7.1.5
def henon_map_4D(x, px, y, py, Qx, Qy, coupling, num_turns):
    """
    Simulates the 4D Henon map for phase space, incorporating coupling between two oscillations.

    Parameters
    ----------
    x, y : float
        Initial positions in two orthogonal dimensions.
    px, py : float
        Initial momenta in two orthogonal dimensions.
    Qx, Qy : float
        Tunes of the map in each dimension.
    coupling : float
        Coupling strength between the dimensions.
    num_turns : int
        Number of turns to simulate.

    Returns
    -------
    x,px,y,py : tuple of ndarray
        Arrays representing the position and momentum of the particle in each dimension at each turn.
    """

    z1_vec = np.nan * np.ones(num_turns) + 1j * np.nan * np.ones(num_turns)
    z2_vec = np.nan * np.ones(num_turns) + 1j * np.nan * np.ones(num_turns)
    z1_vec[0] = x - 1j * px
    z2_vec[0] = y - 1j * py
    for ii in range(num_turns - 1):
        _z1 = z1_vec[ii]
        _z2 = z2_vec[ii]

        _z1_sum = _z1 + np.conjugate(_z1)
        _z2_sum = _z2 + np.conjugate(_z2)
        z1_vec[ii + 1] = np.exp(2 * np.pi * 1j * Qx) * (
            _z1 - 1j / 4 * (_z1_sum**2 - coupling * _z2_sum**2)
        )
        z2_vec[ii + 1] = np.exp(2 * np.pi * 1j * Qy) * (
            _z2 + 1j / 2 * coupling * _z1_sum * _z2_sum
        )

    # Returns x,px,y,py
    return np.real(z1_vec), -np.imag(z1_vec), np.real(z2_vec), -np.imag(z2_vec)


# ---------------------------------------




# ---------------------------------------
def dfft(nu,z):
    return np.array([point_dfft(_f, z) for _f in nu])
#----------------------------------------


def naff_dfft(nu, z, N=None, num_harmonics=1, window_order=2, window_type="hann"):

    assert num_harmonics >= 1, "number_of_harmonics needs to be >= 1"

    # initialisation, creating a copy of the signal since we'll modify it
    # ---------------------
    N  = N_ARANGE(len(z))
    _z = z.copy()
    # ---------------------

    # Computing the window function
    # ---------------------
    window_fun = {"hann": hann}[window_type.lower()]
    w_of_N = window_fun(N, order=window_order)
    # ---------------------


    A_dfft = []
    for _ in range(num_harmonics):

        # Computing frequency and amplitude
        amp, freq = fundamental_frequency(_z*w_of_N)

        # Computing cfft
        _dfft = dfft(nu,_z*w_of_N)

        # Saving results
        A_dfft.append(_dfft)

        # Substraction procedure
        zgs = amp * expArray(2 * np.pi * 1j * freq, len(N))
        _z -= zgs

    return A_dfft
# ---------------------------------------



# ---------------------------------------
def generate_signal(amplitudes, frequencies, N):
    """
    Generates a complex signal from given amplitudes and frequencies over a range of turns.

    Parameters
    ----------
    amplitudes : list or float
        Amplitudes of the components of the signal.
    frequencies : list or float
        Frequencies of the components of the signal.
    N : ndarray
        Array of turn numbers for which the signal is generated.

    Returns
    -------
    x,px : tuple of ndarray
        Arrays representing the position and momentum of the particle at each turn.
    """

    if isinstance(amplitudes, (float, int)):
        amplitudes = [amplitudes]
    if isinstance(frequencies, (float, int)):
        frequencies = [frequencies]

    assert len(amplitudes) == len(
        frequencies
    ), "Amplitudes and frequencies must have the same length"

    signal = sum(
        [
            A * np.exp(1j * (2 * np.pi * (Q) * N))
            for A, Q in zip(amplitudes, frequencies)
        ]
    )
    x = signal.real
    px = -signal.imag

    return x, px


# ---------------------------------------


# ---------------------------------------
def generate_pure_KAM(
    amplitudes, combinations, fundamental_tunes, N, return_frequencies=False
):
    """
    Generates a signal using the Kolmogorov-Arnold-Moser (KAM) theorem, simulating resonances.

    Parameters
    ----------
    amplitudes : list or float
        Amplitudes of the combinations.
    combinations : list of tuples or tuple
        Combination indices for linear combination of the fundamental tunes.
    fundamental_tunes : list
        Fundamental tunes for linear combinations.
    N : ndarray
        Array of turn numbers for signal generation.
    return_frequencies : bool, optional
        If True, also returns the frequencies used for signal generation.

    Returns
    -------
    x,px : tuple
        Arrays representing the position and momentum of the particle at each turn and, optionally, frequencies used for generation.
    """

    if isinstance(amplitudes, (float, int)):
        amplitudes = [amplitudes]
    if isinstance(combinations, (float, int)):
        combinations = [combinations]

    assert len(amplitudes) == len(
        combinations
    ), "amplitudes and resonances must have the same length"

    # Computing the frequencies
    Q_vec = fundamental_tunes + [1]
    n_vec = combinations
    assert (
        len(Q_vec) == np.shape(n_vec)[1]
    ), "combinations should have n+1 indices if n fundamental tunes are provided"
    frequencies = [np.dot(_r, Q_vec) for _r in n_vec]

    # Generating the signal
    signal = sum(
        [
            A * np.exp(1j * (2 * np.pi * (Q) * N))
            for A, Q in zip(amplitudes, frequencies)
        ]
    )
    x = signal.real
    px = -signal.imag

    if return_frequencies:
        return x, px, frequencies
    else:
        return x, px


# ---------------------------------------


