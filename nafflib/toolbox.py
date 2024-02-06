import numpy as np


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
def parse_real_signal(amplitudes, frequencies, conjugate_tol=1e-10, to_pandas=False):
    """
    Parses a real signal from its complex representation, identifying and handling complex conjugates.

    Parameters
    ----------
    amplitudes : ndarray
        Amplitudes of the complex signal.
    frequencies : ndarray
        Frequencies of the complex signal.
    conjugate_tol : float, optional
        Tolerance for identifying complex conjugates. Defaults to 1e-10.
    to_pandas : bool, optional
        If True, returns a pandas DataFrame; otherwise, returns two arrays.

    Returns
    -------
    amplitudes,frequencies : DataFrame or tuple of ndarray
        Amplitudes and frequencies of the real signal components, either as a DataFrame or as two arrays.
    """

    A, Q = amplitudes, frequencies
    phasors = A * np.exp(2 * np.pi * 1j * Q)

    freq = []
    amp = []

    for i in range(len(Q) - 1):

        # Finding closest conjugate pair
        comparisons = phasors[i] + phasors
        pair_idx = np.argmin(np.abs(np.imag(comparisons)))
        pair_A = np.array([A[i], A[pair_idx]])
        pair_Q = np.array([Q[i], Q[pair_idx]])

        # Check if the pair is a complex conjugate, otherwise both freqs are kept
        if np.abs(np.diff(np.abs(pair_Q)))[0] > conjugate_tol:
            freq.append(Q[i])
            amp.append(A[i])
            continue

        # Creating avg amplitude and freq
        real = np.mean(np.real(pair_A))
        imag = np.mean(np.abs(np.imag(pair_A)))
        sign = np.sign(np.imag(pair_A))[pair_Q >= 0]

        if pair_Q[0] == pair_Q[1]:
            # the pair is a copy of itself (DC signal case)
            freq.append(pair_Q[0])
            amp.append(real + 1j * imag)
        else:
            # Complex conjugate found, adding only once
            if np.mean(np.abs(pair_Q)) not in freq:
                if np.any(np.isnan(pair_Q)):
                    freq.append(np.nan)
                    amp.append(np.nan + 1j * np.nan)
                else:
                    freq.append(np.mean(np.abs(pair_Q)))
                    amp.append(2 * (real + sign[0] * 1j * imag))

    if to_pandas:
        import pandas as pd

        return pd.DataFrame({"amplitude": amp, "frequency": freq})
    else:
        return np.array(amp), np.array(freq)


# ---------------------------------------


# ---------------------------------------
def find_linear_combinations(
    frequencies, fundamental_tunes=[], max_harmonic_order=10, r_vec=None, to_pandas=False
):
    """
    Identifies linear combinations of fundamental tunes that closely match the given frequencies.

    Parameters
    ----------
    frequencies : ndarray
        Array of frequencies to analyze.
    fundamental_tunes : list, optional
        List of fundamental tunes for resonance analysis. Length can be 1, 2, or 3.
    max_harmonic_order : int, optional
        Maximum order of harmonics to consider.
    r_vec : list of ndarrays, optional
        List of arrays representing the possible combinations of the fundamental tunes. If not provided, it is generated. (for big arrays, it is better to provide it to avoid memory issues when looping)
    to_pandas : bool, optional
        If True, returns a pandas DataFrame; otherwise, returns lists and an array.

    Returns
    -------
    r_values,error,f_values : DataFrame or tuple
        Resonance values, errors, and corresponding frequencies, either as a DataFrame or separate data structures.
    """

    assert len(fundamental_tunes) in [
        1,
        2,
        3,
    ], "need 1, 2 or 3 fundamental tunes (2D,4D,6D)"

    # Create a 3D array of all possible combinations of r_vec
    idx = max_harmonic_order
    if r_vec is None:
        if len(fundamental_tunes) == 1:
            r1, r2 = np.mgrid[-idx : idx + 1, -idx : idx + 1]
            r_vec = [r1, r2]
        elif len(fundamental_tunes) == 2:
            r1, r2, r3 = np.mgrid[-idx : idx + 1, -idx : idx + 1, -idx : idx + 1]
            r_vec = [r1, r2, r3]
        else:
            r1, r2, r3, r4 = np.mgrid[
                -idx : idx + 1, -idx : idx + 1, -idx : idx + 1, -idx : idx + 1
            ]
            r_vec = [r1, r2, r3, r4]
    else:
        assert len(r_vec) == len(fundamental_tunes) + 1, "r_vec must be of length fundamental_tunes+1"

    # Computing all linear combinations of r1*Qx + r2*Qy + r3*Qz + m
    Q_vec = fundamental_tunes + [1]
    all_combinations = sum([_r * _Q for _r, _Q in zip(r_vec, Q_vec)])

    # Find the closest combination for each frequency
    r_values = []
    f_values = []
    err = []
    for freq in frequencies:

        # Find the index of the closest combination
        closest_idx = np.unravel_index(
            np.argmin(np.abs(freq - all_combinations)), all_combinations.shape
        )

        # Get the corresponding values for r1,r2,r3,r4
        closest_combination = tuple(_r[closest_idx] for _r in r_vec)
        closest_value = all_combinations[closest_idx]

        r_values.append(closest_combination)
        f_values.append(closest_value)
        err.append(np.abs(closest_value - freq))

    if to_pandas:
        import pandas as pd

        return pd.DataFrame({"resonance": r_values, "err": err, "frequency": f_values})
    else:
        return [tuple(_r) for _r in r_values], np.array(err), np.array(f_values)


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
    r_vec = combinations
    assert (
        len(Q_vec) == np.shape(r_vec)[1]
    ), "combinations should have n+1 indices if n fundamental tunes are provided"
    frequencies = [np.dot(_r, Q_vec) for _r in r_vec]

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
