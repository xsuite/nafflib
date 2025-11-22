import numpy as np
import functools

from .windowing import hann
from .optimize import newton_method,expArray


@functools.lru_cache(maxsize=10) 
def N_ARANGE(size):
    assert size < 2e9, "size should be smaller than 2e9 (int32 limit)"
    return np.arange(size,dtype=np.int32)




def _parse_real_signal(amplitudes, frequencies, rel_conjugate_tol=1e-2):
    """
    Parses a real signal from its complex representation, identifying and handling complex conjugates.

    Parameters
    ----------
    amplitudes : ndarray
        Amplitudes of the complex signal.
    frequencies : ndarray
        Frequencies of the complex signal.
    rel_conjugate_tol : float, optional
        Relative tolerance for identifying complex conjugates: remainder/norm. Default is 1e-2.

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
        # Finding closest conjugate pair (comparison should be zero: z + z* = 2*Re(z))
        comparisons = np.abs((phasors[i] + phasors) - (2 * np.real(phasors[i])))
        # The frequency should at least be of opposite sign!
        comparisons[Q * Q[i] > 0] = np.inf
        # Make comparison relative with respect to the phasor under inspection
        comparisons = comparisons / np.abs(phasors[i])
        # --
        pair_idx = np.argmin(comparisons)
        pair_A = np.array([A[i], A[pair_idx]])
        pair_Q = np.array([Q[i], Q[pair_idx]])
        # --

        # Check if the pair is a complex conjugate, otherwise both freqs
        # will make it out of the loop
        if comparisons[pair_idx] > rel_conjugate_tol:
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
            amp.append(pair_A[0])
        else:
            # Complex conjugate found, adding only once
            if np.mean(np.abs(pair_Q)) not in freq:
                if np.any(np.isnan(pair_Q)):
                    freq.append(np.nan)
                    amp.append(np.nan + 1j * np.nan)
                else:
                    freq.append(np.mean(np.abs(pair_Q)))
                    amp.append(2 * (real + sign[0] * 1j * imag))

    return np.array(amp), np.array(freq)


def _fft_f0_estimate(z, force_len=None):
    """
    Estimate the main frequency of a signal using FFT. The signal is cropped or padded to a specified length
    or the closest power of 2 for improved accuracy.

    Parameters
    ----------
    z : ndarray
        The complex signal array.
    force_len : int, optional
        The length to which the signal is cropped or padded for the FFT. If None,
        the length is set to the closest power of 2. Defaults to None.

    Returns
    -------
    tune_est,resolution : tuple
        The estimated main frequency (tune) and the resolution of the estimation.
    """
    # Cropping signal to closest power of 2
    if force_len is None:
        force_len = 2 ** int(np.log2(len(z)))

    # Search for maximum in Fourier spectrum
    z_spectrum = np.fft.fft(z, n=force_len)
    idx_max = np.argmax(np.abs(z_spectrum))

    # Estimation of Tune with FFT
    tune_est = idx_max / force_len
    resolution = 1 / force_len

    return tune_est, resolution



def fundamental_frequency(z):
    """
    Finds the fundamental frequency of a signal using the NAFF algorithm. It applies a window,
    estimates the frequency using FFT, and then refines it with the Newton method.

    Parameters
    ----------
    z : ndarray
        The complex signal array.
    window_order : int, optional
        The order of the windowing function. Defaults to 2.
    window_type : str, optional
        The type of windowing function to use. Defaults to 'hann'.

    Returns
    -------
    amplitude,f0 : tuple
        The amplitude and the fundamental frequency of the signal.
    """

    # Initialisation
    # ---------------------
    N = N_ARANGE(len(z))
    # ---------------------

    # Estimation of the main frequency with an FFT
    f0_est, resolution = _fft_f0_estimate(z)

    # Preparing the estimate for the refinement method
    if f0_est >= 0.5:
        f0_est = -(1.0 - f0_est)

    # Refinement of the tune calulation
    amplitude, f0 = newton_method(z, N, freq_estimate=f0_est, resolution=resolution)

    return amplitude, f0


def naff(z,num_harmonics=1, window_order=2, window_type="hann"):
    """
    Applies the NAFF algorithm to find spectral lines of a signal. It identifies multiple harmonics
    of the signal, removes them, and repeats the process.

    Parameters
    ----------
    z : ndarray
        The complex signal array.
    num_harmonics : int, optional
        The number of harmonics to identify in the signal.
    window_order : int, optional
        The order of the windowing function. Defaults to 2.
    window_type : str, optional
        The type of windowing function to use. Defaults to 'hann'.

    Returns
    -------
    amplitudes,frequencies : tuple or DataFrame
        Harmonic amplitudes and frequencies, either as a DataFrame or as two separate arrays.
    """

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

    frequencies = []
    amplitudes = []
    for _ in range(num_harmonics):

        # Computing frequency and amplitude of the windowed signal
        amp, freq = fundamental_frequency(_z*w_of_N)

        # Saving results
        frequencies.append(freq)
        amplitudes.append(amp)

        # Substraction procedure
        zgs = amp * expArray(2 * np.pi * 1j * freq, len(N))
        _z -= zgs

    return np.array(amplitudes), np.array(frequencies)


def tune(x, px=None,window_order=2, window_type="hann"):
    """
    Computes the tune (fundamental frequency) of a signal or a complex signal formed from x and px.

    Parameters
    ----------
    x : ndarray
        The signal array or the real part of a complex signal.
    px : ndarray, optional
        The imaginary part of the complex signal, if any. If None, x is treated as the full signal.
    window_order : int, optional
        The order of the windowing function. Defaults to 2.
    window_type : str, optional
        The type of windowing function to use. Defaults to 'hann'.

    Returns
    -------
    freq : float
        The tune of the signal.
    """
    # initialisation
    # ---------------------
    if px is not None:
        real_signal = False
        x, px = np.asarray(x), np.asarray(px)
        z = x - 1j * px
    else:
        if np.any(np.imag(np.asarray(x)) != 0):
            real_signal = False
            z = x
        else:
            real_signal = True
            x, px = np.asarray(x), 0
            z = x - 1j * px

    # ---------------------

    # FOR COMPLEX SIGNAL:
    # ---------------------
    if not real_signal:
        amplitudes, frequencies = naff( z,
                                        num_harmonics=1,
                                        window_order=window_order,
                                        window_type=window_type,
                                    )
        freq = frequencies[0]

    # FOR REAL SIGNAL:
    # ---------------------
    else:
        # Looking for 2n+1 many frequencies (for DC) then parsing the complex conjugates
        amplitudes, frequencies = naff(
            z,
            num_harmonics=2 * 1 + 1,
            window_order=window_order,
            window_type=window_type,
        )
        amps, freqs = _parse_real_signal(amplitudes, frequencies)
        freq = freqs[0]
    return freq


def harmonics(
    x, px=None,num_harmonics=1, window_order=2, window_type="hann", to_pandas=False
):
    """
    Identifies harmonics in a signal or a complex signal formed from x and px using the NAFF algorithm.

    Parameters
    ----------
    x : ndarray
        The signal array or the real part of a complex signal.
    px : ndarray, optional
        The imaginary part of the complex signal, if any. If None, x is treated as the full signal.
    num_harmonics : int, optional
        The number of harmonics to identify in the signal.
    window_order : int, optional
        The order of the windowing function. Defaults to 2.
    window_type : str, optional
        The type of windowing function to use. Defaults to 'hann'.
    to_pandas : bool, optional
        If True, returns a pandas DataFrame; otherwise, returns two arrays.

    Returns
    -------
    amplitudes,frequencies : tuple or DataFrame
        Harmonic amplitudes and frequencies, either as a DataFrame or as two separate arrays.
    """
    # initialisation
    # ---------------------
    if px is not None:
        real_signal = False
        x, px = np.asarray(x), np.asarray(px)
        z = x - 1j * px
    else:
        if np.any(np.imag(np.asarray(x)) != 0):
            real_signal = False
            z = x
        else:
            real_signal = True
            x, px = np.asarray(x), 0
            z = x - 1j * px

    # ---------------------

    # FOR COMPLEX SIGNAL:
    # ---------------------
    if not real_signal:
        amplitudes, frequencies = naff(
            z,
            num_harmonics=num_harmonics,
            window_order=window_order,
            window_type=window_type,
        )

    # FOR REAL SIGNAL:
    # ---------------------
    else:
        # Looking for 2n+1 many frequencies (for DC) then parsing the complex conjugates
        amplitudes, frequencies = naff(
            z,
            num_harmonics=2 * num_harmonics + 1,
            window_order=window_order,
            window_type=window_type,
        )
        amplitudes, frequencies = _parse_real_signal(amplitudes, frequencies)
        amplitudes = amplitudes[:num_harmonics]
        frequencies = frequencies[:num_harmonics]

    # Final handling for pandas
    if to_pandas:
        import pandas as pd

        return pd.DataFrame({"amplitude": amplitudes, "frequency": frequencies})
    else:
        return amplitudes, frequencies



def _tune_wrapper(args):
    _x, _px, window_order, window_type = args
    return tune(
        _x, _px,window_order=window_order, window_type=window_type
    )

def multiparticle_tunes(x, px=None,window_order=2, window_type="hann", processes=None):
    """
    Calculates the tunes for multiple particles given their positions and optionally their momenta.
    This function computes the fundamental frequency for each particle individually.

    Parameters
    ----------
    x : ndarray or list of ndarray
        An array or a list of arrays containing the position data of each particle.
    px : ndarray or list of ndarray, optional
        An array or a list of arrays containing the momentum data of each particle.
        If None, only position data is used. Defaults to None.
    window_order : int, optional
        The order of the windowing function used in frequency analysis. Defaults to 2.
    window_type : str, optional
        The type of windowing function to apply. Defaults to 'hann'.
    processes : int, optional
        The number of processes to use in the calculation. If None, the multiprocessing is disabled.
        WARNING: if a small number of particles is used, the overhead of multiprocessing will slow down the calculation.

    Returns
    -------
    frequencies : ndarray
        An array containing the computed tunes for each particle.
    """
    n_particles = np.shape(x)[0]

    # Initializating px
    # --------------------
    if px is None:
        px = n_particles * [None]
    # --------------------
        
    if processes is not None:
        import multiprocessing as mp
        args = [(_x, _px, window_order, window_type) for _x, _px in zip(x, px)]
        
        with mp.Pool(processes=processes) as pool:
            results = pool.map(_tune_wrapper, args)

        freq_i = results

    else:
        freq_i = np.zeros(np.shape(x)[0])
        for ii in range(n_particles):
            freq_i[ii] = tune(
                x[ii], px[ii], window_order=window_order, window_type=window_type
            )

    return freq_i


def _harmonics_wrapper(args):
    _x, _px,num_harmonics, window_order, window_type = args
    return harmonics(
        _x, _px, num_harmonics=num_harmonics, window_order=window_order, window_type=window_type
    )

def multiparticle_harmonics(
    x, px=None,num_harmonics=1, window_order=2, window_type="hann", processes=None):
    """
    Calculates the harmonics for multiple particles given their positions and optionally their momenta.

    Parameters
    ----------
    x : ndarray or list of ndarray
        An array or a list of arrays containing the position data of each particle.
    px : ndarray or list of ndarray, optional
        An array or a list of arrays containing the momentum data of each particle.
        If None, only position data is used. Defaults to None.
    window_order : int, optional
        The order of the windowing function used in frequency analysis. Defaults to 2.
    window_type : str, optional
        The type of windowing function to apply. Defaults to 'hann'.
    processes : int, optional
        The number of processes to use in the calculation. If None, the multiprocessing is disabled.
        WARNING: if a small number of particles is used, the overhead of multiprocessing will slow down the calculation.

    Returns
    -------
    amplitudes,frequencies : tuple or DataFrame
        A list containing the computed harmonics for each particle.
    """
    n_particles = np.shape(x)[0]

    # Initializating px
    # --------------------
    if px is None:
        px = n_particles * [None]

    N = N_ARANGE(len(x[0]))
    # --------------------

    if processes is not None:
        import multiprocessing as mp
        args = [(_x, _px, num_harmonics, window_order, window_type) for _x, _px in zip(x, px)]
        
        with mp.Pool(processes=processes) as pool:
            results = pool.map(_harmonics_wrapper, args)

        amp_i, freq_i = zip(*results)
    
    else:
        freq_i = []
        amp_i = []
        for _x, _px in zip(x, px):
            _A, _Q = harmonics(
                _x,
                _px,
                num_harmonics=num_harmonics,
                window_order=window_order,
                window_type=window_type,
            )
            amp_i.append(_A)
            freq_i.append(_Q)

    return np.array(amp_i), np.array(freq_i)
