import numpy as np

from .windowing import hann
from .optimize import newton_method
from .toolbox import parse_real_signal


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


def fundamental_frequency(z, N=None, window_order=1, window_type="hann"):
    """
    Finds the fundamental frequency of a signal using the NAFF algorithm. It applies a window,
    estimates the frequency using FFT, and then refines it with the Newton method.

    Parameters
    ----------
    z : ndarray
        The complex signal array.
    N : ndarray or None, optional
        The array of turn numbers. If None, it defaults to a range equal to the length of z.
    window_order : int, optional
        The order of the windowing function.
    window_type : str, optional
        The type of windowing function to use. Defaults to 'hann'.

    Returns
    -------
    amplitude,f0 : tuple
        The amplitude and the fundamental frequency of the signal.
    """

    # Initialisation
    # ---------------------
    if N is None:
        N = np.arange(len(z))
    # ---------------------

    # Windowing of the signal
    # ---------------------
    window_fun = {"hann": hann}[window_type.lower()]
    z_w = z * window_fun(N, order=window_order)
    # ---------------------

    # Estimation of the main frequency with an FFT
    f0_est, resolution = _fft_f0_estimate(z_w)

    # Preparing the estimate for the Newton refinement method
    if f0_est >= 0.5:
        f0_est = -(1.0 - f0_est)
    f0_est = f0_est - resolution

    # Refinement of the tune calulation
    amplitude, f0 = newton_method(z_w, N, freq_estimate=f0_est, resolution=resolution)

    return amplitude, f0


def naff(z, num_harmonics=1, window_order=1, window_type="hann", to_pandas=False):
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
        The order of the windowing function.
    window_type : str, optional
        The type of windowing function to use. Defaults to 'hann'.
    to_pandas : bool, optional
        If True, returns a pandas DataFrame; otherwise, returns two arrays.

    Returns
    -------
    amplitudes,frequencies : tuple or DataFrame
        Harmonic amplitudes and frequencies, either as a DataFrame or as two separate arrays.
    """

    assert num_harmonics >= 1, "number_of_harmonics needs to be >= 1"

    # initialisation, creating a copy of the signal since we'll modify it
    # ---------------------
    N = np.arange(len(z))
    _z = z.copy()
    # ---------------------

    frequencies = []
    amplitudes = []
    for _ in range(num_harmonics):

        # Computing frequency and amplitude
        amp, freq = fundamental_frequency(
            _z, N=N, window_order=window_order, window_type=window_type
        )

        # Saving results
        frequencies.append(freq)
        amplitudes.append(amp)

        # Substraction procedure
        zgs = amp * np.exp(2 * np.pi * 1j * freq * N)
        _z -= zgs

    if to_pandas:
        import pandas as pd

        return pd.DataFrame({"amplitude": amplitudes, "frequency": frequencies})
    else:
        return np.array(amplitudes), np.array(frequencies)


def tune(x, px=None, window_order=1, window_type="hann"):
    """
    Computes the tune (fundamental frequency) of a signal or a complex signal formed from x and px.

    Parameters
    ----------
    x : ndarray
        The signal array or the real part of a complex signal.
    px : ndarray, optional
        The imaginary part of the complex signal, if any. If None, x is treated as the full signal.
    window_order : int, optional
        The order of the windowing function.
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
        x, px = np.asarray(x), np.asarray(px)
        z = x - 1j * px
    else:
        if np.any(np.imag(np.asarray(x)) != 0):
            # x is complex!
            z = x
        else:
            x, px = np.asarray(x), 0
            z = x - 1j * px
    N = np.arange(len(z))
    # ---------------------

    amp, freq = fundamental_frequency(
        z, N=N, window_order=window_order, window_type=window_type
    )

    return freq


def harmonics(
    x, px=None, num_harmonics=1, window_order=1, window_type="hann", to_pandas=False
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
        The order of the windowing function.
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
        return naff(
            z,
            num_harmonics=num_harmonics,
            window_order=window_order,
            window_type=window_type,
            to_pandas=to_pandas,
        )

    # FOR REAL SIGNAL:
    # ---------------------
    else:
        # Looking for twice as many frequencies then parsing the complex conjugates
        amplitudes, frequencies = naff(
            z,
            num_harmonics=2 * num_harmonics,
            window_order=window_order,
            window_type=window_type,
        )

        return parse_real_signal(amplitudes, frequencies, to_pandas=to_pandas)


def multiparticle_tunes(x, px=None, window_order=1, window_type="hann"):
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
        The order of the windowing function used in frequency analysis. Defaults to 1.
    window_type : str, optional
        The type of windowing function to apply. Defaults to 'hann'.

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

    freq_i = np.zeros(np.shape(x)[0])
    for ii in range(n_particles):
        freq_i[ii] = tune(
            x[ii], px[ii], window_order=window_order, window_type=window_type
        )

    return freq_i


def multiparticle_harmonics(x, px=None, num_harmonics=1, window_order=1, window_type="hann"):
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
        The order of the windowing function used in frequency analysis. Defaults to 1.
    window_type : str, optional
        The type of windowing function to apply. Defaults to 'hann'.

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
    # --------------------

    freq_i = []
    amp_i  = []
    for _x,_px in zip(x,px):
        _A,_Q = harmonics(
            _x,_px, num_harmonics=num_harmonics, window_order=window_order, window_type=window_type
        )
        amp_i.append(_A)
        freq_i.append(_Q)
        

    return np.array(amp_i),np.array(freq_i)