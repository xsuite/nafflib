import numpy as np

from .windowing import hann
from .optimize import newton_method
from .toolbox import parse_real_signal



def _fft_f0_estimate(z,force_len = None):
    """
    Estimate the main frequency using an FFT. The signal is cropped to the closest power 
    of 2 for improved accuracy.

    Parameters
    ----------
    z : numpy.ndarray
        Complex array of the signal.
    n_forced : int, optional
        Number of turns to use for the FFT. If > len(z), the signal is padded 
        with zeros. Defaults to None.

    Returns
    -------
    tuple of float
        A tuple containing the estimated tune and the resolution.
    """
    # Cropping signal to closest power of 2
    if force_len is None:
        force_len = 2**int(np.log2(len(z)))

    # Search for maximum in Fourier spectrum
    z_spectrum = np.fft.fft(z,n=force_len)
    idx_max    = np.argmax(np.abs(z_spectrum))
    
    # Estimation of Tune with FFT
    tune_est   = idx_max/force_len
    resolution = 1/force_len

    return tune_est,resolution

def fundamental_frequency(z,N = None,window_order = 1,window_type = 'hann'):
    """
    Subroutine of the NAFF algorithm. 
    1. Applies a Hann window
    2. Estimates the fundamental frequency with an FFT 
    3. Refines the estimate with a complex Newton method
    4. Returns the main frequency and the amplitude.
    """

    # Initialisation
    #---------------------
    if N is None:
        N = np.arange(len(z))
    #---------------------

    # Windowing of the signal
    #---------------------
    window_fun = {'hann':hann}[window_type.lower()]
    z_w = z * window_fun(N,order=window_order)
    #---------------------

    # Estimation of the main frequency with an FFT
    f0_est,resolution = _fft_f0_estimate(z_w)

    # Preparing the estimate for the Newton refinement method
    if f0_est >= 0.5:
        f0_est = -(1.0 - f0_est)
    f0_est = f0_est - resolution

    # Refinement of the tune calulation
    amplitude,f0 = newton_method(z_w,N,freq_estimate = f0_est,resolution = resolution)

    return amplitude,f0


def naff(z,num_harmonics = 1,window_order = 1,window_type = 'hann',to_pandas = False):
    """
    Applies the NAFF algorithm to find the spectral lines of a signal.
    """

    assert num_harmonics >=1, 'number_of_harmonics needs to be >= 1'
    
    # initialisation, creating a copy of the signal since we'll modify it
    #---------------------
    N  = np.arange(len(z))
    _z = z.copy()
    #---------------------


    frequencies = []
    amplitudes  = [] 
    for _ in range(num_harmonics):

        # Computing frequency and amplitude
        amp,freq  = fundamental_frequency(_z,N=N,window_order= window_order,
                                                window_type = window_type)

        # Saving results
        frequencies.append(freq)
        amplitudes.append(amp)

        # Substraction procedure
        zgs  = amp * np.exp(2 * np.pi * 1j * freq * N)
        _z   -= zgs

    if to_pandas:
        import pandas as pd
        return pd.DataFrame({'amplitude':amplitudes,'frequency':frequencies})
    else:
        return np.array(amplitudes),np.array(frequencies)


def tune(x,px=None,window_order = 1,window_type = 'hann'):

    # initialisation
    #---------------------
    if px is not None:
        x,px = np.asarray(x),np.asarray(px)
        z  = x - 1j*px
    else:
        if np.any(np.imag(np.asarray(x)) != 0):
            # x is complex! 
            z = x
        else:
            x,px = np.asarray(x),0
            z  = x - 1j*px
    N  = np.arange(len(z))
    #---------------------

    amp,freq  = fundamental_frequency(z,N=N,window_order= window_order,
                                            window_type = window_type)
    
    return np.abs(freq)


def harmonics(x,px = None,num_harmonics = 1,window_order = 1,window_type = 'hann',to_pandas = False):
    # initialisation
    #---------------------
    if px is not None:
        real_signal = False
        x,px = np.asarray(x),np.asarray(px)
        z  = x - 1j*px
    else:
        if np.any(np.imag(np.asarray(x)) != 0):
            real_signal = False
            z = x
        else:
            real_signal = True
            x,px = np.asarray(x),0
            z  = x - 1j*px
    #---------------------
            

    # FOR COMPLEX SIGNAL:
    #---------------------
    if not real_signal:
        return naff(z,  num_harmonics = num_harmonics,
                        window_order = window_order,
                        window_type = window_type,
                        to_pandas = to_pandas)
    
    # FOR REAL SIGNAL:
    #---------------------
    else:
        # Looking for twice as many frequencies then parsing the complex conjugates
        amplitudes,frequencies = naff(z, num_harmonics = 2*num_harmonics,
                                        window_order = window_order,
                                        window_type = window_type)
        
        return parse_real_signal(amplitudes,frequencies,to_pandas=to_pandas)




def multiparticle_tunes(x,px=None,window_order = 1,window_type = 'hann'):

    n_particles = np.shape(x)[0]

    # Initializating px
    #--------------------
    if px is None:
        px = n_particles*[None]
    #--------------------

    freq_i  = np.zeros(np.shape(x)[0])
    for ii in range(n_particles):
        freq_i[ii] = tune(x[ii],px[ii],window_order = window_order,window_type = window_type)
    
    return freq_i



