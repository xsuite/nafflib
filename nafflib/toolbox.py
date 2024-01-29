import numpy as np


#---------------------------------------
def henon_map(x,px,Q,n_turns):
    z_vec = np.nan*np.ones(n_turns) + 1j*np.nan*np.ones(n_turns)
    z_vec[0] = x - 1j*px
    for ii in range(n_turns-1):
        _z = z_vec[ii]

        z_vec[ii+1] =  np.exp(2*np.pi*1j* Q ) * (_z-1j/4 * (_z + np.conjugate(_z))**2)
    return np.real(z_vec),-np.imag(z_vec)
#---------------------------------------



#---------------------------------------
def parse_real_signal(amplitudes,frequencies,conjugate_tol=1e-10,to_pandas = False):
    
    A,Q = amplitudes,frequencies
    phasors = A*np.exp(2*np.pi*1j*Q)

    freq = []
    amp  = []

    for i in range(len(Q)-1):

        # Finding closest conjugate pair
        comparisons = phasors[i] + phasors
        pair_idx = np.argmin(np.abs(np.imag(comparisons)))
        pair_A = np.array([A[i],A[pair_idx]])
        pair_Q = np.array([Q[i],Q[pair_idx]])

        # Check if the pair is a complex conjugate, otherwise both freqs are kept
        if np.abs(np.diff(np.abs(pair_Q)))[0]>conjugate_tol:
            freq.append(Q[i])
            amp.append(A[i])
            continue

        # Creating avg amplitude and freq
        real = np.mean(np.real(pair_A))
        imag = np.mean(np.abs(np.imag(pair_A)))
        sign = np.sign(np.imag(pair_A))[pair_Q>=0]

        if pair_Q[0]==pair_Q[1]:
            # the pair is a copy of itself (DC signal case)
            freq.append(pair_Q[0])
            amp.append(real+1j*imag)
        else:
            # Complex conjugate found, adding only once
            if np.mean(np.abs(pair_Q)) not in freq:
                freq.append(np.mean(np.abs(pair_Q)))
                amp.append(2*(real+sign[0]*1j*imag))
    
    if to_pandas:
        import pandas as pd
        return pd.DataFrame({'amplitude':amp,'frequency':freq})
    else:
        return np.array(amp),np.array(freq)
#--------------------------------------- 



#---------------------------------------
def find_linear_combinations(frequencies,fundamental_tunes = [],max_harmonic_order = 10,to_pandas = False):
    """
    Categorisation of resonances. Returns the linear combinations of the fundamental tunes that are closest to the provided frequencies.
    This should be called after get_harmonics to have a list of frequencies.
    """

    assert len(fundamental_tunes) in [1,2,3], "need 1, 2 or 3 fundamental tunes (2D,4D,6D)"

    # Create a 3D array of all possible combinations of r_vec
    idx = max_harmonic_order
    if len(fundamental_tunes) == 1:
        r1,r2       = np.mgrid[-idx:idx+1, -idx:idx+1]
        r_vec       = [r1,r2]
    elif len(fundamental_tunes) == 2:
        r1,r2,r3    = np.mgrid[-idx:idx+1, -idx:idx+1,-idx:idx+1]
        r_vec       = [r1,r2,r3]
    else:
        r1,r2,r3,r4 = np.mgrid[-idx:idx+1, -idx:idx+1,-idx:idx+1,-idx:idx+1]
        r_vec       = [r1,r2,r3,r4]


    # Computing all linear combinations of r1*Qx + r2*Qy + r3*Qz + m
    Q_vec = fundamental_tunes + [1]
    all_combinations    = sum([_r*_Q for _r,_Q in zip(r_vec,Q_vec)])
    
    # Find the closest combination for each frequency
    r_values = []
    err = []
    for freq in frequencies:

        # Find the index of the closest combination
        closest_idx = np.unravel_index(np.argmin(np.abs(freq - all_combinations)), all_combinations.shape)

        # Get the corresponding values for r1,r2,r3,r4
        closest_combination = tuple(_r[closest_idx] for _r in r_vec)
        closest_value = all_combinations[closest_idx]

        r_values.append(closest_combination)
        err.append(np.abs(closest_value-freq))

    if to_pandas:
        import pandas as pd
        return pd.DataFrame({'resonance':r_values,'err':err,'freq':frequencies})
    else:
        return [tuple(_r) for _r in r_values],np.array(err),np.array(frequencies)
#---------------------------------------

#---------------------------------------
def generate_signal(amplitudes,frequencies,N):
    """
    Generate a signal with the provided amplitudes and frequencies over turns N.
    """

    if isinstance(amplitudes,(float,int)):
        amplitudes = [amplitudes]
    if isinstance(frequencies,(float,int)):
        frequencies = [frequencies]

    assert len(amplitudes) == len(frequencies), "Amplitudes and frequencies must have the same length"

    signal = sum([A*np.exp(1j*(2*np.pi*(Q)*N  )) for A,Q in zip(amplitudes,frequencies)])
    x  =  signal.real
    px = -signal.imag

    return x,px
#---------------------------------------


#---------------------------------------
def generate_pure_KAM(amplitudes,jklm,fundamental_tunes,N,return_frequencies=False):
    """
    Generate a signal with the provided amplitudes and frequencies over turns N.
    """

    if isinstance(amplitudes,(float,int)):
        amplitudes = [amplitudes]
    if isinstance(jklm,(float,int)):
        jklm = [jklm]

    assert len(amplitudes) == len(jklm), "amplitudes and jklm must have the same length"

    # Generating the signal
    Q1,Q2,Q3 = fundamental_tunes
    signal = sum([A*np.exp(2 * np.pi * 1j * (j*Q1 + k*Q2 + l*Q3 + m) * N) for A,(j,k,l,m) in zip(amplitudes,jklm)])
    x  =  signal.real
    px = -signal.imag

    if return_frequencies:
        return x,px,[j*Q1 + k*Q2 + l*Q3 + m for (j,k,l,m) in jklm]
    else:
        return x,px
#---------------------------------------