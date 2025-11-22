import numpy as np
import functools
import numba
import warnings


# Caching large grid of indices to look for linear combinations
#===========================================
@functools.lru_cache(maxsize=10) 
def N_GRID(dim,max_n,max_alias=None):
    assert max_n<100,'max_n should be smaller than 100'
    assert dim <= 3, 'dim should be smaller than 3'
    if max_alias is None:
        max_alias = max_n
    if dim == 1 : 
        _grid   = np.mgrid[ -max_n:max_n+1,-max_alias:max_alias+1]
    elif dim == 2 :
        _grid   = np.mgrid[ -max_n:max_n+1,-max_n:max_n+1,-max_alias:max_alias+1]
    elif dim == 3 :
        _grid   = np.mgrid[ -max_n:max_n+1,-max_n:max_n+1,-max_n:max_n+1,-max_alias:max_alias+1]

    # Flattening the grid
    n_grid      = [_n.ravel().astype(np.int16) for _n in _grid]

    # Preordering the grid with something similar to the LHC (just to speed up argsort!)
    dummy_tunes = [0.31,0.32,-2e-3,1]
    dummy_freqs = sum([_n * _Q for _n, _Q in zip(n_grid, dummy_tunes)])
    ordering    = np.argsort(dummy_freqs) 
    n_grid      = [_n[ordering] for _n in n_grid]
    
    return n_grid
#===========================================



# To find the closest match for each frequency (faster than np.argmin)
#===========================================
@numba.jit(nopython=True)
def _closest_match(freq_sorted,all_sorted,ordering):
    ii=0
    closest_idxs = []
    for freq in freq_sorted:
        found = False
        for jj in range(ii,len(all_sorted)):
            if all_sorted[jj]>freq:
                found = True
                if jj>0:
                    if abs(all_sorted[jj-1]-freq) < abs(all_sorted[jj]-freq):
                        ii = jj-1
                        closest_idxs.append(ordering[jj-1])
                    else:
                        ii = jj
                        closest_idxs.append(ordering[jj])
                else:
                    ii = jj
                    closest_idxs.append(ordering[jj])
                break
        if not found:
            closest_idxs.append(ordering[-1])

    return closest_idxs
#==================================================



def linear_combinations(frequencies,Qvec=[],max_n=10,max_alias=5,return_errors=False,warning_tol = 1e-8,to_pandas=False):
    """
    Identifies linear combinations of fundamental frequencies that best match a given set of frequencies.

    Parameters
    ----------
    frequencies : ndarray
        Array of frequencies to analyze.
    Qvec : list, optional
        List of fundamental tunes for resonance analysis. Length can be 1, 2, or 3.
    max_n : int, optional
        Maximum linear combination "n" value to consider (performance vs accuracy trade-off)
    max_alias : int, optional
        Maximum aliasing value to consider (performance vs accuracy trade-off)
    n_grid : list of ndarrays, optional
        List of arrays representing the possible combinations of the fundamental tunes. If not provided, it is generated. (for big arrays, it is better to provide it to avoid memory issues when looping)
    to_pandas : bool, optional
        If True, returns a pandas DataFrame; otherwise, returns lists and an array.

    Returns
    -------
    n_values,error,f_values : DataFrame or tuple
        L. combination values, errors, and corresponding frequencies, either as a DataFrame or separate data structures.
    """

    # Dealing with nans:
    #----------------------------
    if np.all(np.isnan(frequencies)):
        n_values = [tuple(np.nan for _ in range(len(Qvec))) for __ in range(len(frequencies))]
        errors   = [np.nan for __ in range(len(frequencies))]
        f_values = [np.nan for __ in range(len(frequencies))]
        if to_pandas:
            import pandas as pd
            return pd.DataFrame({"tuple": n_values, "error": errors, "frequency": f_values})
        elif return_errors:
            return n_values,np.array(errors),np.array(f_values)
        else:
            return n_values
    #----------------------------


    # Extracting N_GRID from the cache
    n_grid = N_GRID(dim=len(Qvec),max_n=max_n,max_alias=max_alias)

    # Listing all possible frequencie values
    all_values = sum([_n*_Q for _n,_Q in zip(n_grid,list(Qvec)+[1])])#np.sum([_n*_Q for _n,_Q in zip(n_grid,list(Qvec)+[1])],axis=0)


    # Sorting frequencies
    #----------------------------
    rev_sort = [idx for idx,_ in sorted(enumerate(np.argsort(frequencies)), key=lambda toflip: toflip[1])]
    freq_sort = np.sort(frequencies)
    #----------------------------

    # Cropping for shorter search (ignoring nans)
    #----------------------------
    freq_jump = np.max(np.diff(freq_sort[~np.isnan(freq_sort)]))
    condition = (all_values>np.nanmin(frequencies)-freq_jump) & (all_values<np.nanmax(frequencies)+freq_jump)
    #-
    n_grid      = [_n[condition] for _n in n_grid]
    all_values  = all_values[condition]
    ordering    = np.argsort(all_values)
    #----------------------------

    # Searching the best match for each frequency
    #----------------------------
    closest_found = _closest_match(freq_sort,all_values[ordering],ordering)
    #----------------------------


    # Unpack the results
    n_values = []
    f_values = []
    errors   = []
    #----------------------------
    for freq,closest in zip(frequencies,np.array(closest_found)[rev_sort]):
        # Forcing nan if freq == nan:
        if np.isnan(freq):
            n_tuple = tuple(np.nan for _n in n_grid)
            f_val   = np.nan
            err_val = np.nan
        else:

            # Get the corresponding values for n1,n2,n3,n4
            n_tuple = tuple(_n[closest] for _n in n_grid)
            
            # Get the corresponding frequency
            f_val   = all_values[closest]

            # Get the error
            err_val = np.abs(f_val - freq)

        # Append to the list
        n_values.append(n_tuple)
        f_values.append(f_val)
        errors.append(err_val)
    #----------------------------

    if np.any(np.array(errors)>warning_tol):
        warnings.warn('Large frequency error, consider increasing max_n or max_alias!')

    if to_pandas:
        import pandas as pd
        return pd.DataFrame({"tuple": n_values, "error": errors, "frequency": f_values})
    elif return_errors:
        return n_values,np.array(errors),np.array(f_values)
    else:
        return n_values


