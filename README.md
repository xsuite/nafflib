# nafflib

**P. Belanger, K. Paraschou, G. Sterbini, G. Iadarola et al.**

A Python implementation of the Numerical Analysis of Fundamental Frequencies algorithm (NAFF [1-2]) from **J. Laskar**. This implementation uses a tailor-made optimizer (`nafflib.optimise.newton_method`, from **A. Bazzani, R. Bartolini & F. Schmidt**) to find the frequencies up to machine precision for tracking data. A Hann window is used to help with the convergence (`nafflib.windowing.hann`).

An in-depth description of the NAFF algorithm is provided in P. Belanger's Ph.D. thesis. Additionnal references are given below.

[1] J. Laskar, Introduction to Frequency Map Analysis. http://link.springer.com/10.1007/978-94-011-4673-9_13  
[2] J. Laskar et al., The Measure of Chaos by the Numerical Analysis of the Fundamental Frequencies. Application to the Standard Mapping. https://doi.org/10.1016/0167-2789(92)90028-L   
[3] A. Wolski,Beam Dynamics in High Energy Particle Accelerators, Section 11.5: *A Numerical Method: Frequency Map Analysis*. https://www.worldscientific.com/doi/abs/10.1142/9781783262786_0011 
# Installation
```bash
pip install nafflib
```

# Usage
Examples can be found in the `examples` folder. The harmonics of the data are computed using position-momentum data and the order of the window can be specified by the user. Altough the algorithm works with position-only data, the use of position-momentum is preferred if possible. 



### Tune
The tune of a signal can be obtained from real or complexe signals as suchs:
```python
# Let's assume the following signal
z = x - 1j * px

# The two following calls are equivalent
# --------------------------------------------------
Q = nafflib.tune(z, window_order=2, window_type="hann")
Q = nafflib.tune(x, px, window_order=2, window_type="hann")
# --------------------------------------------------

# Using the position only:
# --------------------------------------------------
Q = nafflib.tune(x, window_order=2, window_type="hann")
# --------------------------------------------------
``` 

### Harmonics
 
Phase space trajectories (x,px),(y,py),(zeta,pzeta) are used to extract the spectral lines of the signal plane-by-plane with the `harmonics()` function. The number of harmonics is specified with the `num_harmonics` argument. Again, the function can be used with position only or position-momentum (preferred) information.

```python
# Let's assume the following signal
z = x - 1j * px

# The two following calls are equivalent
# --------------------------------------------------
amplitudes,frequencies = nafflib.harmonics(z, num_harmonics=5, window_order=2, window_type="hann")
amplitudes,frequencies = nafflib.harmonics(x, px, num_harmonics=5, window_order=2, window_type="hann")
# --------------------------------------------------

# From position only:
# --------------------------------------------------
amplitudes,frequencies = nafflib.harmonics(x, num_harmonics=5, window_order=2, window_type="hann")
# --------------------------------------------------
``` 

### Categorization of harmonics

For KAM trajectories, the harmonics are expected to come as a linear combinations of the fundamental frequencies (3 for a 6D system). 

To properly study the harmonics of a system, it is helpful to retreive these linear combinations as done below: 


```python
# Let's assume the following dictionnary of phase space data:
data = {"x": x, "px": px, "y": y, "py": py, "zeta": zeta, "pzeta": pzeta}

# Let's extract the fundamental frequencies
Q_vec = [
    nafflib.tune(data[f"{plane}"], data[f"p{plane}"]) for plane in ["x", "y", "zeta"]
]

# Let's extract some harmonics
Ax, Qx = nafflib.harmonics(x, px, num_harmonics=5)

# Indexing harmonics
#============================================================================
max_n       = 90  #(high numbers needed in 2D..)
max_alias   = 50
warning_tol = np.inf #Disable warnings
#-------------------------------------
nx      = nafflib.linear_combinations(Qx,   Qvec    = Qvec,
                                            max_n   = max_n,
                                            max_alias   = max_alias, 
                                            warning_tol = warning_tol)
#============================================================================
```





