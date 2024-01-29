# sussix

A Python implementation for the frequency analysis tool SUSSIX (R. Bartolini, F. Schmidt et al.) used to study beam dynamics in particle accelerators,  https://cds.cern.ch/record/702438/. 

Sussix is a generic NAFF approach which benefits from a well-optimized solver (`sussix/optimise.py`) to find the frequencies up to machine precision in a lot of cases. A Hann window is used to help with the convergence (`sussix/windowing.py`)

An insightful description of the NAFF algorithm is provided in the textbook by A. Wolski, section 11.5: *A Numerical Method: Frequency Map Analysis* (https://www.worldscientific.com/doi/abs/10.1142/9781783262786_0011)

# Installation
```bash
pip install sussix
```

# Usage
Examples can be found in the `examples` folder. The spectrum of the data is computed using position-momentum data and the order of the Hann window can be specified by the user. Altough the algorithm works with position-only data, the use of position-momentum is preferred if possible. 



### Tune
The tune of a signal can be obtained from real or complexe signals as suchs:
```python
# Using the position only:
#--------------------------------------------------
sussix.get_tune(x,Hann_order=1)
#--------------------------------------------------

# Or position-momentum
#--------------------------------------------------
sussix.get_tune(x,px,Hann_order=1)
#--------------------------------------------------
``` 

### Spectrum
 
Phase space trajectories (x,px,y,py,zeta,pzeta) are used to extract the spectral lines of the signal with the `get_spectrum()` function. The number of harmonics is specified with the `number_of_harmonics` argument. Again, the function can be used with position only or position-momentum (preferred) information.

```python
# Individual spectrum

# From position only:
#--------------------------------------------------
sussix.get_spectrum(x,number_of_harmonics = 5,Hann_order = 1)
#--------------------------------------------------
# From position-momentum
#--------------------------------------------------
sussix.get_spectrum(x,px,number_of_harmonics = 5,Hann_order = 1)
#--------------------------------------------------

# Or can pass multiple canonical pairs with kwargs
#--------------------------------------------------
sussix.get_spectrum(x     = None,
                    px    = None,
                    y     = None,
                    py    = None,
                    zeta  = None,
                    pzeta = None,
                    number_of_harmonics = 5,
                    Hann_order          = 1)
# ->    returns df where df['x'] has the complex amplitudes (amplitude + phase) 
#       and frequencies in the x plane
#--------------------------------------------------


``` 

### Analysis

The indices (j,k,l,m) of the resonant frequencies can be found using `find_linear_combinations()` as done in the original SUSSIX code. 

The phase space trajectory can also be reconstructed from the frequency content using `generate_signal()`, which  simply sums the phasors optained from sussix. If the spectrum was obtained from position only, the user should discard the px output from `generate_signal`.
```python
x,px = sussix.generate_signal(spectrum.amplitude,spectrum.frequency,np.arange(int(1e4)))

# Or if spectrum comes from position only:
x,_ = sussix.generate_signal(spectrum.amplitude,spectrum.frequency,np.arange(int(1e4)))
```





