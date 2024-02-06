# nafflib

**P. Belanger, K. Paraschou, G. Sterbini, G. Iadarola et al.**

A Python implementation of the Numerical Analysis of Fundamental Frequencies algorithm (NAFF [1-2]) from **J. Laskar**. This implementation uses a tailor-made optimizer (`nafflib.optimise.newton_method`, from **A. Bazzani, R. Bartolini & F. Schmidt**) to find the frequencies up to machine precision for tracking data. A Hann window is used to help with the convergence (`nafflib.windowing.hann`).

An insightful description of the NAFF algorithm is provided in the textbook by A. Wolski [3].

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
Q = nafflib.tune(z, window_order=1, window_type="hann")
Q = nafflib.tune(x, px, window_order=1, window_type="hann")
# --------------------------------------------------

# Using the position only:
# --------------------------------------------------
Q = nafflib.tune(x, window_order=1, window_type="hann")
# --------------------------------------------------
``` 

### Harmonics
 
Phase space trajectories (x,px),(y,py),(zeta,pzeta) are used to extract the spectral lines of the signal plane-by-plane with the `harmonics()` function. The number of harmonics is specified with the `num_harmonics` argument. Again, the function can be used with position only or position-momentum (preferred) information.

```python
# Let's assume the following signal
z = x - 1j * px

# The two following calls are equivalent
# --------------------------------------------------
spectrum = nafflib.harmonics(z, num_harmonics=5, window_order=1, window_type="hann")
spectrum = nafflib.harmonics(x, px, num_harmonics=5, window_order=1, window_type="hann")
# -> where spectrum = (amplitudes,frequencies)
# --------------------------------------------------

# From position only:
# --------------------------------------------------
spectrum = nafflib.harmonics(x, num_harmonics=5, window_order=1, window_type="hann")
# -> where spectrum = (amplitudes,frequencies)
# --------------------------------------------------
``` 

### Categorization of harmonics

For stable motion sufficiently close to integrable invariants of a conservative system, the frequencies are expected to come as a linear combinations of the fundamental tunes (3 for a 6D system). 

To properly study the harmonics of a system, the user should almost always try to unambigously identify the spectral lines, since very close lines can be mistaken for one another and **ordering them by amplitude will definitely lead to the wrong results**. 

Such a categorization of the spectral lines can be done for stable motion from a hamiltonian system like the LHC or any standard map by using the linear combination of fundamental frequencies as a unique ID to follow a given spectral line. See for example the `examples/nb_convergence.ipynb` notebook for such an approach and the `examples/ex_04_regularity_4D.py` for an example of the problems which can arise when ordering spectral lines by amplitude.


```python
# Let's assume the following dictionnary of phase space data:
data = {"x": x, "px": px, "y": y, "py": py, "zeta": zeta, "pzeta": pzeta}

# Let's extract the fundamental frequencies
Q_vec = [
    nafflib.tune(data[f"{plane}"], data[f"p{plane}"]) for plane in ["x", "y", "zeta"]
]

# Let's extract some harmonics
A, Q = nafflib.harmonics(x, px, num_harmonics=5)

# Let's find the linear combination of fundamental tunes (Q_vec)
# -----------------
# Note: max_harmonics_order might need to be set to a higher value to find
#       the proper linear combination of frequencies
# -----------------
categorization = nafflib.find_linear_combinations(
    Q, fundamental_tunes=Q_vec, max_harmonic_order=10
)
# -> where categorization = (r_vec,err,combined_frequency)

```





