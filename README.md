# NAFFlib

A C library for the NAFF algorithm based on [NAFF_cpp](https://github.com/skostogl/NAFF_cpp) from Sofia Kostoglou.
Includes a Python wrapper.

***
```
$ make
```
to compile NAFFlib_c.so (python3) and NAFFlib2_c.so shared objects.

Alternatively, python3 and python2 support can be compiled separately by:
```
$ make py2
$ make py3
```

***
The NAFFlib module can be imported by:

```
import NAFFlib
```

#### Functions in the Python wrapper:
1. ```q = NAFFlib.get_tune(x, order, interpolation)```
where:  
- ```x``` is the (complex or real) input signal in the form of a one-dimensional non-empty numpy array.
- ```order``` is the value of the Hann's window order parameter to be used with 0 corresponding to no window. This variable is optional and its default value is equal to 2.
- ```interpolation``` is a boolean variable denoting whether the 7-point Newton-Cotes integration rule should be used. This variable is optional and its default value is equal to 0. It is not recommended that this is set 1.

returns the single positive frequency ```q``` that is dominant in the Fourier spectrum.


2. ```Q, A, B = NAFFlib.get_tunes(x, N, order, interpolation)```
where:  
- ```x``` is the (complex or real) input signal in the form of a one-dimensional non-empty numpy array.
- ```N``` is the number of frequencies to be found in the signal. This variable is optional and by default is set to 1.
- ```order``` is the value of the Hann's window order parameter to be used with 0 corresponding to no window. This variable is optional and its default value is equal to 2.
- ```interpolation``` is a boolean variable denoting whether the 7-point Newton-Cotes integration rule should be used. This variable is optional and its default value is equal to 0. It is not recommended that this is set 1.

returns three one-dimensional numpy arrays ```Q, A, B``` of size N with the first (being real-valued) containing the most dominant positive frequencies in the Fourier spectrum, the second containing the complex-valued amplitudes of the corresponding frequency and the third containing the complex-valued amplitudes of the negative of the corresponding frequency. It is recommended that this function is used with a real-valued input where the Fourier power spectrum is guaranteed to be an even function.


3. ```Q, A = NAFFlib.get_tunes(x, N, order, interpolation)```
where:  
- ```x``` is the (complex or real) input signal in the form of a one-dimensional non-empty numpy array.
- ```N``` is the number of frequencies to be found in the signal. This variable is optional and by default is set to 1.
- ```order``` is the value of the Hann's window order parameter to be used with 0 corresponding to no window. This variable is optional and its default value is equal to 2.
- ```interpolation``` is a boolean variable denoting whether the 7-point Newton-Cotes integration rule should be used. This variable is optional and its default value is equal to 0. It is not recommended that this is set 1.

returns two one-dimensional numpy arrays ```Q, A, B``` of size ```N``` with the first (being real-valued) containing the most dominant (positive or negative) frequencies in the Fourier spectrum and the second containing the complex-valued amplitudes of the corresponding frequency. It is recommended that this function is used with a complex-valued input where the Fourier power spectrum is not necessarily an even function. It should be emphasized that *positive and negative frequencies are treated separately*. 
