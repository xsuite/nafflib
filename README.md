# NAFFlib
Authors:  
* **Sofia Kostoglou**  
* **Konstantinos Paraschou**
* **Dario Pellegrini** 

Maintained by **Konstantinos Paraschou**.

NAFFlib is a library written in C and wrapped in 
Python which includes an implementation of the NAFF 
algorithm. 
[[1]](https://www.sciencedirect.com/science/article/pii/001910359090084M) 
[[2]](http://jacow.org/ipac2017/papers/thpab044.pdf)

## Installation
### Automatic
```
$ pip install NAFFlib
```
Depending on the OS, you might need to 
```
$ pip3 install NAFFlib
```
for Python3 support.

### Local Installation
```
$ git clone git@github.com:PyCOMPLETE/NAFFlib.git
$ cd NAFFlib
$ make
$ cd ..
```
to compile NAFFlib_c.so (Python3) and NAFFlib2_c.so (Python2) shared objects.
If either of Python2 or Python3 is not supported in the OS, the line
```
$ make
```
should be replaced with either
```
$ make py2
```
for Python2 support or
```
$ make py3
```
for Python3 support.

***
The NAFFlib module can be imported by:

```
import NAFFlib
```

### Functions in NAFFlib:
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

returns two one-dimensional numpy arrays ```Q, A``` of size ```N``` with the first (being real-valued) containing the most dominant (positive or negative) frequencies in the Fourier spectrum and the second containing the complex-valued amplitudes of the corresponding frequency. It is recommended that this function is used with a complex-valued input where the Fourier power spectrum is not necessarily an even function. It should be emphasized that *positive and negative frequencies are treated separately*. 

4. ```Q = NAFFlib.multiparticle_tunes(x, order, interpolation)```
where:  
- ```x``` is an array of (complex or real) input signals in the form of a two-dimensional non-empty numpy array. The first axis should correspond to the id of each different track while the second axis should correspond to the turn number.
- ```order``` is the value of the Hann's window order parameter to be used with 0 corresponding to no window. This variable is optional and its default value is equal to 2.
- ```interpolation``` is a boolean variable denoting whether the 7-point Newton-Cotes integration rule should be used. This variable is optional and its default value is equal to 0. It is not recommended that this is set 1.

returns a one-dimensional numpy array ```Q``` of size ```len(x)``` which contains the single most dominant frequency of the different tracks. 
