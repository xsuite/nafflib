import os, ctypes
modulepath = os.path.dirname(os.path.abspath(__file__))
libpath = os.path.join(modulepath, 'NAFFlib.so')
NAFFlib = ctypes.CDLL(libpath)

import numpy as np
import math
import matplotlib.pyplot as plt

class NAFFlib_struct(ctypes.Structure):
    _fields_ = [("a",ctypes.c_int),
                ("b",ctypes.c_int),
                ("win",ctypes.c_char)]

class c_double_complex(ctypes.Structure):
    _fields_ = [("real", ctypes.c_double),
                ("imag", ctypes.c_double)]

    def __init__(self, z):
        self.real = z.real
        self.imag = z.imag
    
    def to_complex(self):
        return self.real + 1.j* self.imag

def complex_c_pointer(nparray):
    # The pointer will not point to the original array
    # but rather to a copy of the original.
    temp_arr = c_double_complex*len(nparray)
    return temp_arr(*(c_double_complex(z) for z in nparray))

def get_tune(x):
    signal = complex_c_pointer(x)
    N = ctypes.c_int(len(x))
    hanning_order = ctypes.c_double(2.)
    tune = ctypes.c_double(0.)
    NAFFlib.pyget_f1(signal,N,hanning_order,ctypes.byref(tune))
    #print tune
    return tune.value

def get_tune2(x):
    signal = complex_c_pointer(x)
    N = ctypes.c_int(len(x))
    hanning_order = ctypes.c_double(3.)
    tune = ctypes.c_double(0.)
    fft_estimate = 1.*np.argmax(np.fft.fft(x)[:len(x)/2])/len(x)
    #print fft_estimate
    NAFFlib.py_f1(signal,N,hanning_order,ctypes.c_double(fft_estimate), ctypes.byref(tune))
    #print tune
    return tune.value

def get_tunes(x, n):
    signal = complex_c_pointer(x)
    N = ctypes.c_int(len(x))
    hann_order = ctypes.c_double(3.)
    tunes = np.empty(n,dtype=np.float64)
    t_amps = np.empty(n,dtype=np.float64)
    t_namps = np.empty(n,dtype=np.float64)
    amps = complex_c_pointer(t_amps) 
    namps = complex_c_pointer(t_namps) 
    nfreqs = ctypes.c_int(n)
    NAFFlib.get_f_neg(signal, N, hann_order, ctypes.c_void_p(tunes.ctypes.data), amps, namps, nfreqs)
    for i in range(n):
        print tunes[i], abs(amps[i].to_complex()-namps[i].to_complex()), amps[i].to_complex(), namps[i].to_complex()

#a = ctypes.c_int(3)
#data = NAFFlib_struct(5,3,"h")
#NAFFlib.hello(ctypes.byref(a))
#print a.value
#
#NAFFlib.hello2(ctypes.byref(data),a)
#
#A = np.linspace(0,10,3)
#NAFFlib.hi3( A.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), len(A))
#print sum(A)
#
#A = np.array([1+6j,2+7j,3+8j,4+9j,5+10j])
#cA = complex_c_pointer(A)
#NAFFlib.test_func(complex_c_pointer(A),len(A))
#print A
#
#A = np.empty(101)
#cA = complex_c_pointer(A)
#order = 3.
#NAFFlib.taylor_window(cA, len(A),ctypes.c_double(order))
#
#A = np.array([ z.to_complex().real for z in cA])
#plt.plot(A)
##plt.show()

#q=0.1234567891
#i = np.linspace(0,30000,30001)
#x = np.cos(2.*np.pi*q*i,dtype=np.float64)
##x = np.cos(2.*np.pi*q*i,dtype=np.complex128)
##print('{0:.10f}'.format(x[5]) )
##print(np.pi)
#print '%.10f'%get_tune(x)
#print '%.10f'%get_tune2(x)
##print q,1-q
##NAFFlib.hi();
##
##asd = ctypes.c_double(1.234567890123456)
##print asd.value
