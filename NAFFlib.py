import os, ctypes
modulepath = os.path.dirname(os.path.abspath(__file__))
libpath = os.path.join(modulepath, 'NAFFlib.so')
NAFFlib = ctypes.CDLL(libpath)

import numpy as np
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
    hanning_order = ctypes.c_double(3.)
    tune = ctypes.c_double(0.)
    NAFFlib.pyget_f1(signal,N,hanning_order,ctypes.byref(tune))
    print tune
    return tune.value


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

q=0.12345
i = np.linspace(0,10000,10000)
x = np.cos(2.*np.pi*q*i,dtype=np.float64)
print '%.10f'%x[5]
print '%.7f'%get_tune(x)
print q,1-q
NAFFlib.hi();
