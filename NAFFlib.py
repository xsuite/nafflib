import os, ctypes
modulepath = os.path.dirname(os.path.abspath(__file__))
libpath = os.path.join(modulepath, 'cNAFF.so')
cNAFF = ctypes.CDLL(libpath)

import numpy as np
##from cNAFF import cNAFF_hello
import matplotlib.pyplot as plt

class cNAFF_struct(ctypes.Structure):
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

a = ctypes.c_int(3)
data = cNAFF_struct(5,3,"h")
cNAFF.hello(ctypes.byref(a))
print a.value

cNAFF.hello2(ctypes.byref(data),a)

A = np.linspace(0,10,3)
cNAFF.hi3( A.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), len(A))
print sum(A)

A = np.array([1+6j,2+7j,3+8j,4+9j,5+10j])
cA = complex_c_pointer(A)
cNAFF.test_func(complex_c_pointer(A),len(A))
print A

A = np.empty(101)
cA = complex_c_pointer(A)
order = 3.
cNAFF.taylor_window(cA, len(A),ctypes.c_double(order))

A = np.array([ z.to_complex().real for z in cA])
plt.plot(A)
#plt.show()

cNAFF.hi();
