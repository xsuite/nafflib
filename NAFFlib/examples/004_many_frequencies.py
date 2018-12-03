from __future__ import print_function
import sys
sys.path.append('..')
import NAFFlib as NAFF

import numpy as np
np.random.seed(123456)
N=100000
noise_rms=0.8
i = np.linspace(1,N,N)

q_true=0.24783351
q_exp = 0.44323

x=np.empty_like(i,dtype=np.complex128)
print('                  True     frequency is Q_true = {0:.10f}'.format(q_exp))
print('                  True     frequency is Q_true = {0:.10f}'.format(q_true))
print('                  True     frequency is Q_true = {0:.10f}'.format(q_true*1.1))
print('                  True     frequency is Q_true = {0:.10f}'.format(q_true*1.1**2))

x  =   1*np.cos(2*np.pi*q_true*i)              \
    +  0.5*np.cos(2*np.pi*1.1*q_true*i)        \
    +  0.3j*np.cos(2*np.pi*1.1**2*q_true*i)    \
    +  np.exp(1j*2*np.pi*q_exp*i)              
#    +  0j + np.random.normal(0,noise_rms,N)    \
#    +  1j*np.random.normal(0,noise_rms,N)      
#q = NAFF.get_tunes(x,5,3)
q, Ap, An = NAFF.get_tunes(x,5)
print('Tune, (Positive frequency\'s amplitude), (Negative frequency\'s amplitude)') 
for i in range(len(q)):
    print('{0:.6f}, {1:+.6f}{2:+.6f}i, {3:+.6f}{4:+.6f}i'.format(q[i], Ap[i].real, Ap[i].imag, An[i].real, An[i].imag))
