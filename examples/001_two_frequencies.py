from __future__ import print_function
import sys
sys.path.append('..')
sys.path.append('.')
import nafflib as NAFF

import numpy as np

np.random.seed(123456)
N=100000
noise_rms=0.8
i = np.linspace(1,N,N)

q_true=0.24783351

x=np.empty_like(i,dtype=np.complex128)
print('                  True     frequency is Q_true = {0:.10f}'.format(q_true))
print('                  True     frequency is Q_true = {0:.10f}'.format(q_true*1.1))
print('                  True     frequency is Q_true = {0:.10f}'.format(q_true*1.1**2))

x  =   1*np.cos(2*np.pi*q_true*i)              \
    +  0.5*np.cos(2*np.pi*1.1*q_true*i)        \
    +  0.3j*np.cos(2*np.pi*1.1**2*q_true*i)    \
    +  np.exp(1j*2*np.pi*0.44323*i)              
#    +  0j + np.random.normal(0,noise_rms,N)    \
#    +  1j*np.random.normal(0,noise_rms,N)      
NAFF.get_tunes(x,5)
#print('(real           ) Estimated frequency is Q_hat = {0:.10f}'.format(NAFF.get_tune(x)))

#x  = 1*np.cos(2*np.pi*q_true*i) + 1j*np.sin(2*np.pi*q_true*i)
#print('(complex        ) Estimated frequency is Q_hat = {0:.10f}'.format(NAFF.get_tune(x)))
#
#x  = 1*np.cos(2*np.pi*q_true*i) + np.random.normal(0,noise_rms,N) +0j*np.sin(2*np.pi*q_true*i)
#print('(real    + noise) Estimated frequency is Q_hat = {0:.10f}'.format(NAFF.get_tune(x)))
#
#x  = 1*np.cos(2*np.pi*q_true*i) + np.random.normal(0,noise_rms,N) +1j*np.sin(2*np.pi*q_true*i) + 1j*np.random.normal(0,noise_rms,N)
#print('(complex + noise) Estimated frequency is Q_hat = {0:.10f}'.format(NAFF.get_tune(x)))
