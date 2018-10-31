from __future__ import print_function
import sys
sys.path.append('..')
import NAFFlib as NAFF
import numpy as np

np.random.seed(123456)
N=100
noise_rms=0.1
i = np.linspace(1,N,N)

q_true=6.24783351

x=np.empty_like(i,dtype=np.complex128)
print('                  True     frequency is Q_true = {0:.10f}'.format(q_true))

x  = 1*np.cos(2*np.pi*q_true*i) + 0j*np.sin(2*np.pi*q_true*i)
q = NAFF.get_tune(x)
print('(real           ) Estimated frequency is Q_hat = {0:.10f}'.format(q))

x  = 1*np.cos(2*np.pi*q_true*i) + 1j*np.sin(2*np.pi*q_true*i)
q = NAFF.get_tune(x)
print('(complex        ) Estimated frequency is Q_hat = {0:.10f}'.format(q))

x  = 1*np.cos(2*np.pi*q_true*i) + np.random.normal(0,noise_rms,N) +0j*np.sin(2*np.pi*q_true*i)
q = NAFF.get_tune(x)
print('(real    + noise) Estimated frequency is Q_hat = {0:.10f}'.format(q))

x  = 1*np.cos(2*np.pi*q_true*i) + np.random.normal(0,noise_rms,N) +1j*np.sin(2*np.pi*q_true*i) + 1j*np.random.normal(0,noise_rms,N)
q = NAFF.get_tune(x)
print('(complex + noise) Estimated frequency is Q_hat = {0:.10f}'.format(q))

print()
order = 1
interp = 1
print('Order = {0}, Interpolate = {1}'.format(order,interp))
x=np.empty_like(i,dtype=np.complex128)
print('                  True     frequency is Q_true = {0:.10f}'.format(q_true))

x  = 1*np.cos(2*np.pi*q_true*i) + 0j*np.sin(2*np.pi*q_true*i)
q = NAFF.get_tune(x, order, interp)
print('(real           ) Estimated frequency is Q_hat = {0:.10f}'.format(q))

x  = 1*np.cos(2*np.pi*q_true*i) + 1j*np.sin(2*np.pi*q_true*i)
q = NAFF.get_tune(x, order, interp)
print('(complex        ) Estimated frequency is Q_hat = {0:.10f}'.format(q))

x  = 1*np.cos(2*np.pi*q_true*i) + np.random.normal(0,noise_rms,N) +0j*np.sin(2*np.pi*q_true*i)
q = NAFF.get_tune(x, order, interp)
print('(real    + noise) Estimated frequency is Q_hat = {0:.10f}'.format(q))

x  = 1*np.cos(2*np.pi*q_true*i) + np.random.normal(0,noise_rms,N) +1j*np.sin(2*np.pi*q_true*i) + 1j*np.random.normal(0,noise_rms,N)
q = NAFF.get_tune(x, order, interp)
print('(complex + noise) Estimated frequency is Q_hat = {0:.10f}'.format(q))
