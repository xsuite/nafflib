import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import nafflib


# Simple FFT Wrapper
# -------------------
def get_FFT(x):
    x = np.array(x)
    turns = np.arange(len(x))

    # Cropping signal to closest power of 2
    Nt = len(x)
    crop_at = 2 ** int(np.log2(Nt))

    spectrum = np.fft.fft(x[:crop_at])
    freq = np.fft.fftfreq(turns[:crop_at].shape[-1])

    idx = np.argmax(np.abs(spectrum))
    Qx = freq[idx]
    return freq[freq > 0], np.abs(spectrum)[freq > 0]


# -------------------

# ===========================
# Simple signal: https://en.wikipedia.org/wiki/Phase_modulation, should come as a sum of harmonic peaks
Qx, Qzeta = [0.31, 0.0018]
N = np.arange(int(1e5))

a = 0.2
phase = a * np.sin(2 * np.pi * Qzeta * N)
# -------------------
signal = 10 * np.exp(2 * np.pi * 1j * (Qx * N + phase))
x, px = np.real(signal), -np.imag(signal)
# ===========================


# Finding harmonics
A, Q = nafflib.harmonics(x, px, window_order=4, num_harmonics=10)

# Plotting
# NOTE: FFT amplitude is NOT expected to be the same as phasor amplitude
# =====================
plt.figure(figsize=(15, 4))
freq, s_fft = get_FFT(x)
for _Q in Q:
    plt.axvline(_Q - Qx, color="k", alpha=0.2)
plt.plot(freq - Qx, s_fft / np.max(s_fft), label="fft")
plt.plot(
    Q - Qx,
    np.abs(A) / np.max(np.abs(A)),
    "o",
    mfc="none",
    color="k",
    alpha=0.5,
    label="phasor expansion",
)
plt.xlim(-10 * Qzeta, 10 * Qzeta)
plt.yscale("log")
plt.legend()
plt.ylabel("Normalized amplitude")
plt.xlabel("Freq (Q - Qx)")
plt.show()
# =====================
