import numpy as np
from nafflib.naff import _fft_f0_estimate
import nafflib


# -----
# Henon map tune
Q_h = 0.2064898024701758
# -----
example_signals = {}
for x_start, label in zip([0.1, 0.3, 0.51], ["low_J", "mid_J", "high_J"]):
    example_signals[label] = nafflib.henon_map(x_start, 0.35 * x_start, Q_h, int(3e4))


def test_newton_method():

    for label, signal in example_signals.items():
        # Extracting signal
        x, px = signal
        z = x - 1j * px
        N = np.arange(len(z))
        A0, Q0 = nafflib.harmonics(
            z, num_harmonics=1, window_order=2, window_type="hann"
        )
        A0, Q0 = A0[0], Q0[0]
        # Windowing of the signal, testing order and N_turns
        # ---------------------
        amplitudes = []
        frequencies = []
        window_fun = nafflib.hann
        for order in [1, 2, 3, 4, 5]:
            for N_max, tol_Q, tol_A in zip(
                [1467, 2311, 5012, 10056, 26432],
                [1e-8, 1e-8, 1e-10, 1e-12, 1e-14],
                [1e-5, 1e-5, 1e-6, 1e-7, 1e-8],
            ):
                z_w = z[:N_max] * window_fun(N[:N_max], order=order)

                # Pretty much the script from fundamental_frequency, but here we vary the order and N_max
                # =======================================
                # Estimation of the main frequency with an FFT
                f0_est, resolution = _fft_f0_estimate(z_w)

                # Preparing the estimate for the Newton refinement method
                if f0_est >= 0.5:
                    f0_est = -(1.0 - f0_est)

                # Refinement of the tune calulation
                amplitude, f0 = nafflib.optimize.newton_method(
                    z_w, N[:N_max], freq_estimate=f0_est, resolution=resolution
                )
                # =======================================

                # print(f'{np.abs(Q0-f0):.5e}')
                assert np.isclose(
                    Q0, f0, atol=tol_Q, rtol=0
                ), f"Tune estimation failed with order = {order} and N_max = {N_max} for particle@{label}"

                # print(f'{np.abs(np.real(A0)-np.real(amplitude)):.5e}')
                assert np.isclose(
                    np.real(A0), np.real(amplitude), atol=tol_A, rtol=0
                ), f"Amp. estimation failed with order = {order} and N_max = {N_max} for particle@{label}"

                assert np.isclose(
                    np.imag(A0), np.imag(amplitude), atol=tol_A, rtol=0
                ), f"Amp. estimation failed with order = {order} and N_max = {N_max} for particle@{label}"
