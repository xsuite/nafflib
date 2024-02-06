import numpy as np
import nafflib


# -----
# Henon map tune
Q_h = 0.2064898024701758
# -----
example_signals = {}
for x_start, label in zip([0.1, 0.3, 0.51], ["low_J", "mid_J", "high_J"]):
    example_signals[label] = nafflib.henon_map(x_start, 0.35 * x_start, Q_h, int(3e4))


def test_hann():

    for label, signal in example_signals.items():

        # Extracting signal
        x, px = signal
        z = x - 1j * px
        N = np.arange(len(z))
        Q_tune = nafflib.tune(z)

        # Windowing of the signal, testing order and N_turns
        # ---------------------
        window_fun = nafflib.hann
        for order in [1, 2, 3, 4, 5]:
            for N_max in [1467, 2311, 5012, 10056, 26432]:
                z_w = z[:N_max] * window_fun(N[:N_max], order=order)
                f0_est, resolution = nafflib._fft_f0_estimate(z_w)
                # print(f'{np.abs(Q_tune-f0_est):.5e}')
                assert np.isclose(
                    Q_tune, f0_est, atol=1e-3, rtol=0
                ), f"Tune estimation failed with order = {order} and N_max = {N_max} for particle@{label}"
