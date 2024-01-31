import numpy as np
import nafflib


# ========================
# -----
# Henon parameters
num_turns = int(1e5)
points = np.array([0.1, 0.3, 0.5])
slope = 0.35

Q_list = [0.2064898024701758, 0.3761365735491556, 0.1261960823639152]
# -----
example_signals = {}

for idx_map, Q_map in enumerate(Q_list):
    for idx_p, point in enumerate(points):
        x, px = nafflib.henon_map(point, slope * point, Q_map, num_turns)
        example_signals[f"map_{idx_map} : Q"] = Q_map
        example_signals[f"map_{idx_map} : part_{idx_p}"] = x, px
# ========================


def test_map():
    n_harm = 20
    N = np.arange(num_turns)
    for idx_map, Q_map in enumerate(Q_list):
        for idx_p, point in enumerate(points):
            Q_map = example_signals[f"map_{idx_map} : Q"]
            x, px = example_signals[f"map_{idx_map} : part_{idx_p}"]

            # Finding harmonics
            A, Q = nafflib.harmonics(
                x, px, num_harmonics=n_harm, window_order=2, window_type="hann"
            )

            # Finding linear combinations (high order for this 1D system!)
            r, err_r, _ = nafflib.find_linear_combinations(
                Q, fundamental_tunes=[Q[0]], max_harmonic_order=30
            )

            # Reconstructing signal and comparing
            x_r, px_r = nafflib.generate_signal(A, Q, N)

            # Checking results
            assert np.allclose(
                err_r, 0, atol=1e-9, rtol=0
            ), f"Linear combinations not found, map@{idx_map}, part@{idx_p}"
            assert np.allclose(
                x, x_r, atol=1e-2, rtol=0
            ), f"X-Phasor tracking exceeded tolerance, map@{idx_map}, part@{idx_p}"
            assert np.allclose(
                px, px_r, atol=1e-2, rtol=0
            ), f"PX-Phasor tracking exceeded tolerance, map@{idx_map}, part@{idx_p}"

            # Testing again on reconstruction
            A_r, Q_r = nafflib.harmonics(
                x_r, px_r, num_harmonics=n_harm, window_order=2, window_type="hann"
            )
            x_rr, px_rr = nafflib.generate_signal(A_r, Q_r, N)

            # Checking results
            assert np.allclose(
                x_r, x_rr, atol=1e-10, rtol=0
            ), f"Reconstruction, X-Phasor tracking exceeded tolerance, map@{idx_map}, part@{idx_p}"
            assert np.allclose(
                px_r, px_rr, atol=1e-10, rtol=0
            ), f"Reconstruction, PX-Phasor tracking exceeded tolerance, map@{idx_map}, part@{idx_p}"
