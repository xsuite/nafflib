import numpy as np
import nafflib


def cauchy(x, loc, scale):
    # from https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.cauchy.html
    return 1 / (np.pi * scale * (1 + ((x - loc) / scale) ** 2))


def dense_spectrum(Q0, param):
    # Creating dummy signal with dense spectral lines
    # ==============================================
    N = np.arange(int(1e5))
    Q0 = 0.31025793875089835
    Qs = 0.002
    dQ = Qs / param
    Jx = 0.5 * (10**2)

    n_bands_Qs = 1
    n_bands_dQ = 5
    i, j = np.arange(-n_bands_dQ, n_bands_dQ + 1), np.arange(
        -n_bands_Qs, n_bands_Qs + 1
    )
    Ai, Aj = cauchy(i / np.max(i), 0, 0.05), cauchy(j / np.max(j), 0, 0.05)
    Ai, Aj = Ai / np.max(Ai), Aj / np.max(Aj)

    amplitudes = np.array(
        [
            [(np.sqrt(2 * Jx) * _Ai * _Aj) for _i, _Ai in zip(i, Ai)]
            for _j, _Aj in zip(j, Aj)
        ]
    ).flatten()
    frequencies = np.array(
        [[Q0 + _j * Qs + _i * dQ for _i, _Ai in zip(i, Ai)] for _j, _Aj in zip(j, Aj)]
    ).flatten()

    Q_vec = [Q0, Qs, dQ]

    return (
        amplitudes,
        frequencies,
        Q_vec,
        nafflib.generate_signal(amplitudes, frequencies, N),
    )
    # ==============================================


# -----
# Dummy signal main frequency
Q0 = 0.31025793875089835
# -----
example_signals = {}
for num, label in zip([1, 2, 3], ["param_1", "param_2", "param_3"]):
    # set of irrationnal numbers for spectrum spacing
    param = 5 * num * (np.pi / 3) ** num
    A, Q, Q_vec, signal = dense_spectrum(Q0, param)
    example_signals[label] = signal
    example_signals[f"{label}:A"] = A
    example_signals[f"{label}:Q"] = Q
    example_signals[f"{label}:Q_vec"] = Q_vec


def test_dense_spectrum():

    for num, label in zip([1, 2, 3], ["param_1", "param_2", "param_3"]):

        # Extracting signal
        fundamental_tunes = [
            Q0,
        ]
        x, px = example_signals[label]
        A, Q, Q_vec = (
            example_signals[f"{label}:A"],
            example_signals[f"{label}:Q"],
            example_signals[f"{label}:Q_vec"],
        )

        # nafflib
        A_found, Q_found = nafflib.harmonics(
            x, px, num_harmonics=len(Q), window_order=2, window_type="hann"
        )

        # Sorting lines and compiling errors
        # -----------------------------------
        r, _, _ = nafflib.find_linear_combinations(Q, fundamental_tunes=Q_vec)
        r_found, _, _ = nafflib.find_linear_combinations(
            Q_found, fundamental_tunes=Q_vec
        )

        errors_Q = []
        errors_A = []
        for res, _A, _Q in zip(r, A, Q):
            found_idx = r_found.index(res)
            errors_Q.append(Q_found[found_idx] - _Q)
            errors_A.append(np.abs(A_found[found_idx]) - np.abs(_A))
        # -----------------------------------

        # Comparing results
        assert np.allclose(
            errors_Q, 0, atol=1e-9, rtol=0
        ), f"Frequency mismatch for param@{label}"
        assert np.allclose(
            errors_A, 0, atol=1e-5, rtol=0
        ), f"Amplitude mismatch for param@{label}"
