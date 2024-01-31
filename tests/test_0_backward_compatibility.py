import numpy as np
import nafflib


# -----
# Henon map tune
Q_h = 0.2064898024701758
# -----
example_signals = {}
for x_start, label in zip([0.1, 0.3, 0.51], ["low_J", "mid_J", "high_J"]):
    example_signals[label] = nafflib.henon_map(x_start, 0.35 * x_start, Q_h, int(3e4))


def test_isinstance():
    for label, signal in example_signals.items():
        # Extracting signal
        x, px = signal
        z = x - 1j * px

        # Testing _fft_f0_estimate
        output = nafflib._fft_f0_estimate(z)
        assert isinstance(output, tuple), "_fft_f0_estimate output is not a tuple"
        assert all(
            isinstance(item, float) for item in output
        ), "Not all elements in _fft_f0_estimate output are floats"

        # Testing fundamental_frequency
        amplitude, frequency = nafflib.fundamental_frequency(z)
        assert isinstance(
            amplitude, complex
        ), "Amplitude from fundamental_frequency is not complex"
        assert isinstance(
            frequency, float
        ), "Frequency from fundamental_frequency is not a float"

        # Testing naff
        amplitudes, frequencies = nafflib.naff(z, num_harmonics=10)
        assert isinstance(
            amplitudes, np.ndarray
        ), "Amplitudes from naff without pandas is not a numpy array"
        assert isinstance(
            frequencies, np.ndarray
        ), "Frequencies from naff without pandas is not a numpy array"
        assert all(
            isinstance(item, complex) for item in amplitudes
        ), "Not all amplitudes are complex"
        assert all(
            isinstance(item, float) for item in frequencies
        ), "Not all amplitudes are complex"

        # Testing tune
        frequency = nafflib.tune(z)
        assert isinstance(frequency, float), "Output from tune is not a float"

        # Testing harmonics
        amplitudes, frequencies = nafflib.harmonics(z, num_harmonics=10)
        assert isinstance(
            amplitudes, np.ndarray
        ), "Amplitudes from harmonics without pandas is not a numpy array"
        assert isinstance(
            frequencies, np.ndarray
        ), "Frequencies from harmonics without pandas is not a numpy array"
        assert all(
            isinstance(item, complex) for item in amplitudes
        ), "Not all amplitudes are complex"
        assert all(
            isinstance(item, float) for item in frequencies
        ), "Not all amplitudes are complex"

        # Testing harmonics position only
        amplitudes, frequencies = nafflib.harmonics(x, num_harmonics=10)
        assert isinstance(
            amplitudes, np.ndarray
        ), "Amplitudes from harmonics without pandas is not a numpy array"
        assert isinstance(
            frequencies, np.ndarray
        ), "Frequencies from harmonics without pandas is not a numpy array"
        assert all(
            isinstance(item, complex) for item in amplitudes
        ), "Not all amplitudes are complex"
        assert all(
            isinstance(item, float) for item in frequencies
        ), "Not all amplitudes are complex"

    # Testing multiparticle_tunes
    multiparticles = [signal[0] - 1j * signal[1] for signal in example_signals.values()]
    output = nafflib.multiparticle_tunes(multiparticles)
    assert isinstance(
        output, np.ndarray
    ), "Output from multiparticle_tunes is not a numpy array"
    assert (
        output.dtype == float
    ), "Elements in the output of multiparticle_tunes are not floats"


# Testing old package
# Should work on both old and revamped versions!
def test_NAFFlib():
    # ---------------------------------------
    def henon_map(x, px, Q, n_turns):
        z_vec = np.nan * np.ones(n_turns) + 1j * np.nan * np.ones(n_turns)
        z_vec[0] = x - 1j * px
        for ii in range(n_turns - 1):
            _z = z_vec[ii]

            z_vec[ii + 1] = np.exp(2 * np.pi * 1j * Q) * (
                _z - 1j / 4 * (_z + np.conjugate(_z)) ** 2
            )
        return np.real(z_vec), -np.imag(z_vec)

    # ---------------------------------------

    # Henon map tune
    # -----
    Q_h = 0.2064898024701758
    # -----
    example_signals = {}
    for x_start, label in zip([0.1, 0.3, 0.51], ["low_J", "mid_J", "high_J"]):
        example_signals[label] = henon_map(x_start, 0.35 * x_start, Q_h, int(1e4))

    for label, signal in example_signals.items():
        # Extracting signal
        x, px = signal
        z = x - 1j * px

        # Constructing exact phasors
        Q, A = nafflib.get_tunes_all(z, 20)
        _exact = sum(
            [
                _A * np.exp(1j * (2 * np.pi * (_Q) * np.arange(len(x))))
                for _A, _Q in zip(A, Q)
            ]
        )
        x_r = _exact.real
        px_r = -_exact.imag

        # Finding spectrum again
        Q_r, A_r = nafflib.get_tunes_all(x_r - 1j * px_r, 10)

        # Testing
        assert np.allclose(
            Q[0], Q_h, atol=1e-1, rtol=0
        ), f"Found wrong tune, particle@{label}"
        assert np.allclose(
            Q[:10], Q_r, atol=1e-8, rtol=0
        ), f"Expected tolerance not met, on Q, particle@{label}"
        assert np.allclose(
            A[:10], A_r, atol=1e-5, rtol=0
        ), f"Expected tolerance not met, on A, particle@{label}"

        # Checking get_tune function as well
        assert np.isclose(
            nafflib.get_tune(z), Q[0], atol=1e-10, rtol=0
        ), f"Expected tolerance not met from get_tune, particle@{label}"
