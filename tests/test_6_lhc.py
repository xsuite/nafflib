import numpy as np
import nafflib


# =================================
num_turns = int(3e4)
files = [
    nafflib.__path__[0] + f"/../tests/data/LHC_particle_{s}_momentum_{i}sigma.csv"
    for s in ["on", "off"]
    for i in [1, 3, 5]
]


# -------------------
def read_csv(filename):
    filecontent = np.genfromtxt(
        filename,
        delimiter=",",
        skip_header=1,
        converters={col: lambda s: complex(s.decode()) for col in [1, 3, 5]},
        unpack=True,
    )
    data = {}
    data["Ax"] = filecontent[1]
    data["Qx"] = filecontent[2]
    data["Ay"] = filecontent[3]
    data["Qy"] = filecontent[4]
    data["Azeta"] = filecontent[5]
    data["Qzeta"] = filecontent[6]
    return data


# -------------------
example_signals = []
for file in files:
    data = read_csv(file)
    data["name"] = file
    for plane in ["x", "y", "zeta"]:
        z, pz = nafflib.generate_signal(
            data[f"A{plane}"], data[f"Q{plane}"], np.arange(num_turns)
        )
        data[f"{plane}"] = z
        data[f"p{plane}"] = pz
    example_signals.append(data)
# =================================


def test_map():
    n_harm = 100
    cmp_harm = 5
    N = np.arange(num_turns)
    for data in example_signals:
        Q_vec = [
            nafflib.tune(data[f"{plane}"], data[f"p{plane}"], window_order=4)
            for plane in ["x", "y", "zeta"]
        ]

        for plane in ["x", "y", "zeta"]:
            A, Q = data[f"A{plane}"], data[f"Q{plane}"]
            z_r, pz_r = data[f"{plane}"], data[f"p{plane}"]

            # Finding harmonics
            A_r, Q_r = nafflib.harmonics(
                z_r, pz_r, num_harmonics=n_harm, window_order=4, window_type="hann"
            )

            # Reconstructing signal
            z_rr, pz_rr = nafflib.generate_signal(A_r, Q_r, N)

            # Comparing tracking results
            assert np.allclose(
                z_r, z_rr, atol=1e-3, rtol=0
            ), f'X-Phasor tracking exceeded tolerance, in {data["name"]}, {plane} plane'
            assert np.allclose(
                pz_r, pz_rr, atol=1e-3, rtol=0
            ), f'PX-Phasor tracking exceeded tolerance, in {data["name"]}, {plane} plane'

            # Looking for matching harmonics
            # Sorting lines and compiling errors
            # -----------------------------------
            r, _, _ = nafflib.find_linear_combinations(
                Q[: int(2 * cmp_harm)], fundamental_tunes=Q_vec
            )
            r_found, _, _ = nafflib.find_linear_combinations(
                Q_r[:cmp_harm], fundamental_tunes=Q_vec, max_harmonic_order=10
            )

            errors_Q = []
            errors_A = []
            for res, _A, _Q in zip(r_found, A_r, Q_r):
                match_idx = r.index(res)
                errors_Q.append(np.abs(Q[match_idx] - _Q))
                errors_A.append(np.abs(np.abs(A[match_idx]) - np.abs(_A)))
            # -----------------------------------

            assert np.allclose(
                errors_Q, 0, atol=1e-11, rtol=0
            ), f'Q-tolerance not met, in {data["name"]}, {plane} plane'
            assert np.allclose(
                errors_A, 0, atol=1e-10, rtol=0
            ), f'A-tolerance not met, in {data["name"]}, {plane} plane'
