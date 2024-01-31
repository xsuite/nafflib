import matplotlib.pyplot as plt
import numpy as np
import nafflib

# ------------
num_turns = int(1e5)
x_points = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
px_point = 0.35 * x_points

# Henon parameters:
# ---
coupling = 0.1
Qx_list = [0.205, (3 - np.sqrt(5)) / 2]
Qx = Qx_list[1]
Qy = np.sqrt(2) - 1
# ------------


# Choosing number of harmonics
num_harmonics = 50

# Choosing turns to plot
N_start = 0
N_stop = 100

fig, axes = plt.subplots(1, 2, figsize=(18, 6))
plt.suptitle(f"Henon map 4D, (Qx,Qy) = ({Qx},{Qy})")
for x0, px0 in zip(x_points, px_point):

    # Tracking
    x, px, y, py = nafflib.henon_map_4D(
        x0, px0, x0, px0, Qx=Qx, Qy=Qy, coupling=coupling, num_turns=num_turns
    )
    # Saving in a dict
    dct = {"x": x, "px": px, "y": y, "py": py}

    # Plotting plane-by-plane
    for plane, ax in zip(["x", "y"], axes):
        plt.sca(ax)

        # Extracting harmonics
        A, Q = nafflib.harmonics(
            dct[plane], dct[f"p{plane}"], num_harmonics, window_order=4
        )

        # Reconstructing signal
        z_r, pz_r = nafflib.generate_signal(A, Q, np.arange(num_turns))

        # Plotting
        plt.plot(
            dct[plane][N_start:N_stop],
            dct[f"p{plane}"][N_start:N_stop],
            "o",
            mfc="none",
            color="C4",
            alpha=0.9,
        )
        plt.plot(z_r[N_start:N_stop], pz_r[N_start:N_stop], ".", color="C0", alpha=0.9)


# Adding labels
for plane, ax in zip(["x", "y", "zeta"], axes):
    plt.sca(ax)
    plt.plot(np.nan, np.nan, "o", mfc="none", color="C4", alpha=0.9, label="Henon")
    plt.plot(np.nan, np.nan, ".", color="C0", alpha=0.9, label="nafflib")
    plt.legend()

    plt.axis("equal")
    plt.xlabel(rf"$\tilde {plane}/\sqrt{{\varepsilon_{plane}}}$")
    plt.ylabel(rf"$\tilde p_{plane}/\sqrt{{\varepsilon_{plane}}}$")

    plt.xlim(-1, 1)
    plt.ylim(-1, 1)

plt.show()
