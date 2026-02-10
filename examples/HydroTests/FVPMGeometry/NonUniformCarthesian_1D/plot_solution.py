import matplotlib.pyplot as plt
import swiftsimio as sw
import numpy as np
import argparse
import sys


def plot_mass_profile():
    parser = argparse.ArgumentParser(
        description="Plot mass profile from SWIFT snapshots."
    )
    parser.add_argument(
        "files", nargs="+", help="Snapshot files to plot (e.g. snapshot_0000.hdf5)"
    )
    parser.add_argument(
        "-o", "--output", type=str, default="mass_profile.png", help="Output filename."
    )
    parser.add_argument(
        "--var",
        type=str,
        default="metal_mass",
        choices=["metal_mass_fraction", "metal_mass", "gas_mass"],
        help="Variable to plot."
    )
    args = parser.parse_args()

    fig, ax = plt.subplots(figsize=(10, 6))

    for filename in args.files:
        try:
            data = sw.load(filename)
        except Exception as e:
            print(f"Could not load {filename}: {e}")
            continue

        x = data.gas.coordinates[:, 0].v
        m = data.gas.masses.v
        Z = data.gas.metal_mass_fractions.v
        m_Z = m*Z
        t = data.metadata.time.v

        # Sort by x to ensure the line plot connects neighbors correctly
        idx = np.argsort(x)

        ax.plot(
            x[idx],
            m_Z[idx],
            marker="o",
            markersize=3,
            linestyle="-",
            alpha=0.7,
            label=f"t={t:.3f}",
        )

    ax.axvline(x=0.5, color="k", linestyle="--", alpha=0.5)
    ax.set_xlabel("$x$")
    ax.set_ylabel("Metal mass fraction")
    ax.set_title("Advected Metal Mass Profile")
    ax.legend()
    ax.grid(True, alpha=0.2)

    plt.savefig(args.output, dpi=200)
    print(f"Saved plot to {args.output}")


if __name__ == "__main__":
    plot_mass_profile()
