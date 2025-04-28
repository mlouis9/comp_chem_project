#!/usr/bin/env python3
"""Plot KE, PE, temperature, and conserved energy from a CP2K *.ener file."""

import sys
import pandas as pd
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
def load_ener(path):
    """Return a tidy DataFrame from a CP2K *.ener file."""
    cols = [
        "Step", "Time_fs", "Kin_au", "Temp_K",
        "Pot_au", "Cons_au", "Wall_s",
    ]
    df = pd.read_csv(
        path,
        delim_whitespace=True,
        comment="#",          # skip header lines that start with '#'
        names=cols,
        engine="python",
    )
    df["Time_ps"] = df["Time_fs"] / 1000.0   # fs → ps for nicer x-axis
    return df

# ----------------------------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("usage: python plot_energy.py path/to/file.ener")
    ener_file = sys.argv[1]

    df = load_ener(ener_file)

    fig, axes = plt.subplots(4, 1, sharex=True, figsize=(6, 8))

    axes[0].plot(df["Time_ps"], df["Kin_au"])
    axes[0].set_ylabel("KE (a.u.)")

    axes[1].plot(df["Time_ps"], df["Pot_au"])
    axes[1].set_ylabel("PE (a.u.)")

    axes[2].plot(df["Time_ps"], df["Temp_K"])
    axes[2].axhline(1250, color="red", linestyle="--", label="Target T")
    axes[2].set_ylabel("T (K)")

    axes[3].plot(df["Time_ps"], df["Cons_au"])
    axes[3].set_ylabel("Cons. E (a.u.)")
    axes[3].set_xlabel("Time (ps)")

    for ax in axes:
        ax.grid(alpha=0.3)

    fig.tight_layout()
    plt.savefig("energy_time_series.png", dpi=300)
    plt.legend()
    plt.show()
