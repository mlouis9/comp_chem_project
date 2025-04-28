#!/usr/bin/env python3
"""
Compute g(r) for U–U, Cl–Cl and U–Cl from a CP2K XYZ trajectory.

Usage
-----
python rdf_ucl3.py trajectory.xyz 12.38
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

# ------------- user-tunable parameters ------------- #
DR          = 0.02      # Å, radial bin width
START_SKIPS = 100      # discard this many initial frames
STRIDE      = 1         # analyse every STRIDE-th frame afterwards
R_MAX       = None      # if None → 0.5*box_length
# --------------------------------------------------- #

# ---------- minimalist I/O ---------- #
def read_xyz(path):
    """Yield (symbols[N], coords[N,3]) from an XYZ trajectory."""
    with open(path) as fh:
        while True:
            headline = fh.readline()
            if not headline:
                break
            n = int(headline)
            fh.readline()                       # comment line
            s, xyz = [], []
            for _ in range(n):
                sym, *nums = fh.readline().split()
                s.append(sym)
                xyz.append([float(x) for x in nums])
            yield np.array(s), np.array(xyz, float)

def minimum_image(vecs, box_len):
    return vecs - box_len * np.rint(vecs / box_len)

# ---------- RDF driver ---------- #
def accumulate_histograms(traj_iter, box_len, dr, r_max,
                          start_skips, stride):
    r_max  = r_max or box_len * 0.5
    nbins  = int(np.floor(r_max / dr))
    volume = box_len ** 3
    shell_vol = 4*np.pi/3 * (
        (np.arange(1, nbins+1) * dr)**3 -
        (np.arange(nbins)     * dr)**3
    )
    hists   = {k: np.zeros(nbins) for k in ("UU", "ClCl", "UCl")}
    nframes = defaultdict(int)

    # Will be filled when the first analysed frame is met
    n_u = n_cl = None

    for fidx, (sym, pos) in enumerate(traj_iter):
        if fidx < start_skips or (fidx - start_skips) % stride:
            continue

        idx_u  = np.where(sym == "U")[0]
        idx_cl = np.where(sym == "Cl")[0]

        if n_u is None:
            n_u, n_cl = len(idx_u), len(idx_cl)

        species_pairs = {
            "UU"  : (idx_u.copy(),  idx_u.copy(),  True),
            "ClCl": (idx_cl.copy(), idx_cl.copy(), True),
            "UCl" : (idx_u.copy(),  idx_cl.copy(), False),
        }

        for tag, (iset, jset, identical) in species_pairs.items():
            if identical:
                for i in iset:
                    js = jset[jset > i]        # avoid double counts
                    vecs = minimum_image(pos[js] - pos[i], box_len)
                    d = np.linalg.norm(vecs, axis=1)
                    bins = (d // dr).astype(int)
                    np.add.at(hists[tag], bins[bins < nbins], 1)
            else:                               # U–Cl: full cartesian product
                vecs = minimum_image(
                    pos[np.ix_(jset)] - pos[iset][:, None], box_len
                ).reshape(-1, 3)
                d = np.linalg.norm(vecs, axis=1)
                bins = (d // dr).astype(int)
                np.add.at(hists[tag], bins[bins < nbins], 1)

            nframes[tag] += 1

    if n_u is None:     # nothing analysed
        raise RuntimeError(
            "No frame left after START_SKIPS / STRIDE. "
            "Reduce START_SKIPS or check the trajectory length."
        )

    # --- normalisation --- #
    g = {}
    npairs = {
        "UU"  : n_u * (n_u - 1) / 2,
        "ClCl": n_cl * (n_cl - 1) / 2,
        "UCl" : n_u * n_cl,
    }
    for tag, hist in hists.items():
        norm = npairs[tag] * nframes[tag] * shell_vol / volume
        g[tag] = hist / norm

    r = (np.arange(nbins) + 0.5) * dr
    return r, g

# ---------- main ---------- #
if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit("usage: python rdf_ucl3.py traj.xyz box_length_Å")

    traj_path, Lbox = sys.argv[1], float(sys.argv[2])

    r, g = accumulate_histograms(
        read_xyz(traj_path),
        Lbox,
        DR,
        R_MAX,
        START_SKIPS,
        STRIDE,
    )

    plt.figure(figsize=(6, 4))
    plt.plot(r, g["UU"],   label=r"$g_{\mathrm{UU}}(r)$")
    plt.plot(r, g["ClCl"], label=r"$g_{\mathrm{ClCl}}(r)$")
    plt.plot(r, g["UCl"],  label=r"$g_{\mathrm{UCl}}(r)$")
    plt.xlim(0, r[-1])
    plt.xlabel(r"$r\;[\mathrm{\AA}]$")
    plt.ylabel(r"$g(r)$")
    plt.legend()
    plt.tight_layout()
    plt.savefig("rdf_ucl3.png", dpi=300)
    plt.show()
