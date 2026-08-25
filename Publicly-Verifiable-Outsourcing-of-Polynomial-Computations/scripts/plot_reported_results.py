#!/usr/bin/env python3
"""Replot numerical coordinates reported in the paper.

The CSV files are transcribed from the PGFPlots coordinates in the paper. They
reproduce the reported plots; they are not raw timing logs from a fresh C run.
"""
from pathlib import Path
import csv
import os

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / "build" / "matplotlib"))

import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "reported"
OUT = ROOT / "results" / "plots"
OUT.mkdir(parents=True, exist_ok=True)
def read_csv(name):
    with (DATA / name).open(newline="") as f: return list(csv.DictReader(f))
def plot_groups(rows, group_keys, xkey, ykey, xlabel, ylabel, outname, title):
    fig, ax = plt.subplots()
    groups = {}
    for r in rows: groups.setdefault(tuple(r[k] for k in group_keys), []).append(r)
    for key, vals in groups.items():
        vals.sort(key=lambda r: float(r[xkey]))
        label = ", ".join(f"{k}={v}" for k, v in zip(group_keys, key))
        ax.plot([float(v[xkey]) for v in vals], [float(v[ykey]) for v in vals], marker="o", label=label)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.set_title(title); ax.legend(); fig.tight_layout(); fig.savefig(OUT / outname, dpi=200); plt.close(fig)
plot_groups(read_csv("figure2_low_degree.csv"), ["scheme"], "m", "log2_Rt_star", "m", "log2(Rt*)", "figure2_low_degree.png", "Figure 2 reported data")
plot_groups(read_csv("figure3_thresholds.csv"), ["scheme","t"], "m", "log2_Rt_star", "m", "log2(Rt*)", "figure3_thresholds.png", "Figure 3 reported data")
rows4 = read_csv("figure4_high_degree.csv")
for d in sorted({r["d"] for r in rows4}, key=int):
    plot_groups([r for r in rows4 if r["d"] == d], ["scheme"], "m", "log2_Rt_star", "m", "log2(Rt*)", f"figure4_d{d}.png", f"Figure 4 reported data, d={d}")
plot_groups(read_csv("figure5_gamma1_vs_pi3.csv"), ["d","t"], "m", "Tc_Gamma1_over_Tc_Pi3", "m", "Tc(Gamma1)/Tc(Pi3)", "figure5_gamma1_vs_pi3.png", "Figure 5 reported data")
print(f"Wrote plots to {OUT}")
