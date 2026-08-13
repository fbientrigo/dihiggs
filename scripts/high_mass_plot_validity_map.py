#!/usr/bin/env python3
"""Diagnostic scatter plot: mH2 vs mA (Delta_heavy annotated), accepted vs
rejected attempted points, colored by rejection_stage. Reads
results/high_mass_valid_points/high_mass_attempted_points.csv (written by
scripts/high_mass_physical_point_search.py) and writes
mass_validity_map.png/.pdf into the same directory. Not a publication figure;
a diagnostic for this task only.
"""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = ROOT / "results/high_mass_valid_points"
ATTEMPTED_CSV = RESULTS_DIR / "high_mass_attempted_points.csv"

STAGE_COLORS = {
    "accepted": "#2ca02c",
    "positivity": "#d62728",
    "unitarity": "#ff7f0e",
    "perturbativity": "#9467bd",
    "numerical": "#8c564b",
    "construction": "#7f7f7f",
    "cli_invocation": "#000000",
}


def main():
    rows = []
    with ATTEMPTED_CSV.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            rows.append(row)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left panel: mH2 vs mA, colored by rejection_stage (accepted=green on top)
    ax = axes[0]
    by_stage = {}
    for row in rows:
        try:
            mh2 = float(row["mH_input_GeV"])
            ma = float(row["mA_input_GeV"])
        except (KeyError, ValueError):
            continue
        stage = row.get("rejection_stage", "?")
        accepted = row.get("accepted") == "1"
        key = "accepted" if accepted else stage
        by_stage.setdefault(key, {"x": [], "y": []})
        by_stage[key]["x"].append(mh2)
        by_stage[key]["y"].append(ma)

    order = [k for k in STAGE_COLORS if k != "accepted" and k in by_stage] + (
        ["accepted"] if "accepted" in by_stage else []
    )
    for key in order:
        data = by_stage[key]
        color = STAGE_COLORS.get(key, "#1f77b4")
        size = 28 if key == "accepted" else 4
        alpha = 1.0 if key == "accepted" else 0.15
        zorder = 10 if key == "accepted" else 1
        ax.scatter(data["x"], data["y"], s=size, alpha=alpha, color=color,
                   label=f"{key} (n={len(data['x'])})", zorder=zorder,
                   edgecolors="black" if key == "accepted" else "none", linewidths=0.5)
    ax.set_xlabel("m_H2 [GeV]")
    ax.set_ylabel("m_A = m_Hp [GeV]")
    ax.set_title("PHYSICAL_POINT_SCAN: attempted points by rejection stage")
    ax.legend(loc="upper left", fontsize=8, markerscale=2)
    ax.grid(alpha=0.3)

    # Right panel: per mass_region acceptance rate bar chart
    ax2 = axes[1]
    region_counts = {}
    for row in rows:
        region = row.get("mass_region", "?")
        accepted = row.get("accepted") == "1"
        c = region_counts.setdefault(region, {"total": 0, "accepted": 0})
        c["total"] += 1
        c["accepted"] += int(accepted)
    regions_sorted = sorted(region_counts.keys(),
                            key=lambda r: region_counts[r].get("_mh2", 0) if False else r)
    # sort by embedded mass number if present (R150_..., R200_...)
    def sort_key(r):
        digits = "".join(ch for ch in r if ch.isdigit())
        return int(digits) if digits else 0
    regions_sorted = sorted(region_counts.keys(), key=sort_key)
    fracs = [100.0 * region_counts[r]["accepted"] / region_counts[r]["total"] for r in regions_sorted]
    counts_acc = [region_counts[r]["accepted"] for r in regions_sorted]
    bars = ax2.bar(regions_sorted, fracs, color="#2ca02c")
    for bar, n in zip(bars, counts_acc):
        ax2.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                 str(n), ha="center", va="bottom", fontsize=8)
    ax2.set_ylabel("% attempted points accepted (theory_ok & no forbidden cascade)")
    ax2.set_title("Acceptance rate by mass region (bar label = accepted count)")
    ax2.tick_params(axis="x", rotation=45)
    ax2.grid(alpha=0.3, axis="y")

    fig.tight_layout()
    png = RESULTS_DIR / "mass_validity_map.png"
    pdf = RESULTS_DIR / "mass_validity_map.pdf"
    fig.savefig(png, dpi=150)
    fig.savefig(pdf)
    print(f"Wrote {png}")
    print(f"Wrote {pdf}")


if __name__ == "__main__":
    main()
