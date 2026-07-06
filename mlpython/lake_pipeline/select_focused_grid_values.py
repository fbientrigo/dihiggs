#!/usr/bin/env python3
"""Select actual existing grid values for the focused-grid pipeline.

This campaign's grid is *diagonal*: each ``mA_target`` is bound one-to-one to a
single ``tan_beta`` and a single ``lambda6`` (see ``focused_grid_inventory``).
So the real unit of selection is the co-occurring
``(mA_target, tan_beta, lambda6)`` combo, driven by ``mA_target``.

Selection rules:

* ``--mA-targets`` picks the driving axis (default: 5 of the 6 available).
* ``tan_beta`` / ``lambda6`` are taken from co-occurrence for those targets.
* ``--tan-betas`` / ``--lambda6-values`` act only as an *optional narrowing
  filter*.  When a narrowing filter drops a selected ``mA_target``, that is
  reported loudly as a warning (and, if everything is dropped, the script fails
  clearly listing the available combos).

Only actual existing values are used; nothing is interpolated.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import polars as pl

DEFAULT_MA_TARGETS = [300.0, 350.0, 400.0, 450.0, 500.0]
# Relative tolerance for matching tan_beta (carries float noise, e.g.
# 10.00000000000001); mA_target and lambda6 are stored exactly.
TAN_BETA_RTOL = 1e-6
ABS_TOL = 1e-9


def parse_float_list(text: str | None) -> list[float] | None:
    if text is None:
        return None
    return [float(x) for x in text.replace(",", " ").split()]


def approx_in(value: float, candidates: list[float], rtol: float) -> bool:
    return any(
        abs(value - c) <= (rtol * abs(c) + ABS_TOL) for c in candidates
    )


def load_combos(input_path: Path) -> list[dict]:
    lf = pl.scan_parquet(input_path)
    cols = ["mA_target", "tan_beta", "lambda6", "lambda7"]
    df = lf.select(cols).unique().sort(cols).collect()
    return df.to_dicts()


def select(
    combos: list[dict],
    mA_targets: list[float],
    tan_betas: list[float] | None,
    lambda6_values: list[float] | None,
    lambda7: float,
) -> tuple[list[dict], list[str]]:
    warnings: list[str] = []

    available_mA = sorted({c["mA_target"] for c in combos})
    for v in mA_targets:
        if not approx_in(v, available_mA, ABS_TOL):
            raise SystemExit(
                f"[select] requested mA_target={v} not found. "
                f"Available: {available_mA}"
            )

    available_tb = sorted({c["tan_beta"] for c in combos})
    available_l6 = sorted({c["lambda6"] for c in combos})
    if tan_betas is not None:
        for v in tan_betas:
            if not approx_in(v, available_tb, TAN_BETA_RTOL):
                raise SystemExit(
                    f"[select] requested tan_beta={v} not found anywhere. "
                    f"Available: {available_tb}"
                )
    if lambda6_values is not None:
        for v in lambda6_values:
            if not approx_in(v, available_l6, ABS_TOL):
                raise SystemExit(
                    f"[select] requested lambda6={v} not found anywhere. "
                    f"Available: {available_l6}"
                )

    selected: list[dict] = []
    for v in mA_targets:
        combo = next(
            (c for c in combos
             if abs(c["mA_target"] - v) <= ABS_TOL
             and abs(c["lambda7"] - lambda7) <= ABS_TOL),
            None,
        )
        if combo is None:
            warnings.append(
                f"mA_target={v} has no combo with lambda7={lambda7}; skipped."
            )
            continue
        if tan_betas is not None and not approx_in(
            combo["tan_beta"], tan_betas, TAN_BETA_RTOL
        ):
            warnings.append(
                f"mA_target={v} dropped: its bound tan_beta={combo['tan_beta']} "
                f"is not in requested --tan-betas={tan_betas}."
            )
            continue
        if lambda6_values is not None and not approx_in(
            combo["lambda6"], lambda6_values, ABS_TOL
        ):
            warnings.append(
                f"mA_target={v} dropped: its bound lambda6={combo['lambda6']} "
                f"is not in requested --lambda6-values={lambda6_values}."
            )
            continue
        selected.append(combo)
    return selected, warnings


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="Path to silver_all.parquet")
    p.add_argument("--output-dir", required=True,
                   help="Directory for selected_values.json")
    p.add_argument("--mA-targets", default=None,
                   help="Comma list; default 300,350,400,450,500")
    p.add_argument("--tan-betas", default=None,
                   help="Optional narrowing filter (comma list)")
    p.add_argument("--lambda6-values", default=None,
                   help="Optional narrowing filter (comma list)")
    p.add_argument("--lambda7", type=float, default=0.0)
    p.add_argument("--lambda1-target", type=float, default=1.0)
    p.add_argument("--variation-idx", type=int, default=0)
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    input_path = Path(args.input)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    mA_targets = parse_float_list(args.mA_targets) or list(DEFAULT_MA_TARGETS)
    tan_betas = parse_float_list(args.tan_betas)
    lambda6_values = parse_float_list(args.lambda6_values)

    combos_all = load_combos(input_path)
    selected, warnings = select(
        combos_all, mA_targets, tan_betas, lambda6_values, args.lambda7
    )

    for w in warnings:
        print(f"[select][warn] {w}", file=sys.stderr)

    payload = {
        "input_parquet": str(input_path),
        "lambda7": args.lambda7,
        "lambda1_target": args.lambda1_target,
        "variation_idx": args.variation_idx,
        "requested": {
            "mA_targets": mA_targets,
            "tan_betas": tan_betas,
            "lambda6_values": lambda6_values,
        },
        "mA_targets": sorted({c["mA_target"] for c in selected}),
        "tan_beta_values": sorted({c["tan_beta"] for c in selected}),
        "lambda6_values": sorted({c["lambda6"] for c in selected}),
        "combos": sorted(
            selected, key=lambda c: (c["mA_target"], c["lambda6"], c["tan_beta"])
        ),
        "warnings": warnings,
    }

    out_path = output_dir / "selected_values.json"
    with out_path.open("w") as f:
        json.dump(payload, f, indent=2, sort_keys=True)
        f.write("\n")

    print(f"[select] selected {len(selected)} combo(s):")
    for c in payload["combos"]:
        print(
            f"  mA_target={c['mA_target']:g} tan_beta={c['tan_beta']:g} "
            f"lambda6={c['lambda6']:g} lambda7={c['lambda7']:g}"
        )
    print(f"[select] wrote {out_path}")

    if not selected:
        print(
            "[select] ERROR: no co-occurring combos match the requested "
            "narrowing filters.\n"
            f"[select] available combos: {combos_all}",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
