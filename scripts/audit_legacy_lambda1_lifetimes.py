#!/usr/bin/env python3
"""Read-only audit of legacy fixed/15 lambda1 CSV lifetime integrity."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

HBARC_GEV_MM = 1.973269804e-13
HEADER = (
    "m_phi,mA,alpha,beta,lambda6,lambda7,m12,sin_ba,tan_beta,"
    "positivity_ok,unitarity_ok,perturbativity_ok,width_bb,width_tautau,"
    "width_WW,width_ZZ,width_gaga,width_Zga,width_gg,width_hh,total_width,"
    "br_gaga,lam1,computed_lam1,lam2,computed_lam2,lam3,lam4,lam5"
).split(",")
REPLAY_FIELDS = ("m_phi", "mA", "sin_ba", "tan_beta", "lam1", "lambda6", "lambda7")
EXCLUDED = ("recomputed", "backup", "recovered", "artifact", "v2")


def classification(path: Path) -> str:
    text = str(path).lower()
    if "supported" in text:
        return "autoresearch-supported"
    if "debug" in text or "smoke" in text:
        return "diagnostic-discard"
    return "legacy-unclassified"


def eligible(row: dict[str, str]) -> bool:
    try:
        return all(math.isfinite(float(row[name])) for name in REPLAY_FIELDS)
    except (KeyError, ValueError):
        return False


def candidates(roots: list[Path]) -> list[Path]:
    found: set[Path] = set()
    for root in roots:
        paths = [root] if root.is_file() else root.rglob("*.csv")
        for path in paths:
            if any(term in part.lower() for part in path.parts for term in EXCLUDED):
                continue
            try:
                with path.open(newline="") as stream:
                    if next(csv.reader(stream), []) == HEADER:
                        found.add(path.resolve())
            except (OSError, UnicodeDecodeError):
                continue
    return sorted(found, key=str)


def audit(roots: list[Path]) -> dict:
    paths = candidates(roots)
    totals = Counter()
    positive_widths: list[float] = []
    classes = Counter()
    checksums = {"all": hashlib.sha256()}
    for path in paths:
        kind = classification(path)
        with path.open(newline="") as stream:
            for row in csv.DictReader(stream):
                totals["rows"] += 1
                classes[kind] += 1
                replayable = eligible(row)
                totals["recoverable"] += replayable
                try:
                    width = float(row["total_width"])
                except ValueError:
                    totals["invalid_widths"] += 1
                    continue
                if width == 0.0:
                    totals["zero_widths"] += 1
                    totals["autoresearch_discards"] += 1
                elif math.isfinite(width) and width > 0.0:
                    positive_widths.append(width)
        checksums.setdefault(kind, hashlib.sha256())
        digest_record = f"{path.resolve()}\0{hashlib.sha256(path.read_bytes()).hexdigest()}\n".encode()
        checksums["all"].update(digest_record)
        checksums[kind].update(digest_record)
    minimum = min(positive_widths) if positive_widths else None
    return {
        "audit_schema_version": "dihiggs.lambda1-lifetime-audit.v1",
        "selection": {
            "exact_header": HEADER,
            "excluded_path_terms": list(EXCLUDED),
            "historical_mh_gev": 125.0,
        },
        "roots": [str(root.resolve()) for root in roots],
        "summary": {
            "files": len(paths),
            "rows": totals["rows"],
            "zero_widths": totals["zero_widths"],
            "minimum_positive_width_gev": minimum,
            "apparent_maximum_ctau_mm": HBARC_GEV_MM / minimum if minimum else None,
            "recoverable": totals["recoverable"],
            "deterministic_replay_eligible": totals["recoverable"],
            "autoresearch_discards": totals["autoresearch_discards"],
            "ranking_boundary_rows": positive_widths.count(minimum) if minimum else 0,
            "campaign_classification_rows": dict(sorted(classes.items())),
        },
        "checksums": {
            "algorithm": "sha256(path NUL file_sha256 LF), sorted by path",
            "all_selected_files": checksums["all"].hexdigest(),
            "by_campaign_classification": {
                name: digest.hexdigest() for name, digest in sorted(checksums.items()) if name != "all"
            },
        },
    }


def markdown(report: dict) -> str:
    summary = report["summary"]
    minimum = summary["minimum_positive_width_gev"]
    lifetime = summary["apparent_maximum_ctau_mm"]
    lines = [
        "# Historical lambda1 lifetime serialization audit",
        "",
        "This read-only audit selects only the exact legacy fixed/15 schema and excludes",
        "recomputed, backup, recovered, artifact, and v2 paths. The historical evaluator",
        "used explicit-in-code `mh = 125.0 GeV`; this is provenance, not a v2 default.",
        f"Audited roots: {', '.join(report['roots'])}.",
        "",
        "| Metric | Value |",
        "|---|---:|",
        f"| Files | {summary['files']} |",
        f"| Rows | {summary['rows']} |",
        f"| Serialized zero widths | {summary['zero_widths']} |",
        f"| Minimum positive width (GeV) | {minimum:.17e} |" if minimum else "| Minimum positive width (GeV) | n/a |",
        f"| Apparent maximum ctau (mm) | {lifetime:.17e} |" if lifetime else "| Apparent maximum ctau (mm) | n/a |",
        f"| Recoverable coordinate rows | {summary['recoverable']} |",
        f"| Deterministic replay eligible | {summary['deterministic_replay_eligible']} |",
        f"| Autoresearch `width > 0` discards | {summary['autoresearch_discards']} |",
        f"| Positive-width ranking boundary rows | {summary['ranking_boundary_rows']} |",
        "",
        "## Campaign classification",
        "",
        "| Classification | Rows |",
        "|---|---:|",
    ]
    lines.extend(
        f"| {name} | {count} |"
        for name, count in summary["campaign_classification_rows"].items()
    )
    lines += [
        "",
        "## Interpretation",
        "",
        "A serialized zero cannot recover a width from the CSV alone. Replay is eligible only",
        "when all legacy lambda1 coordinates are finite and the historical `mh` provenance is",
        "known. Re-run only zero-width or positive-width ranking-boundary rows; never a complete grid.",
        "Aggregate source checksums and campaign classifications are in the JSON manifest.",
        "",
    ]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("roots", nargs="+", type=Path)
    parser.add_argument("--json", required=True, type=Path)
    parser.add_argument("--markdown", required=True, type=Path)
    args = parser.parse_args()
    report = audit(args.roots)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.markdown.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    args.markdown.write_text(markdown(report))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
