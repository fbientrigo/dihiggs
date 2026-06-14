"""
test_resume_stale_csv.py — Stale-CSV validation in resume.should_skip().

A non-empty output CSV may only short-circuit a task (skip / reuse) when the
CSV's recorded grid signature matches the current effective grid + fixed
parameters.  A leftover CSV from a different grid must NOT cause a silent skip.
"""

from __future__ import annotations

import json
from pathlib import Path

from dihiggs.app.orchestrator.resume import should_skip


def _make_csv_with_meta(tb_dir: Path, grid_sig: str | None, *, event: str = "done") -> Path:
    tb_dir.mkdir(parents=True, exist_ok=True)
    csv_path = tb_dir / "scan_tb_10000.csv"
    csv_path.write_text("mphi,lam1,ctau_m\n200.0,1.0,5.39e-4\n", encoding="utf-8")
    if grid_sig is not None:
        meta = {"event": event, "grid_signature": grid_sig}
        (tb_dir / "scan_meta.json").write_text(json.dumps(meta), encoding="utf-8")
    return csv_path


def test_matching_grid_signature_skips(tmp_path: Path) -> None:
    csv_path = _make_csv_with_meta(tmp_path / "tb_10000", "SIG_A")
    skip, reason = should_skip(
        output_csv=csv_path,
        grid_sig="SIG_A",
        force=False,
        done_signatures=set(),
    )
    assert skip is True
    assert reason == "output_csv_exists"


def test_changed_grid_signature_does_not_skip_stale_csv(tmp_path: Path) -> None:
    # Same CSV path, but the current grid signature differs from the one
    # recorded in scan_meta.json: the CSV is stale and must be re-run.
    csv_path = _make_csv_with_meta(tmp_path / "tb_10000", "SIG_A")
    skip, reason = should_skip(
        output_csv=csv_path,
        grid_sig="SIG_B",
        force=False,
        done_signatures=set(),
    )
    assert skip is False
    assert reason == "stale_csv_grid_signature_mismatch"


def test_missing_meta_does_not_skip(tmp_path: Path) -> None:
    csv_path = _make_csv_with_meta(tmp_path / "tb_10000", None)
    skip, reason = should_skip(
        output_csv=csv_path,
        grid_sig="SIG_A",
        force=False,
        done_signatures=set(),
    )
    assert skip is False
    assert reason == "stale_csv_meta_missing_or_not_done"


def test_meta_not_done_does_not_skip(tmp_path: Path) -> None:
    # Partial/failed run: meta exists with matching signature but event != done.
    csv_path = _make_csv_with_meta(tmp_path / "tb_10000", "SIG_A", event="fail")
    skip, reason = should_skip(
        output_csv=csv_path,
        grid_sig="SIG_A",
        force=False,
        done_signatures=set(),
    )
    assert skip is False
    assert reason == "stale_csv_meta_missing_or_not_done"


def test_force_never_skips(tmp_path: Path) -> None:
    csv_path = _make_csv_with_meta(tmp_path / "tb_10000", "SIG_A")
    skip, _ = should_skip(
        output_csv=csv_path,
        grid_sig="SIG_A",
        force=True,
        done_signatures=set(),
    )
    assert skip is False
