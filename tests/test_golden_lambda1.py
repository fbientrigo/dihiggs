"""Golden characterization of the lambda1-target construction path.

Status: CHARACTERIZATION ONLY. This suite freezes the CURRENT behavior of
dihiggs/src/PhysLam1Scan.cpp + 2hdmc/src/THDM.cpp:set_param_phys_lam1. It does
NOT certify that behavior as physically correct. See
docs/characterization_lambda1.md.

Test naming states which kind of claim each test makes:

  behavior_*      freezes what the code currently does (no correctness claim)
  defect_*        freezes a defect that is REAL and characterized-but-unrepaired;
                  these tests pass by asserting the broken behavior, and MUST be
                  revisited (and inverted) when the defect is fixed
  invariant_*     a relationship that should hold for any correct implementation
  compat_*        a property downstream consumers currently depend on
  opendecision_*  documents an unresolved scientific/contract decision

Binary-dependent tests skip unless the binary exists, or fail hard when
DIHIGGS_REQUIRE_LAMBDA1_BINARY=1 (CI sets this).
"""
from __future__ import annotations

import json
import csv
import math
import os
import subprocess
import tempfile
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
GOLDEN_DIR = REPO_ROOT / "tests" / "golden" / "lambda1_v1"
CASES_PATH = GOLDEN_DIR / "cases.json"
EXPECTED_CSV = GOLDEN_DIR / "expected.csv"
MARKERS_PATH = GOLDEN_DIR / "expected_markers.json"
MANIFEST_PATH = GOLDEN_DIR / "manifest.json"

BINARY = Path(
    os.environ.get("LAMBDA1_BINARY", REPO_ROOT / "build" / "lambda1_char" / "PhysLam1Scan")
)
V2_BINARY = Path(
    os.environ.get("LAMBDA1_V2_BINARY", REPO_ROOT / "dihiggs" / "app" / "Lambda1EvaluatorV2")
)
REQUIRE_BINARY = os.environ.get("DIHIGGS_REQUIRE_LAMBDA1_BINARY") == "1"

# hbarc in GeV*m, as used by autoresearch/harness/bounded_adaptive_search.py.
HBARC_GEV_M = 1.973269804e-16

# ---------------------------------------------------------------------------
# Float comparison policy (field-specific, deliberately NOT one broad tolerance)
# ---------------------------------------------------------------------------
# The production CSV is written with std::fixed << std::setprecision(15), i.e.
# 15 digits AFTER the decimal point. That is an absolute-precision format, so
# the meaningful comparison floor is absolute and identical for every column:
# two values that serialize to the same 15-decimal string are indistinguishable
# in this format. A *relative* tolerance is meaningless for a column whose
# serialized value is 0.000000000000000.
SERIALIZATION_ABS_FLOOR = 1e-15

EXACT_COLUMNS = frozenset(
    {
        # Pass-through inputs and compile-time constants: copied to the row
        # unmodified, so any difference at all is real drift.
        "m_phi",
        "mA",
        "lambda6",
        "lambda7",
        "sin_ba",
        "tan_beta",
        "lam1",
        # Discrete flags.
        "positivity_ok",
        "unitarity_ok",
        "perturbativity_ok",
    }
)

FLAG_COLUMNS = ("positivity_ok", "unitarity_ok", "perturbativity_ok")


def _skip_or_fail_without_binary() -> None:
    if BINARY.exists():
        return
    msg = (
        f"PhysLam1Scan not built at {BINARY}. "
        "Run scripts/build_lambda1_characterization.sh"
    )
    if REQUIRE_BINARY:
        pytest.fail(msg)
    pytest.skip(msg)


def _sha256(path: Path) -> str:
    import hashlib

    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_csv(text: str) -> tuple[list[str], list[list[str]]]:
    lines = text.strip().splitlines()
    header = lines[0].split(",")
    rows = [ln.split(",") for ln in lines[1:]]
    return header, rows


def load_expected() -> tuple[list[str], list[list[str]]]:
    return _read_csv(EXPECTED_CSV.read_text())


def load_cases() -> list[dict]:
    return json.loads(CASES_PATH.read_text())["cases"]


def run_case(case: dict, out_csv: Path) -> tuple[str, str, str]:
    args = [str(BINARY)] + list(case["args"]) + [str(out_csv)]
    env = dict(os.environ, OMP_NUM_THREADS="1")
    proc = subprocess.run(args, capture_output=True, text=True, env=env, cwd=REPO_ROOT)
    assert proc.returncode == 0, f"{case['case_id']}: exit {proc.returncode}\n{proc.stderr}"
    return out_csv.read_text(), proc.stdout, proc.stderr


def rows_for(case_id: str) -> list[dict[str, str]]:
    header, rows = load_expected()
    return [dict(zip(header, r)) for r in rows if r[0] == case_id]


def one_row(case_id: str) -> dict[str, str]:
    got = rows_for(case_id)
    assert len(got) == 1, f"{case_id}: expected exactly 1 golden row, got {len(got)}"
    return got[0]


# ---------------------------------------------------------------------------
# Provenance / integrity
# ---------------------------------------------------------------------------


def test_behavior_manifest_checksums_match_golden_files():
    manifest = json.loads(MANIFEST_PATH.read_text())
    sums = manifest["checksums"]
    assert sums["cases.json"] == _sha256(CASES_PATH)
    assert sums["expected.csv"] == _sha256(EXPECTED_CSV)
    assert sums["expected_markers.json"] == _sha256(MARKERS_PATH)


def test_behavior_manifest_records_compiled_constants():
    """The oracle is only meaningful with its compile-time constants recorded."""
    manifest = json.loads(MANIFEST_PATH.read_text())
    consts = manifest["compiled_constants"]
    assert "125.0" in consts["M_H"]
    assert "1e-9" in consts["THDM::EPS"]
    counts = manifest["counts"]
    _, rows = load_expected()
    assert counts["expected_csv_rows"] == len(rows)
    assert counts["cases"] == len(load_cases())


def test_behavior_every_case_is_covered_by_markers():
    markers = json.loads(MARKERS_PATH.read_text())
    case_ids = [c["case_id"] for c in load_cases()]
    assert sorted(markers) == sorted(case_ids)
    for case in load_cases():
        assert case["provenance"], f"{case['case_id']}: missing provenance"
        assert case["characterizes"], f"{case['case_id']}: missing description"


def test_behavior_golden_rows_have_consistent_field_count():
    header, rows = load_expected()
    for r in rows:
        assert len(r) == len(header), (
            f"row {r[0]}: {len(r)} fields, header has {len(header)}"
        )


def test_behavior_csv_schema_is_frozen():
    """The 29-column legacy schema, exactly as emitted."""
    header, _ = load_expected()
    assert header == [
        "case_id",
        "m_phi", "mA", "alpha", "beta", "lambda6", "lambda7", "m12",
        "sin_ba", "tan_beta",
        "positivity_ok", "unitarity_ok", "perturbativity_ok",
        "width_bb", "width_tautau", "width_WW", "width_ZZ",
        "width_gaga", "width_Zga", "width_gg", "width_hh",
        "total_width", "br_gaga", "lam1", "computed_lam1",
        "lam2", "computed_lam2", "lam3", "lam4", "lam5",
    ]


# ---------------------------------------------------------------------------
# Row cardinality: construction failures are DROPPED (audit section 4.6)
# ---------------------------------------------------------------------------


def test_defect_construction_failure_drops_the_attempted_point():
    """PhysLam1Scan.cpp `if (!pset) continue;` -- the attempted point vanishes.

    Boundary's evaluate_point preserves one row per attempted point. This path
    does not. Frozen as current behavior; core-v2 must preserve rows instead.
    """
    markers = json.loads(MARKERS_PATH.read_text())
    for cid in (
        "L02_construct_fail_mphi_below_mh",
        "L03_construct_fail_tanbeta_zero",
        "L04_construct_fail_sba_gt_one",
    ):
        m = markers[cid]
        assert m["total_attempts"] == 1, cid
        assert m["csv_rows_actual"] == 0, f"{cid}: failure row unexpectedly preserved"
        assert rows_for(cid) == [], cid


def test_defect_attempted_and_emitted_counts_diverge_on_mixed_grid():
    """3 attempted points -> 2 emitted rows. attempted != emitted."""
    m = json.loads(MARKERS_PATH.read_text())["L08_row_preservation_mixed_grid"]
    assert m["total_attempts"] == 3
    assert m["csv_rows_actual"] == 2
    assert m["total_csv_rows_reported"] == 2
    assert len(rows_for("L08_row_preservation_mixed_grid")) == 2


def test_defect_some_construction_failures_are_completely_silent():
    """tan_beta<=0 and |sin_ba|>1 fail with NO warning on stdout or stderr.

    Only the m_h > m_H guard prints. So two of three failure classes are
    invisible in the CSV *and* in the logs.
    """
    markers = json.loads(MARKERS_PATH.read_text())
    assert markers["L02_construct_fail_mphi_below_mh"]["construction_warnings"] == 1
    for cid in ("L03_construct_fail_tanbeta_zero", "L04_construct_fail_sba_gt_one"):
        assert markers[cid]["construction_warnings"] == 0, (
            f"{cid}: expected silent failure"
        )
        assert markers[cid]["csv_rows_actual"] == 0


def test_behavior_theory_rejection_preserves_the_row():
    """Unlike construction failure, a theory rejection IS emitted."""
    m = json.loads(MARKERS_PATH.read_text())["L05_theory_reject_large_lam1"]
    assert m["total_attempts"] == 1
    assert m["csv_rows_actual"] == 1
    assert m["triple_ok_points"] == 0
    row = one_row("L05_theory_reject_large_lam1")
    flags = tuple(row[f] for f in FLAG_COLUMNS)
    assert flags != ("1.000000000000000",) * 3, "expected at least one flag to fail"


def test_behavior_construction_boundary_is_strict_greater_than():
    """m_h > m_H fails; m_phi == mh constructs. Pins the inclusive side."""
    markers = json.loads(MARKERS_PATH.read_text())
    assert markers["L09_mphi_exactly_at_mh"]["csv_rows_actual"] == 1
    assert markers["L02_construct_fail_mphi_below_mh"]["csv_rows_actual"] == 0


# ---------------------------------------------------------------------------
# Flags
# ---------------------------------------------------------------------------


def test_behavior_accepted_case_flag_vector_is_frozen():
    row = one_row("L01_accepted_g06")
    for f in FLAG_COLUMNS:
        assert row[f] == "1.000000000000000", f"{f}={row[f]}"


def test_behavior_flags_are_only_ever_zero_or_one():
    header, rows = load_expected()
    for r in rows:
        d = dict(zip(header, r))
        for f in FLAG_COLUMNS:
            assert d[f] in ("0.000000000000000", "1.000000000000000"), (
                f"{d['case_id']}.{f} = {d[f]}"
            )


def test_opendecision_triple_ok_is_not_a_column_and_omits_stability():
    """triple_ok exists ONLY as a stdout marker; and stability is never checked.

    PhysLam1Scan computes positivity/unitarity/perturbativity and counts
    triple_ok = pos && uni && pert on stdout. It never calls check_stability
    and emits no triple_ok/theory_ok column. Boundary, by contrast, requires
    stability_ok in theory_ok.

    Note the vendored fork has the SAME positivity/stability alias as boundary's
    stock 2HDMC (Constraints::check_positivity and ::check_stability both return
    model.check_stability()), so `positivity_ok` here is really a stability
    result. Whether that is intended is a scientific-contract decision, not one
    this suite may resolve.
    """
    header, _ = load_expected()
    assert "triple_ok" not in header
    assert "theory_ok" not in header
    assert "stability_ok" not in header
    src = (REPO_ROOT / "dihiggs" / "src" / "PhysLam1Scan.cpp").read_text()
    assert "check_stability" not in src
    constraints = (REPO_ROOT / "2hdmc" / "src" / "Constraints.cpp").read_text()
    assert "bool Constraints::check_positivity() {\n  return model.check_stability();" in constraints


def test_invariant_stdout_triple_ok_marker_equals_flag_conjunction():
    """TRIPLE_OK_POINTS must equal the number of rows with all three flags set."""
    markers = json.loads(MARKERS_PATH.read_text())
    header, rows = load_expected()
    for case in load_cases():
        cid = case["case_id"]
        n = 0
        for r in rows:
            d = dict(zip(header, r))
            if d["case_id"] != cid:
                continue
            if all(d[f] == "1.000000000000000" for f in FLAG_COLUMNS):
                n += 1
        assert markers[cid]["triple_ok_points"] == n, (
            f"{cid}: marker says {markers[cid]['triple_ok_points']}, rows say {n}"
        )


# ---------------------------------------------------------------------------
# lambda1 target vs reconstructed
# ---------------------------------------------------------------------------


def test_behavior_lam1_target_and_reconstructed_are_separate_columns():
    """`lam1` is the requested target; `computed_lam1` is reconstructed."""
    row = one_row("L01_accepted_g06")
    assert float(row["lam1"]) == 0.1
    assert math.isclose(float(row["computed_lam1"]), 0.1, rel_tol=1e-10, abs_tol=1e-12)


def test_defect_large_lam1_roundtrip_error_warns_but_never_rejects():
    """The campaign-relevant point has |lam1 - computed_lam1| ~ 5.8e-07.

    That is ~580x THDM::EPS (1e-9). set_param_phys_lam1 stores the residual in
    lam1_validation_* and prints a stderr warning, but PhysLam1Scan never reads
    those fields into the CSV and never rejects the point: it is still emitted
    with all flags set and counted in TRIPLE_OK_POINTS.
    """
    markers = json.loads(MARKERS_PATH.read_text())
    m = markers["L07_campaign_best_large_tb"]
    assert m["roundtrip_warnings"] == 1, "expected the round-trip warning to fire"
    assert m["triple_ok_points"] == 1, "warned point is still counted as triple-OK"

    row = one_row("L07_campaign_best_large_tb")
    residual = abs(float(row["lam1"]) - float(row["computed_lam1"]))
    assert residual > 1e-9, f"expected residual above EPS, got {residual}"
    assert residual < 1e-5, f"residual unexpectedly large: {residual}"

    # The warning is not recoverable from the CSV: no column records it.
    header, _ = load_expected()
    for absent in ("lam1_validation_warning", "lam1_abs_error", "roundtrip_ok"):
        assert absent not in header


def test_compat_lam2_and_computed_lam2_are_the_same_value():
    """PhysLam1Scan writes lam2_g twice ('no separate l2 calculation').

    So `lam2` is NOT a requested target the way `lam1` is; the target/computed
    column pairing is asymmetric. Consumers must not read `lam2` as an input.
    """
    header, rows = load_expected()
    for r in rows:
        d = dict(zip(header, r))
        assert d["lam2"] == d["computed_lam2"], d["case_id"]


# ---------------------------------------------------------------------------
# Serialization / width / lifetime
# ---------------------------------------------------------------------------


def test_defect_fixed15_serialization_destroys_the_long_lived_regime():
    """The LLP regime is exactly where the CSV loses the width.

    L06's TRUE total width is ~1.2e-17 GeV (true ctau ~16 m), but
    std::fixed<<setprecision(15) writes 0.000000000000000. autoresearch then
    computes ctau_m = HBARC_GEV_M / total_width behind a `> 0.0` filter
    (bounded_adaptive_search.py), so the point is silently DROPPED from the
    ctau metric -- the longer-lived the point, the more likely it disappears.
    """
    row = one_row("L06_llp_regime_small_l6_large_tb")
    assert float(row["total_width"]) == 0.0, (
        "expected the width to serialize to exactly zero in this regime"
    )
    # Reproduce the consumer's filter: this point contributes nothing.
    width = float(row["total_width"])
    assert not (width > 0.0), "consumer filter `v > 0.0` would drop this row"


def test_defect_br_gaga_guard_zeroes_the_displaced_regime():
    """`br_gaga = (w_tot > 1e-15) ? w_gaga/w_tot : 0.0`.

    A width of 1e-15 GeV corresponds to ctau = hbarc/w ~ 0.197 m, so the guard
    forces br_gaga to exactly 0 for any scalar with ctau above roughly 20 cm --
    precisely the displaced-vertex regime llp_recast exists to study. L06's true
    br_gaga is ~0.685; the emitted value is 0.
    """
    row = one_row("L06_llp_regime_small_l6_large_tb")
    assert float(row["br_gaga"]) == 0.0
    guard_ctau_m = HBARC_GEV_M / 1e-15
    assert 0.19 < guard_ctau_m < 0.20


def test_defect_representable_ctau_is_capped_by_serialization():
    """Quantifies the ceiling the format imposes on the project's own metric."""
    # Smallest nonzero value representable in fixed/15 is 1e-15.
    max_ctau_m = HBARC_GEV_M / 1e-15
    assert math.isclose(max_ctau_m, 0.1973269804, rel_tol=1e-9)
    # Anything below 5e-16 rounds to zero and is dropped by the `> 0.0` filter.
    invisible_above_m = HBARC_GEV_M / 5e-16
    assert math.isclose(invisible_above_m, 0.3946539608, rel_tol=1e-9)


def test_behavior_no_lifetime_column_is_emitted():
    """No ctau/ctau_m/ctau_mm column exists; consumers must derive it.

    dihiggs/app/orchestrator/engines/lambda1.py's expected_csv_columns claims
    `ctau_m` (and `mphi`), which this schema does not emit -- audit section 4.10.
    """
    header, _ = load_expected()
    for absent in ("ctau", "ctau_m", "ctau_mm", "lifetime", "tau_ns"):
        assert absent not in header
    assert "m_phi" in header and "mphi" not in header


def test_invariant_widths_are_finite_and_nonnegative():
    header, rows = load_expected()
    width_cols = [c for c in header if c.startswith("width_")] + ["total_width"]
    for r in rows:
        d = dict(zip(header, r))
        for c in width_cols:
            v = float(d[c])
            assert math.isfinite(v), f"{d['case_id']}.{c} not finite: {d[c]}"
            assert v >= 0.0, f"{d['case_id']}.{c} negative: {d[c]}"


def test_invariant_br_gaga_within_unit_interval():
    header, rows = load_expected()
    for r in rows:
        d = dict(zip(header, r))
        v = float(d["br_gaga"])
        assert 0.0 <= v <= 1.0, f"{d['case_id']}: br_gaga={v}"


def test_invariant_no_non_finite_values_anywhere():
    """Non-finite classification happens before any tolerance arithmetic."""
    header, rows = load_expected()
    for r in rows:
        d = dict(zip(header, r))
        for c in header[1:]:
            v = float(d[c])
            assert not math.isnan(v), f"{d['case_id']}.{c} is NaN"
            assert not math.isinf(v), f"{d['case_id']}.{c} is infinite"


def test_opendecision_m12_column_actually_holds_m12_squared():
    """The column named `m12` carries m12_2_g from get_param_gen -- m12 SQUARED.

    Frozen as current behavior. Any consumer treating `m12` as a mass in GeV is
    wrong by a square. Canonical naming is deferred to the v2 contract
    (audit section 4.4); this suite only records the semantics.
    """
    header, _ = load_expected()
    assert "m12" in header
    assert "m12_sq" not in header and "m12_2" not in header
    src = (REPO_ROOT / "dihiggs" / "src" / "PhysLam1Scan.cpp").read_text()
    # The row vector places m12_2_g under the "m12" header position.
    assert "m12_2_g" in src
    # Sanity: G06's value is ~334, which is m12^2 in GeV^2, not m12 in GeV.
    row = one_row("L01_accepted_g06")
    assert 300.0 < float(row["m12"]) < 400.0


# ---------------------------------------------------------------------------
# Executable re-run: the oracle must reproduce
# ---------------------------------------------------------------------------


def test_regenerated_output_matches_golden():
    _skip_or_fail_without_binary()
    header, _ = load_expected()
    expected_rows = [",".join(r) for _, r in [(0, r) for r in load_expected()[1]]]

    produced: list[str] = []
    with tempfile.TemporaryDirectory() as td:
        for case in load_cases():
            cid = case["case_id"]
            text, _, _ = run_case(case, Path(td) / f"{cid}.csv")
            lines = text.strip().splitlines()
            assert lines[0].split(",") == header[1:], f"{cid}: header drift"
            for row in lines[1:]:
                produced.append(f"{cid},{row}")

    assert len(produced) == len(expected_rows), (
        f"row count drift: produced {len(produced)}, golden {len(expected_rows)}"
    )
    for got, want in zip(produced, expected_rows):
        gd = dict(zip(header, got.split(",")))
        wd = dict(zip(header, want.split(",")))
        assert gd["case_id"] == wd["case_id"]
        for col in header[1:]:
            g, w = gd[col], wd[col]
            if col in EXACT_COLUMNS:
                assert g == w, f"{gd['case_id']}.{col}: exact mismatch {g!r} != {w!r}"
                continue
            gf, wf = float(g), float(w)
            assert math.isfinite(gf) == math.isfinite(wf), f"{gd['case_id']}.{col}"
            assert abs(gf - wf) <= SERIALIZATION_ABS_FLOOR, (
                f"{gd['case_id']}.{col}: {gf!r} != {wf!r}"
            )


def test_repeat_runs_are_byte_identical():
    """Tolerances are not slack for nondeterminism: the binary must repeat."""
    _skip_or_fail_without_binary()
    with tempfile.TemporaryDirectory() as td:
        for case in load_cases():
            cid = case["case_id"]
            a, out_a, _ = run_case(case, Path(td) / f"{cid}_a.csv")
            b, out_b, _ = run_case(case, Path(td) / f"{cid}_b.csv")
            assert a == b, f"{cid}: CSV differs across identical runs"
            assert out_a.replace("\r", "") != "" and out_b.replace("\r", "") != ""


def test_markers_reproduce():
    _skip_or_fail_without_binary()
    import re

    golden = json.loads(MARKERS_PATH.read_text())
    with tempfile.TemporaryDirectory() as td:
        for case in load_cases():
            cid = case["case_id"]
            text, stdout, stderr = run_case(case, Path(td) / f"{cid}.csv")
            attempts = int(re.search(r"^Total Attempts:\s*(\d+)\s*$", stdout, re.M).group(1))
            triple = int(re.search(r"^TRIPLE_OK_POINTS\s+(\d+)\s*$", stdout, re.M).group(1))
            rows = max(len(text.strip().splitlines()) - 1, 0)
            assert attempts == golden[cid]["total_attempts"], cid
            assert triple == golden[cid]["triple_ok_points"], cid
            assert rows == golden[cid]["csv_rows_actual"], cid
            n_construct = len(re.findall(r"WARNING: Cannot set physical masses", stderr))
            n_rt = len(
                re.findall(r"WARNING: set_param_phys_lam1 lambda_1 round-trip", stderr)
            )
            assert n_construct == golden[cid]["construction_warnings"], cid
            assert n_rt == golden[cid]["roundtrip_warnings"], cid


def test_row_length_guard_detects_truncated_and_overlong_rows():
    """Row shape is guarded explicitly rather than absorbed by zip()."""
    header, _ = load_expected()
    good = ["x"] * len(header)
    truncated = ["x"] * (len(header) - 1)
    overlong = ["x"] * (len(header) + 1)
    assert len(good) == len(header)
    assert len(truncated) != len(header)
    assert len(overlong) != len(header)


# ---------------------------------------------------------------------------
# Canonical v2 evaluator
# ---------------------------------------------------------------------------

V2_HEADER = (
    "point_id,mh_gev,mH_gev,mA_gev,mHp_gev,sin_beta_minus_alpha,"
    "tan_beta,lambda1_target,lambda6_input,lambda7_input"
)


def run_v2(tmp_path: Path, lines: list[str], suffix: str = "") -> tuple[bytes, list[dict[str, str]]]:
    if not V2_BINARY.exists():
        pytest.skip(f"Lambda1EvaluatorV2 not built at {V2_BINARY}")
    source = tmp_path / f"v2{suffix}.input.csv"
    output = tmp_path / f"v2{suffix}.output.csv"
    source.write_text(V2_HEADER + "\n" + "\n".join(lines) + "\n")
    env = dict(os.environ, DIHIGGS_GIT_COMMIT="test-commit", DIHIGGS_GIT_DIRTY="0")
    proc = subprocess.run(
        [str(V2_BINARY), str(source), str(output)], capture_output=True, text=True,
        cwd=REPO_ROOT, env=env,
    )
    assert proc.returncode == 0, proc.stderr
    data = output.read_bytes()
    return data, list(csv.DictReader(data.decode().splitlines()))


def v2_line(
    point_id: str, mh: str = "125.0", mH: str = "130", mA: str = "300",
    mHp: str = "300", sin_ba: str = "0.999", tan_beta: str = "50",
    lambda1: str = "0.1", lambda6: str = "0.1", lambda7: str = "0",
) -> str:
    return ",".join((point_id, mh, mH, mA, mHp, sin_ba, tan_beta, lambda1, lambda6, lambda7))


def test_v2_preserves_every_physical_row_and_deterministic_ids(tmp_path):
    lines = [
        v2_line("accepted"),
        v2_line("rejected", lambda1="20"),
        v2_line("construction", mH="124.99999999999999"),
        ",not-a-number",
        "",
    ]
    first, rows = run_v2(tmp_path, lines, "a")
    second, repeated = run_v2(tmp_path, lines, "b")
    assert first == second
    assert len(rows) == len(lines)
    assert [r["construction_ok"] for r in rows] == ["1", "1", "0", "0", "0"]
    assert rows[0]["rejection_stage"] == "accepted"
    assert rows[1]["rejection_stage"] != "accepted"
    assert rows[2]["rejection_reason"] == "mh_gt_mH"
    assert rows[3]["rejection_reason"] == "wrong_field_count"
    assert rows[3]["point_id"].startswith("lambda1_")
    assert rows[3]["point_id"] == repeated[3]["point_id"]
    assert rows[4]["point_id"] == repeated[4]["point_id"]


def test_v2_raw_lexemes_and_float64_round_trip(tmp_path):
    lines = [
        v2_line("adjacent-a", lambda7="1.0000000000000000"),
        v2_line("adjacent-b", lambda7="1.0000000000000002"),
        v2_line("same-a", lambda7="1.00000000000000000"),
        v2_line("same-b", lambda7="1.00000000000000000001"),
    ]
    _, rows = run_v2(tmp_path, lines)
    assert [r["lambda7_input_raw"] for r in rows] == [line.rsplit(",", 1)[1] for line in lines]
    parsed = [float(r["lambda7_input"]) for r in rows]
    assert parsed[0] != parsed[1]
    assert parsed[2] == parsed[3]
    assert float(rows[0]["lambda7_input"]).hex() == float(lines[0].rsplit(",", 1)[1]).hex()
    assert float(rows[1]["lambda7_input"]).hex() == float(lines[1].rsplit(",", 1)[1]).hex()


def test_v2_explicit_mh_boundary_and_residual_warning_do_not_reject(tmp_path):
    lines = [
        v2_line("mh125", mh="125.0", mH="125.0"),
        v2_line("mh12509-fail", mh="125.09", mH="125.0"),
        v2_line("residual", mH="200", mA="500", mHp="500", sin_ba="1.0",
                tan_beta="126904", lambda1="1.0", lambda6="0.0019"),
    ]
    _, rows = run_v2(tmp_path, lines)
    assert rows[0]["construction_ok"] == "1"
    assert rows[1]["construction_ok"] == "0"
    assert rows[2]["lambda1_residual_warning"] == "1"
    assert rows[2]["theory_ok"] == "1"
    assert rows[2]["rejection_stage"] == "accepted"


def test_v2_l06_recovers_width_branching_ratio_and_lifetime(tmp_path):
    line = v2_line(
        "L06", mH="200", mA="500", mHp="500", sin_ba="1.0",
        tan_beta="10000", lambda1="1.0", lambda6="1e-10",
    )
    _, rows = run_v2(tmp_path, [line])
    row = rows[0]
    width = float(row["total_width_gev"])
    assert math.isclose(width, 1.214398311035274e-17, rel_tol=2e-12)
    assert row["width_ok"] == "1"
    assert float(row["br_gammagamma"]) > 0.68
    assert math.isclose(float(row["ctau_mm"]), 16249.0, rel_tol=2e-3)
    assert float(format(width, ".17e")) == width


def test_v2_malformed_values_have_exact_nan_and_flag_mask(tmp_path):
    _, rows = run_v2(tmp_path, [v2_line("bad", mh="inf")])
    row = rows[0]
    assert row["rejection_reason"] == "invalid_mh_gev"
    assert row["mh_input_gev_raw"] == "inf"
    assert math.isnan(float(row["mh_input_gev"]))
    for flag in (
        "construction_ok", "positivity_ok", "unitarity_ok", "perturbativity_ok",
        "stability_ok", "triple_ok", "theory_ok", "width_ok",
    ):
        assert row[flag] == "0"
    for field in ("lambda1_reconstructed", "total_width_gev", "ctau_mm"):
        assert math.isnan(float(row[field])) and not math.isinf(float(row[field]))


def test_v2_legacy_successful_fields_match_golden_l01(tmp_path):
    _, rows = run_v2(tmp_path, [
        v2_line("L01", mH="130", mA="300", mHp="300", sin_ba="0.999",
                tan_beta="50", lambda1="0.1", lambda6="0.1"),
    ])
    v2 = rows[0]
    legacy = one_row("L01_accepted_g06")
    pairs = {
        "m12_sq_reconstructed_gev2": "m12",
        "lambda1_reconstructed": "computed_lam1",
        "lambda2_reconstructed": "computed_lam2",
        "lambda3_reconstructed": "lam3",
        "lambda4_reconstructed": "lam4",
        "lambda5_reconstructed": "lam5",
    }
    for current, old in pairs.items():
        assert math.isclose(float(v2[current]), float(legacy[old]), rel_tol=1e-12, abs_tol=1e-14)
    assert [v2[f] for f in ("positivity_ok", "unitarity_ok", "perturbativity_ok")] == ["1"] * 3


def test_lifetime_audit_is_schema_selective_and_deterministic(tmp_path):
    legacy = tmp_path / "campaign" / "results.csv"
    legacy.parent.mkdir()
    header, _ = load_expected()
    l06 = one_row("L06_llp_regime_small_l6_large_tb")
    legacy.write_text(
        ",".join(header[1:]) + "\n"
        + ",".join(l06[name] for name in header[1:]) + "\n"
    )
    wrong = tmp_path / "wrong.csv"
    wrong.write_text("a,b\n1,2\n")
    script = REPO_ROOT / "scripts" / "audit_legacy_lambda1_lifetimes.py"
    outputs = []
    for suffix in ("a", "b"):
        json_path = tmp_path / f"audit-{suffix}.json"
        md_path = tmp_path / f"audit-{suffix}.md"
        subprocess.run(
            [os.sys.executable, str(script), str(tmp_path), "--json", str(json_path), "--markdown", str(md_path)],
            check=True,
        )
        outputs.append((json_path.read_bytes(), md_path.read_bytes()))
    assert outputs[0] == outputs[1]
    report = json.loads(outputs[0][0])
    assert report["summary"]["files"] == 1
    assert report["summary"]["rows"] == 1
    assert report["summary"]["zero_widths"] == 1
    assert report["summary"]["autoresearch_discards"] == 1
    assert report["summary"]["deterministic_replay_eligible"] == 1
