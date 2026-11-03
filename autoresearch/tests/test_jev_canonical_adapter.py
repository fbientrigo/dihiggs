from __future__ import annotations

import csv
import json
from pathlib import Path
import pytest

from dihiggs.app.orchestrator.lambda1_v2 import (
    INPUT_HEADER,
    Lambda1Fixed,
    SCHEMA_VERSION as CANONICAL_SCHEMA,
    cartesian_rows,
)
from scripts.jev_adapter import (
    CANONICAL_ENGINE,
    CANONICAL_PRODUCER,
    CANONICAL_SCHEMA,
    CanonicalPlan,
    ctau_m_to_mm,
    ctau_mm_to_m,
    execute_canonical_plan,
    extract_campaign_outcome,
    extract_subrun_outcome,
    materialize_canonical_plan,
    plan_template_v2,
    resolve_template_grid,
    update_scientific_state,
    validate_canonical_executable,
)
from scripts.jev_policy import build_scientific_state
from scripts.jev_decide_once import baseline_state
from dihiggs.app.orchestrator.physics_authority import active_mass_gev, convention_provenance


def _sample_cfg() -> dict:
    return {
        "baseline_best_ctau_mm": 0.5389973717365447,
        "target_ctau_mm": 1.0779947434730894,
        "baseline_best_ctau_m": 0.0005389973717365447,
        "target_ctau_m": 0.0010779947434730894,
        "max_iterations": 3,
        "max_walltime_minutes": 60,
        "max_points_total": 9000,
        "min_points_per_subrun": 1000,
        "max_points_per_subrun": 3000,
        "allowed_actions": ["run_template", "summarize_only", "stop"],
        "forbidden_actions": [
            "edit_cpp",
            "edit_scoring",
            "edit_triple_ok",
            "delete_artifacts",
            "widen_envelope",
            "broad_scan",
        ],
        "envelope": {
            "sin_ba": [1.0, 1.0],
            "lambda7": [0.0, 0.0],
            "lambda1": [1.0, 1.0],
            "lambda6": [0.0003, 0.0060],
            "tan_beta": [30000, 500000],
            "mphi": [160, 350],
            "mA": [300, 500],
        },
    }


# ---------------------------------------------------------------------------
# 1. Canonical template -> subrun mapping & 2. Exact point preservation
# ---------------------------------------------------------------------------


def test_refine_template_exact_point_preservation_and_decomposition(tmp_path: Path) -> None:
    cfg = _sample_cfg()
    plan = plan_template_v2(
        template_name="ctau_mphi_mA_refine",
        cfg=cfg,
        iteration=1,
        cycle_id="cycle_test",
        iteration_id="iter01_test",
        output_root=tmp_path / "jev_harness",
    )

    # 3 mA × 3 lambda6 = 9 subruns
    assert plan.subrun_count == 9
    assert len(plan.subruns) == 9

    # Each subrun: 3 tan_beta × 40 mH × 1 lambda1 = 120 points
    for sub in plan.subruns:
        assert sub.row_count == 120
        assert len(sub.tan_beta_values) == 3
        assert len(sub.mH_values) == 40
        assert sub.lambda1_values == [1.0]

    # Exactly 1080 total points
    assert plan.total_row_count == 1080
    assert plan.total_requested_points == 1080

    # Verify rows generated per subrun sum to 1080 unique rows
    all_point_ids = set()
    for sub in plan.subruns:
        rows = cartesian_rows(
            fixed=sub.fixed,
            mH_values=sub.mH_values,
            lambda1_values=sub.lambda1_values,
            tan_beta_values=sub.tan_beta_values,
        )
        assert len(rows) == 120
        for r in rows:
            assert list(r.keys()) == list(INPUT_HEADER)
            all_point_ids.add(r["point_id"])

    assert len(all_point_ids) == 1080


def test_template_subrun_mapping_for_all_allowlisted_templates(tmp_path: Path) -> None:
    cfg = _sample_cfg()

    # mphi_lower_edge_probe: 1 mA × 3 lambda6 = 3 subruns, 5 tb × 80 mH = 400 pts -> 1200 total
    plan_low = plan_template_v2("mphi_lower_edge_probe", cfg, output_root=tmp_path)
    assert plan_low.subrun_count == 3
    assert plan_low.total_row_count == 1200
    for s in plan_low.subruns:
        assert s.row_count == 400

    # mA_upper_edge_probe: 3 mA × 1 lambda6 = 3 subruns, 3 tb × 120 mH = 360 pts -> 1080 total
    plan_ma = plan_template_v2("mA_upper_edge_probe", cfg, output_root=tmp_path)
    assert plan_ma.subrun_count == 3
    assert plan_ma.total_row_count == 1080
    for s in plan_ma.subruns:
        assert s.row_count == 360

    # control_box: 1 mA × 1 lambda6 = 1 subrun, 1 tb × 1000 mH = 1000 pts -> 1000 total
    plan_box = plan_template_v2("control_box", cfg, output_root=tmp_path)
    assert plan_box.subrun_count == 1
    assert plan_box.total_row_count == 1000

    # sensitivity_probe (default variable=tan_beta): 1 mA × 1 lambda6 = 1 subrun
    plan_sens = plan_template_v2("sensitivity_probe", cfg, output_root=tmp_path)
    assert plan_sens.subrun_count == 1
    assert plan_sens.total_row_count >= 1000


# ---------------------------------------------------------------------------
# 3. Correct fixed-group decomposition
# ---------------------------------------------------------------------------


def test_subrun_fixed_parameters_are_constant_per_subrun(tmp_path: Path) -> None:
    cfg = _sample_cfg()
    plan = plan_template_v2("ctau_mphi_mA_refine", cfg, output_root=tmp_path)

    seen_fixed_pairs = set()
    for sub in plan.subruns:
        pair = (sub.fixed.mA, sub.fixed.lambda6)
        assert pair not in seen_fixed_pairs, "Each subrun must represent a unique (mA, lambda6) pair"
        seen_fixed_pairs.add(pair)
        assert sub.fixed.sin_ba == 1.0
        assert sub.fixed.lambda7 == 0.0
        assert sub.fixed.mHp == sub.fixed.mA
        assert sub.fixed.mh == active_mass_gev()


# ---------------------------------------------------------------------------
# 4. Canonical schema and producer identification
# ---------------------------------------------------------------------------


def test_canonical_schema_and_producer_identity(tmp_path: Path) -> None:
    cfg = _sample_cfg()
    plan = plan_template_v2("ctau_mphi_mA_refine", cfg, output_root=tmp_path)
    assert plan.producer == CANONICAL_PRODUCER
    assert plan.producer == "Lambda1EvaluatorV2"
    assert plan.producer_schema == CANONICAL_SCHEMA
    assert plan.producer_schema == "dihiggs.lambda1.v2"
    assert plan.engine == CANONICAL_ENGINE
    assert plan.engine == "lambda1_v2"


# ---------------------------------------------------------------------------
# 5. Explicit ctau unit handling and factor-1000 regression prevention
# ---------------------------------------------------------------------------


def test_explicit_ctau_unit_conversion_and_regression_prevention() -> None:
    ctau_mm = 0.5389973717365447
    ctau_m = 0.0005389973717365447

    assert abs(ctau_m_to_mm(ctau_m) - ctau_mm) < 1e-12
    assert abs(ctau_mm_to_m(ctau_mm) - ctau_m) < 1e-12

    # Factor-1000 check: if a value in mm is incorrectly used as m, it's 1000x larger
    assert ctau_m_to_mm(1.0) == 1000.0
    assert ctau_mm_to_m(1000.0) == 1.0

    # Verify build_scientific_state unit encoding
    cfg = _sample_cfg()
    worker_state = {"best_ctau_mm": ctau_mm, "executed_points_total": 100}
    state = build_scientific_state(cfg, worker_state)
    assert state["objective"]["unit"] == "mm"
    assert state["objective"]["best_ctau_mm"] == ctau_mm
    assert abs(state["objective"]["best_ctau_m"] - ctau_m) < 1e-12

    # Backward compatibility: state with only best_ctau_m correctly derives best_ctau_mm
    legacy_worker_state = {"best_ctau_m": ctau_m, "executed_points_total": 100}
    legacy_state = build_scientific_state(cfg, legacy_worker_state)
    assert legacy_state["objective"]["unit"] == "mm"
    assert abs(legacy_state["objective"]["best_ctau_mm"] - ctau_mm) < 1e-12
    assert legacy_state["objective"]["best_ctau_m"] == ctau_m


# ---------------------------------------------------------------------------
# 6. Dynamic paths: no hardcoded stale paths
# ---------------------------------------------------------------------------


def test_dynamic_paths_contain_no_stale_roots(tmp_path: Path) -> None:
    cfg = _sample_cfg()
    custom_root = tmp_path / "custom_runs" / "jev_harness"
    plan = plan_template_v2("control_box", cfg, output_root=custom_root)
    assert plan.output_root == custom_root

    summary = materialize_canonical_plan(plan, custom_root / "plan_dir")
    plan_json = json.loads(Path(summary["plan_path"]).read_text())

    forbidden_patterns = ["/home/fabi", "/home/fabian/dihiggs", "scripts/out", "PhysScanWithFixings", "lambda1_legacy"]
    plan_str = json.dumps(plan_json)
    for pattern in forbidden_patterns:
        assert pattern not in plan_str, f"Stale pattern {pattern!r} found in canonical plan"


# ---------------------------------------------------------------------------
# 7. Plan-only behavior (no C++ execution)
# ---------------------------------------------------------------------------


def test_plan_only_materializes_without_evaluator_binary(tmp_path: Path) -> None:
    cfg = _sample_cfg()
    plan_dir = tmp_path / "plancheck_cycle01_iter01"
    plan = plan_template_v2("ctau_mphi_mA_refine", cfg, output_root=tmp_path)

    # Use a nonexistent executable path; dry_run must succeed without executing C++
    dummy_exec = tmp_path / "nonexistent_Lambda1EvaluatorV2"
    summary = materialize_canonical_plan(plan, plan_dir, dry_run=True, executable=dummy_exec)

    assert summary["plancheck_ok"] is True
    assert summary["planned_points"] == 1080
    assert summary["subrun_count"] == 9
    assert Path(summary["plan_path"]).exists()

    # Each subrun directory should have input_v2.csv and run_manifest.json
    for sub in summary["subruns"]:
        input_csv = Path(sub["input_csv"])
        assert input_csv.exists()
        lines = input_csv.read_text().splitlines()
        assert len(lines) == 121  # 1 header + 120 rows
        assert lines[0] == ",".join(INPUT_HEADER)


# ---------------------------------------------------------------------------
# 8. Invalid/noncanonical executor rejection
# ---------------------------------------------------------------------------


def test_validate_canonical_executable_rejects_legacy(tmp_path: Path) -> None:
    legacy_binary = tmp_path / "PhysScanWithFixings"
    legacy_binary.write_text("#!/bin/sh\nexit 0\n")
    legacy_binary.chmod(0o755)

    with pytest.raises(ValueError, match="Non-canonical executor rejected"):
        validate_canonical_executable(legacy_binary)

    nonexistent = tmp_path / "Lambda1EvaluatorV2"
    with pytest.raises(FileNotFoundError):
        validate_canonical_executable(nonexistent)


# ---------------------------------------------------------------------------
# 9. Allowlisted actions & stop condition
# ---------------------------------------------------------------------------


def test_resolve_template_grid_rejects_unallowlisted_or_forbidden() -> None:
    cfg = _sample_cfg()
    with pytest.raises(ValueError, match="unknown template"):
        resolve_template_grid("unallowlisted_scan", cfg["envelope"])

    with pytest.raises(ValueError, match="forbidden override"):
        resolve_template_grid("control_box", cfg["envelope"], overrides={"command": "rm -rf"})

    with pytest.raises(ValueError, match="override not allowed"):
        resolve_template_grid("control_box", cfg["envelope"], overrides={"unknown_key": 123})


# ---------------------------------------------------------------------------
# 10. Outcome extraction with synthetic canonical CSV
# ---------------------------------------------------------------------------


def test_outcome_extraction_with_synthetic_canonical_csv(tmp_path: Path) -> None:
    csv_path = tmp_path / "output_v2.csv"
    header = (
        "schema_version,evaluator_commit,evaluator_dirty,evaluator_api,point_id,"
        "mh_input_gev_raw,mh_input_gev,mH_input_gev_raw,mH_input_gev,"
        "mA_input_gev_raw,mA_input_gev,mHp_input_gev_raw,mHp_input_gev,"
        "sin_beta_minus_alpha_input_raw,sin_beta_minus_alpha_input,"
        "tan_beta_input_raw,tan_beta_input,lambda1_target_raw,lambda1_target,"
        "lambda6_input_raw,lambda6_input,lambda7_input_raw,lambda7_input,"
        "construction_ok,yukawa_type_installed,rejection_stage,rejection_reason,"
        "lambda1_reconstructed,lambda2_reconstructed,lambda3_reconstructed,"
        "lambda4_reconstructed,lambda5_reconstructed,lambda6_reconstructed,lambda7_reconstructed,"
        "lambda1_abs_residual,lambda1_residual_warning,m12_sq_reconstructed_gev2,tan_beta_reconstructed,"
        "positivity_ok,unitarity_ok,perturbativity_ok,stability_ok,triple_ok,theory_ok,"
        "width_bb_gev,width_cc_gev,width_tautau_gev,width_WW_gev,width_ZZ_gev,"
        "width_gammagamma_gev,width_Zgamma_gev,width_gg_gev,width_hh_gev,total_width_gev,width_unaccounted_gev,"
        "br_bb,br_cc,br_tautau,br_WW,br_ZZ,br_gammagamma,br_Zgamma,br_gg,br_hh,"
        "width_ok,ctau_mm"
    )

    # 3 rows:
    # Row 1: valid, triple_ok=1, ctau_mm=0.55
    # Row 2: valid, triple_ok=1, ctau_mm=0.88 (best)
    # Row 3: theory rejected, positivity=0, triple_ok=0, ctau_mm=nan
    rows = [
        "dihiggs.lambda1.v2,c1,no,api,p1,125.2,125.2,200,200,500,500,500,500,1,1,124416,124416,1,1,0.0019,0.0019,0,0,1,1,accepted,none,1,0.2,7,-3,-3,0.0019,0,1e-7,0,100,124416,1,1,1,1,1,1,1e-13,1e-14,1e-14,0,0,1e-15,1e-15,1e-14,0,3e-13,1e-16,0.6,0.03,0.07,0,0,0.01,0.005,0.2,0,1,0.55",
        "dihiggs.lambda1.v2,c1,no,api,p2,125.2,125.2,210,210,500,500,500,500,1,1,124416,124416,1,1,0.0019,0.0019,0,0,1,1,accepted,none,1,0.2,7,-3,-3,0.0019,0,1e-7,0,100,124416,1,1,1,1,1,1,1e-13,1e-14,1e-14,0,0,1e-15,1e-15,1e-14,0,2e-13,1e-16,0.6,0.03,0.07,0,0,0.01,0.005,0.2,0,1,0.88",
        "dihiggs.lambda1.v2,c1,no,api,p3,125.2,125.2,220,220,500,500,500,500,1,1,124416,124416,1,1,0.0019,0.0019,0,0,1,1,positivity,check_positivity_false,1,0.2,7,-3,-3,0.0019,0,1e-7,0,100,124416,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,nan",
    ]
    csv_path.write_text(header + "\n" + "\n".join(rows) + "\n", encoding="utf-8")

    outcome = extract_subrun_outcome(csv_path)
    assert outcome["total_evaluated_rows"] == 3
    assert outcome["construction_ok_rows"] == 3
    assert outcome["triple_ok_rows"] == 2
    assert outcome["valid_rows"] == 2
    assert outcome["best_ctau_mm"] == 0.88
    assert abs(outcome["best_ctau_m"] - 0.00088) < 1e-12
    assert outcome["best_point"]["point_id"] == "p2"
    assert outcome["best_point"]["mH_gev"] == 210.0
    assert outcome["failures_by_stage"]["positivity"] == 1
    assert outcome["failures_by_reason"]["check_positivity_false"] == 1

    # Campaign aggregation
    camp_outcome = extract_campaign_outcome([csv_path])
    assert camp_outcome["total_evaluated_rows"] == 3
    assert camp_outcome["valid_rows"] == 2
    assert camp_outcome["best_ctau_mm"] == 0.88
    assert camp_outcome["best_point"]["point_id"] == "p2"


# ---------------------------------------------------------------------------
# 11. Scientific state update from outcome
# ---------------------------------------------------------------------------


def test_scientific_state_update_on_improvement() -> None:
    cfg = _sample_cfg()
    initial_state = {
        "status": "jev_baseline",
        "best_ctau_mm": 0.50,
        "best_ctau_m": 0.0005,
        "executed_points_total": 0,
        "iterations_completed": 0,
        "no_improvement_count": 0,
        "rejected_points_total": 0,
        "history": [],
    }
    outcome = {
        "total_evaluated_rows": 1080,
        "valid_rows": 1000,
        "triple_ok_fraction": 0.95,
        "best_ctau_mm": 0.75,
        "best_ctau_m": 0.00075,
        "best_point": {"point_id": "best_01", "mH_gev": 200.0},
    }

    updated = update_scientific_state(
        current_state=initial_state,
        config=cfg,
        outcome=outcome,
        iteration=1,
        cycle_id="cycle01",
        template_name="ctau_mphi_mA_refine",
    )

    assert updated["executed_points_total"] == 1080
    assert updated["rejected_points_total"] == 80
    assert updated["iterations_completed"] == 1
    assert updated["best_ctau_mm"] == 0.75
    assert updated["best_ctau_m"] == 0.00075
    assert updated["no_improvement_count"] == 0
    assert updated["last_iteration_summary_for_prompt"]["improved"] is True
    assert len(updated["history"]) == 1


# ---------------------------------------------------------------------------
# 10. Physics authority enforcement and provenance
# ---------------------------------------------------------------------------


def test_plan_embeds_active_physics_authority_provenance(tmp_path: Path) -> None:
    cfg = _sample_cfg()
    active_mh = active_mass_gev()
    expected = convention_provenance(source_commit=None, requested_mass_gev=active_mh)
    plan = plan_template_v2("ctau_mphi_mA_refine", cfg, output_root=tmp_path)
    assert plan.fixed_controls["mh_gev"] == active_mh
    assert plan.fixed_controls["mass_convention_id"] == expected["convention_id"]
    assert plan.mass_convention["convention_id"] == expected["convention_id"]
    assert plan.mass_convention["m_h_GeV"] == active_mh
    assert plan.mass_convention["m_h_GeV_text"] == expected["m_h_GeV_text"]
    assert plan.mass_convention["schema_version"] == expected["schema_version"]

    report = materialize_canonical_plan(plan, tmp_path / "plan_dir", dry_run=True)
    assert report["mass_convention"]["convention_id"] == expected["convention_id"]
    assert report["mass_convention"]["m_h_GeV"] == active_mh

    plan_json = json.loads((tmp_path / "plan_dir" / "canonical_plan.json").read_text())
    assert plan_json["mass_convention"]["convention_id"] == expected["convention_id"]
    assert plan_json["mass_convention"]["m_h_GeV"] == active_mh


def test_mismatched_explicit_mass_override_fails_closed(tmp_path: Path) -> None:
    cfg = _sample_cfg()
    active_mh = active_mass_gev()
    # Any explicit mass that conflicts with the active authority must fail closed.
    # Pick a value guaranteed to differ from whatever is currently active.
    conflicting_mh = 125.13 if active_mh != 125.13 else 125.20

    with pytest.raises(ValueError, match="Mismatched explicit Higgs mass"):
        plan_template_v2(
            "ctau_mphi_mA_refine",
            cfg,
            overrides={"mh": conflicting_mh},
            output_root=tmp_path,
        )

    # In cfg directly
    cfg_with_bad_mass = dict(cfg, mh_gev=conflicting_mh)
    with pytest.raises(ValueError, match="Mismatched explicit Higgs mass"):
        plan_template_v2(
            "ctau_mphi_mA_refine",
            cfg_with_bad_mass,
            output_root=tmp_path,
        )


# ---------------------------------------------------------------------------
# 11. Sample worker config: baseline provenance and factor-1000 guard
# ---------------------------------------------------------------------------


def test_sample_worker_config_baseline_is_canonical_replay_not_historical() -> None:
    """Loads scripts/jev_worker_config.sample.json and checks the baseline/target
    it carries: target = 2x baseline, the baseline matches its own recorded replay
    provenance (not the unrecoverable historical value), and best_ctau_mm/m are
    always derived from mm (never a stale stored m), which catches a factor-1000
    unit error.
    """
    import json as _json
    from pathlib import Path as _Path

    cfg_path = _Path(__file__).resolve().parents[2] / "scripts" / "jev_worker_config.sample.json"
    cfg = _json.loads(cfg_path.read_text())

    baseline_mm = cfg["baseline_best_ctau_mm"]
    target_mm = cfg["target_ctau_mm"]
    prov = cfg["baseline_provenance"]

    assert target_mm == pytest.approx(2.0 * baseline_mm)
    assert prov["kind"] == "canonical_replay_of_declared_reference_anchor"
    assert prov["historical_best"] is False
    assert baseline_mm != prov["historical_source_ctau_mm"]
    assert baseline_mm != prov["superseded_unverified_ctau_mm"]

    # Factor-1000 guard: a plausible ctau_mm baseline for this search sits in a
    # narrow band; a value 1000x too large or too small would fall outside it.
    assert 0.01 < baseline_mm < 100

    # No stored *_m keys should exist in the new config; state derivation must
    # compute best_ctau_m itself, from mm, never from a separately stored m.
    assert "baseline_best_ctau_m" not in cfg
    assert "target_ctau_m" not in cfg

    state = baseline_state(cfg)
    assert state["best_ctau_m"] == pytest.approx(state["best_ctau_mm"] / 1000.0)
    assert state["best_ctau_mm"] == pytest.approx(baseline_mm)

    built = build_scientific_state(cfg, state)
    assert built["objective"]["best_ctau_m"] == pytest.approx(built["objective"]["best_ctau_mm"] / 1000.0)

