#!/usr/bin/env python3
"""Attach the completed actual-calculator scan to the campaign reports."""
from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
CAMPAIGN = HERE.parent
REPO = CAMPAIGN.parents[1]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


summary = json.loads((CAMPAIGN / "calculator_scan_summary.json").read_text())
records = [json.loads(line) for line in (CAMPAIGN / "calculator_scan_agreement.jsonl").read_text().splitlines() if line]
assert summary["status"] == "COMPLETE" and len(records) == 3510
status_counts = {
    convention: {
        mode: dict(Counter(row["constructor_status"] for row in records
                           if row["mass_convention"] == convention and row["mode"] == mode))
        for mode in ("legacy", "explicit", "analytic")
    }
    for convention in ("mh125.13_T1", "mh125.20_active")
}
failure_counts = dict(Counter((row["mode"], row["exact_failure_reason"])
                              for row in records if row["constructor_status"] != "PASS"))
route_call_wall_sum = sum(float(row.get("calculator_elapsed_seconds", 0.0)) for row in records)
scan_evidence = {
    "status": "COMPLETE",
    "pilot_points": summary["pilot_point_count"],
    "route_precision_records": summary["route_precision_record_count"],
    "precision_digits": [50, 80, 120],
    "modes": ["legacy", "explicit", "analytic"],
    "status_counts": status_counts,
    "exact_failure_reason_counts": {
        f"{mode}: {reason}": count for (mode, reason), count in failure_counts.items()
    },
    "actual_backend_spectrum_residual_max": summary["maximum_actual_backend_spectrum_residual_by_route_precision"],
    "summed_route_call_wall_seconds": route_call_wall_sum,
    "execution": {
        "workers": 1,
        "batch_size_points": 10,
        "wrapper_waiver": "User explicitly waived lab-health/safe-run for this bounded scan.",
        "enforced_process_limits": "prlimit RLIMIT_AS=9 GiB, RLIMIT_CPU=900 s; OMP_NUM_THREADS=1; OPENBLAS_NUM_THREADS=1",
        "resource_profile": "calculator_scan_profile.json is a one-point process profile; route-call time sum is not end-to-end wall time",
    },
    "prelock_collision_audit": {
        "input_records": 5292,
        "unique_records": 3510,
        "identical_duplicate_rows_removed": 1782,
        "scientific_payload_differences": 0,
        "input_snapshot": "snapshots/concurrent_writer_collision_20261008T/calculator_scan_agreement_5292_lines.jsonl",
        "deduplicated_snapshot": "snapshots/concurrent_writer_collision_20261008T/calculator_scan_agreement_deduplicated.jsonl",
        "sha256_manifest": "snapshots/concurrent_writer_collision_20261008T/SHA256SUMS.txt",
        "interpretation": "Pre-lock observations remain under their original identity and are not relabeled as v4 records.",
    },
    "identity": summary["scan_identity"],
    "jsonl": "calculator_scan_agreement.jsonl",
    "jsonl_sha256": sha(CAMPAIGN / "calculator_scan_agreement.jsonl"),
    "summary": "calculator_scan_summary.json",
    "checkpoint": "calculator_scan_checkpoint.json",
    "interpretation": "Constructor-agreement PASS/FAIL describes numerical agreement with the supplied family after CP-even state and phase mapping. FAIL is not a physical-validity classification; all points still require separate physical constraints and certified extrema work.",
}

legacy_regression = json.loads((CAMPAIGN / "legacy_decay_regression.json").read_text())
scan_evidence["default_decay_regression"] = {
    "status": legacy_regression["status"],
    "anchors": {name: {"status": item["status"],
                        "max_absolute_difference": item["max_absolute_difference"],
                        "mismatch_count": item["mismatch_count"]}
                for name, item in legacy_regression["anchors"].items()},
    "evidence": "legacy_decay_regression.json",
}

for name in ("T4_integration_report.json", "numerical_precision_audit.json"):
    path = CAMPAIGN / name
    obj = json.loads(path.read_text())
    obj["calculator_scan_agreement"] = scan_evidence
    if name == "T4_integration_report.json":
        obj["production_dispatch"] = {
            "status": "PASS", "modes": ["legacy", "explicit", "analytic"],
            "default_mode": "legacy",
            "evidence": "test_orchestrator/test_t4_cli.py covers dry-run no-backend, tampering rejection, and actual analytic dispatch",
            "legacy_schemas": "non-T4 engines retain their existing defaults and output schemas",
        }
        obj["cp_even_state_mapping"] = {
            "status": "PASS",
            "method": "minimum requested-mass assignment with target Higgs-basis phase convention; raw and mapped indices/projections/couplings are retained",
            "degenerate_state_policy": "INCONCLUSIVE; widths NOT_TESTED",
            "evidence": "reordered-state and signed-coupling tests plus pilot-wide records",
        }
        analytic_rows = [x for x in records if x.get("mode") == "analytic" and x.get("constructor_ok")]
        obj["native_diphoton_component_audit"] = {
            "status": "PASS" if analytic_rows and all(x.get("independent_diphoton_components", {}).get("status") == "PASS" for x in analytic_rows) else "INCONCLUSIVE",
            "method": "independent mpmath reconstruction of per-fermion, W, and charged-Higgs terms, signed mapped amplitudes, running masses, and the diagnostic loop width",
            "qcd_off_policy": "legacy production width remains NOT_TESTED when 2HDMC disables the channel with QCD off; the separate loop diagnostic is still compared",
        }
        appendix = [x.get("independent_diphoton_components", {}).get("appendix_LO_vs_running_mass_QCD_width") for x in analytic_rows]
        appendix = [x for x in appendix if x]
        if appendix:
            obj["diphoton_amplitude_crosscheck"] = {
                "status": "INCONCLUSIVE",
                "reason": "native 2HDMC amplitudes and widths pass component reconstruction; the leading-order pole-scheme appendix comparison remains a scheme diagnostic",
                "appendix_LO_vs_running_mass_QCD_width_relative_difference": appendix[0]["relative_difference"],
                "staged_differences": "appendix pole-mass LO, native running-mass before top QCD, and top-QCD width shift are recorded per route",
            }
    path.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n")

status_path = CAMPAIGN / "T2_T4_STATUS.md"
status = status_path.read_text()
partial_start = status.index("The work is partial.")
precision_start = status.index("The precision audit", partial_start)
status = status[:partial_start] + (
    "The work is partial. The pilot-wide actual-calculator comparison completed all 390 points "
    "through legacy, explicit, and analytic routes at 50/80/120-digit input precision. Analytic "
    "constructor agreement passed all 1170 route/precision records; legacy and explicit routes "
    "failed their recorded per-observable family-agreement checks, with exact reasons and residuals "
    "preserved. The independent native diphoton component reconstruction passed; the leading-order "
    "appendix versus running-mass/QCD width comparison remains `INCONCLUSIVE` as a scheme diagnostic. "
    "Full stationary-point enumeration and certified boundary refinement remain incomplete. The "
    "additive `--engine t4` CLI supports legacy, explicit, and analytic modes while preserving existing "
    "engine defaults and schemas. Nondegenerate reversed CP-even states are mapped to calculator "
    "indices with raw/mapped phase provenance; exactly degenerate states remain `INCONCLUSIVE` with "
    "widths `NOT_TESTED`. Displaced-photon selection and experimental exclusions remain `NOT_TESTED`.\n\n"
) + status[precision_start:]
status = status.replace(
    "Heavy execution remains `BLOCKED`: `lab-health` and `safe-run` were unavailable.",
    "The user explicitly waived the unavailable `lab-health` and `safe-run` wrappers for this campaign. The calculator-agreement scan ran serially with `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, `prlimit --as=9663676416 --cpu=900`, a nonblocking `fcntl.flock` held across all output reads/writes, and 10-point resumable checkpoints. Full certification and boundary-refinement jobs remain unrun because those scientific phases are incomplete."
)
status = status.replace(
    "Full stationary-point enumeration, certified boundary refinement, production-CLI exposure, displaced-photon selection, and experimental exclusions remain incomplete or `NOT_TESTED`. The analytic calculator API reports `INCONCLUSIVE` for `mphi <= mh` until it can map reordered or degenerate CP-even states to decay-table indices.",
    "Full stationary-point enumeration and certified boundary refinement remain incomplete. The additive `--engine t4` production CLI supports legacy, explicit, and analytic modes. Displaced-photon selection and experimental exclusions remain `NOT_TESTED`. Reordered nondegenerate CP-even states are mapped to calculator decay indices with explicit phase provenance; exactly degenerate states remain `INCONCLUSIVE` with widths `NOT_TESTED`."
)
status = status.replace(
    "Original library/build outputs were left intact.",
    "The canonical `2hdmc/lib/lib2HDMC.a` and additive `2hdmc/analytic-bridge` were rebuilt with all five library translation units for ABI consistency. The existing `DihiggsPointV2Evaluator` and `Lambda1EvaluatorV2` executables retained their baseline SHA256 values; existing legacy calculator executables were not relinked. The campaign-local relinked evaluator was compared at historical 450/500 GeV and an active-mass anchor, with exact CSV equality."
)
new_status_section = (
    "\n\nThe completed actual-calculator comparison is recorded in "
    "`calculator_scan_agreement.jsonl` (3510 records: 390 points × 3 precisions × 3 routes), "
    "with replay identity/checksum in `calculator_scan_checkpoint.json` and aggregation in "
    "`calculator_scan_summary.json`. The T1 and active-mass slices remain separate. The "
    "scan's constructor-agreement `FAIL` statuses report residuals to the stable-family "
    "target; they do not mean the point is physically invalid. Analytic states use the "
    "requested CP-even eigenstate and phase map; explicit and legacy failures retain "
    "their measured soft-term/quartic/scale residuals. Earlier passes and interrupted "
    "pre-audit attempts are preserved under `snapshots/` or labeled `.pre_`/`.aborted_`; "
    "only the current identity-matching output is summarized. Certified extrema "
    "enumeration and boundary refinement remain incomplete.\n"
)
status = status.split("\n\nThe completed actual-calculator comparison is recorded in ",1)[0] + new_status_section
status += (
    "\nThe locked v4 scan identity is `" + summary["scan_identity"]["fingerprint"] + "`; "
    "the previous pre-lock v3 dataset remains separately checksummed under "
    "`snapshots/pre_scan_lock_v3_20261008T/`. The current runner acquires the lock before "
    "loading or mutating scan files and rejects a concurrent writer with the owner PID and "
    "scan identity. A fresh 390-point v4 pass and a strict no-new-calculator-call replay both "
    "completed successfully.\n"
)
status += (
    "\nThe reproducible v4 scan and strict replay commands (run separately; both acquire the "
    "same nonblocking lock) are:\n\n"
    "```sh\n"
    "prlimit --as=9663676416 --cpu=900 -- env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \\\n    /tmp/dihiggs-t2-t4-venv/bin/python campaigns/t2_t4_20261008T043858Z/scripts/run_t4_scan.py --batch-size 10\n"
    "```\n\n"
    "The first invocation produced the current v4 records; the second validated all record "
    "identities, inputs, schemas, and payload checksums and emitted the summary without additional "
    "calculator calls. The lock helper's concurrent-writer behavior is covered by "
    "`test_scan_lock_blocks_a_concurrent_writer_and_records_owner`.\n"
)
status_path.write_text(status)

report_path = CAMPAIGN / "T4_integration_report.md"
report = report_path.read_text()
scope_start = report.index("## Scope and route contract")
calculator_start = report.index("## Calculator agreement")
scope_section = (
    "## Scope and route contract\n\n"
    "The [orchestration API](../../dihiggs/app/orchestrator/analytic_bridge.py) defaults to `mode=\"legacy\"`. "
    "The additive production CLI supports `--engine t4 --t4-mode legacy|explicit|analytic --point-json ... --output ...`; "
    "existing engines retain their prior default and schemas. The bridge uses `THDM::set_param_phys` for legacy, "
    "`THDM::set_param_phys_lam1` for explicit, and `THDM::set_param_gen_exact` for analytic construction. Explicit "
    "and analytic modes match the requested electroweak scale through `SM::set_GF`; legacy keeps existing SM inputs. "
    "Nondegenerate reversed CP-even states are mapped to actual calculator indices and the target phase convention; "
    "raw and mapped states/couplings are retained. Exactly degenerate states return `INCONCLUSIVE` without decay calls.\n\n"
    "The anchor comparisons use (mh=125.13 GeV, mphi=150 GeV, mA=mHp=200 GeV, tanβ=10^6, "
    "lambda1_ref=1.2, X=1, v=246 GeV). Full decimal and binary64 inputs, observable residuals, "
    "tolerances, spectra, and output spreads are retained in `T4_integration_report.json`.\n\n"
)
report = report[:scope_start] + scope_section + report[calculator_start:]
decay_start = report.index("## Decays and limitations")
pilot_start = report.index("## Pilot-wide calculator agreement")
decay_section = (
    "## Decays and limitations\n\n"
    "For the analytic anchor, total width is 4.1119941e-15 GeV, diphoton branching fraction is "
    "1.6567e-4, and cτ is 47.99 mm. The JSON records partial widths and branching fractions for the "
    "supported channels. Native diphoton diagnostics export per-fermion, W, and charged-Higgs terms, "
    "signed amplitudes, running-mass inputs, QCD settings, and loop width. Independent component and "
    "width reconstruction passes. The appendix pole-mass LO width differs from the native running-mass "
    "plus top-QCD width by about 8.2% at this anchor; that scheme comparison remains `INCONCLUSIVE`. "
    "When QCD-off disables the legacy production diphoton channel, its width is `NOT_TESTED` and a "
    "separate loop diagnostic is retained. Decay widths and lifetime are diagnostic, not validated "
    "phenomenology.\n\n"
    "Displaced-photon selection and experimental exclusions are `NOT_TESTED`. No benchmark is listed "
    "in `validated_benchmarks.json`. The campaign rebuilt the canonical 2HDMC library and additive "
    "bridge with all five library translation units; existing `DihiggsPointV2Evaluator` and "
    "`Lambda1EvaluatorV2` executables retained their baseline SHA256 values, and existing legacy "
    "calculator executables were not relinked. A campaign-local evaluator was compared against the "
    "baseline evaluator at historical 450/500 GeV and an active-mass anchor, with exact CSV equality.\n\n"
)
report = report[:decay_start] + decay_section + report[pilot_start:]
new_report_section = (
    "\n\n## Pilot-wide calculator agreement\n\n"
    "The campaign now contains an identity-bound serial comparison over all 300 T1 "
    "and 90 active-mass pilot points at 50, 80, and 120 digit construction precision "
    "through legacy, explicit, and analytic calculator routes (3510 records total). "
    "Each successful constructor record includes all λ₁…λ₇ and Z₁…Z₇, backend double "
    "inputs, residuals, the four complete actual-backend scattering spectra and "
    "independent spectra, signed couplings, and widths/branching fractions/lifetime. "
    "The exact constructor failure reason is retained per route/precision record. "
    "The complete results are in `calculator_scan_agreement.jsonl`; status and "
    "maximum spectrum residuals are in `calculator_scan_summary.json`. These "
    "constructor-agreement labels are not physical-validation statuses. The scan "
    "shows analytic states assigned and phase-mapped from the actual backend spectrum; "
    "constructor-agreement statuses include per-observable closure checks. The native "
    "diphoton audit compares each fermion, W, and charged-Higgs term and the running-mass/QCD "
    "loop width. Default evaluator CSV rows are unchanged at historical450, historical500, "
    "and an active-mass representative. These results do not certify global extrema or "
    "physical boundaries.\n"
)
report = report.split("\n\n## Pilot-wide calculator agreement",1)[0] + new_report_section
report_path.write_text(report)

# Refresh hashes for all campaign artifacts, excluding the manifest itself.
manifest_path = CAMPAIGN / "reproducibility_manifest.json"
manifest = json.loads(manifest_path.read_text())
manifest["campaign_file_sha256"] = {
    str(path.relative_to(CAMPAIGN)): sha(path)
    for path in sorted(CAMPAIGN.rglob("*"))
    if path.is_file() and path != manifest_path and "__pycache__" not in path.parts and path.suffix != ".pyc"
}
manifest["calculator_scan_agreement"] = scan_evidence
manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
print(json.dumps({"scan_status": "COMPLETE", "records": len(records),
                  "status_counts": status_counts, "summed_route_call_wall_seconds": route_call_wall_sum,
                  "manifest_artifacts": len(manifest["campaign_file_sha256"])}, indent=2))
