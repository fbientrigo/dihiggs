# DiHiggs documentation

This site documents the current **2HDMC / DiHiggs evaluation core**. The active code
is centered on the row-preserving v2 evaluators and the modular orchestrator.
Historical ML, adaptive-search, and replay components remain in the repository only
where their status is explicitly documented.

## Start here

```{toctree}
:maxdepth: 2
:caption: Current use

quickstart
cli/orchestrator
artifacts
hpc
code_reference
documentation_map
```

## Physics conventions

```{toctree}
:maxdepth: 2
:caption: Physics authority

theory/2hdm_conventions
PROJECT_PHYSICS_AUTHORITY
PHYSICS_AUTHORITY_DISCREPANCY_LEDGER
```

## Canonical evaluators and contracts

```{toctree}
:maxdepth: 2
:caption: Contracts

contracts/canonical_evaluators_v2
OPERATING_STATUS_V2
HIGH_MASS_H2_CONTRACT
contracts/replay_safe_scan_output_contract
UFO_GENERICIZATION_REQUIREMENTS
DOWNSTREAM_INTERFACE_GAP_REPORT
```

## Architecture, verification, and frozen history

```{toctree}
:maxdepth: 2
:caption: Project records

WORKSPACE_ARCHITECTURE_2026
COMPUTE_SCALE_PLAN
inter_repo_communication_assessment
scan_harness_architecture_v2
characterization_lambda1
verification/dihiggs_point_v2_verification_v1
verification/lambda1_v2_yukawa_fix_v1
campaigns/yukawa_order_recompute_v1/README
campaigns/yukawa_order_recompute_v2/README
handoffs/REPO_TO_TEX_NOTES_V3
migration/data_contract_v0.1_draft
legacy_and_experimental
autoresearch_frozen
REPOSITORY_CLOSURE_V2
audits/closure_2026-07/REPOSITORY_CLOSURE_INVENTORY_2026-07
audits/closure_2026-07/boundary_alignment_report
audits/closure_2026-07/postmerge_verification
audits/closure_2026-07/pr58_phase3_review
audits/closure_2026-07/pr60_merge_gate_review
audits/lambda1_lifetime_audit_v1
campaigns/yukawa_order_recompute_v1/PHASE_A_INVENTORY
site/SITE_SPEC
site/DESIGN_SPEC
```
