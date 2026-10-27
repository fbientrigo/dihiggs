# Physics authority and downstream migration ledger

Status date: 2026-10-21
Authority: [`conventions/physics_conventions.yaml`](../conventions/physics_conventions.yaml)  
Authority snapshot: `dihiggs@087ba246e9a9fa34b16249b52629a62703893b5c`, SHA-256 `4f3a938ed54dbc779d25f691ce1dd06fe17af69751105dc3b69ccb6b7d92fff6`.

The canonical YAML is owned by `dihiggs`. Downstream repositories carry byte-identical caches and `conventions/physics_conventions.lock.yaml` provenance sidecars. A cache cannot become authoritative because it is newer or locally edited. This ledger tracks migration status and historical conflicts; it does not authorize changing frozen outputs.

## Current convention

The active v3 value for new production is the decimal string `m_h = "125.20"` GeV, convention ID `mh_125p20_active_2026`. This is an adopted project convention, not an external measurement attribution. The source repository's committed YAML bytes are the checksum authority. The current workspace pins that source commit.

The `physical_2HDM` trilinear is signed. A v3 consumer uses `g_physical_hphiphi_GeV`, verifies consistency with the raw 2HDMC fields and magnitude compatibility field, and maps `GHphiphi = -g_physical_hphiphi_GeV`. The two effective prescriptions have unresolved UFO mappings and remain unsupported in sign-sensitive production. Missing or inconsistent metadata fails closed.

## Migration status

| Repository | Current status | Remaining boundary |
|---|---|---|
| `dihiggs` | Authority owner; v3 canonical YAML and evaluator provenance | Preserve the canonical file and its attestation; reconcile further changes only through its owner. |
| `dihiggs_boundary` | v3 controlled cache, sidecar, and checksum lock; new adaptive/profiled handoffs use core provenance | The `evaluate_point_v1` / Atlas v0 runner and its 125.09 HiggsSignals reference are replay-only. Do not change its golden outputs. |
| `dihiggs_hep_cross` | v3 cache and sidecar; constants loader fails closed; driver/banner use signed physical mapping for provenance-bearing v3 points | Effective prescriptions remain unsupported; no MadGraph results are promoted. |
| `dihiggs_ufo` | v3 controlled cache and provenance sidecar | Frozen Pack B builders remain historical; add no v3/effective mapping claim without a separately validated UFO change. |
| `dihiggs_llp_recast` | v3 controlled cache and provenance sidecar | New sample manifests must retain incoming point convention provenance; frozen response results remain historical. |

## Historical preservation

The old migration inventories below were generated against earlier branch snapshots and must not be treated as current instructions. They are retained in Git history only; their assertions that the authority should return to 125.13 or that 125.20 was a PDG measurement are superseded.

| Convention ID | Value | Scope |
|---|---:|---|
| `mh_125p13_pdg_2026_historical_replay` | 125.13 | Completed adaptive M2 campaign and explicit replay fixtures |
| `mh_125p13_high_mass_v1` | 125.13 | Historical high-mass and 150 GeV benchmark chain |
| `mh_125p09_lambda1_boundary_legacy` | 125.09 | Legacy Boundary characterization/replay |
| `mh_125p0_ufo_autoresearch_legacy` | 125.0 | Frozen upstream UFO and autoresearch replay paths |

These values keep their original bytes, hashes, and labels. Boundary v0 and old UFO builder paths are not silently migrated to v3. The current MadGraph result registry remains empty until validation and explicit human promotion.

## Acceptance rule for future migrations

Every downstream cache must hash to the exact canonical source bytes and have a sidecar naming the source repository, full commit, path, checksum, schema, and copy date. Point handoffs must include the v3 schema, source commit, and source SHA-256. Consumers reject missing, stale, or inconsistent provenance. New campaign assumptions such as alignment, `lambda7`, and mass degeneracy remain explicit inputs.
