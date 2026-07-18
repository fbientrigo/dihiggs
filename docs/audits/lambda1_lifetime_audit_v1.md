# Historical lambda1 lifetime serialization audit

This read-only audit selects only the exact legacy fixed/15 schema and excludes
recomputed, backup, recovered, artifact, and v2 paths. The historical evaluator
used explicit-in-code `mh = 125.0 GeV`; this is provenance, not a v2 default.
Audited roots: /home/fabi/dihiggs/runs.

| Metric | Value |
|---|---:|
| Files | 2020 |
| Rows | 2020 |
| Serialized zero widths | 0 |
| Minimum positive width (GeV) | 3.29838342800000007e-06 |
| Apparent maximum ctau (mm) | 5.98253613345525178e-08 |
| Recoverable coordinate rows | 2020 |
| Deterministic replay eligible | 2020 |
| Autoresearch `width > 0` discards | 0 |
| Positive-width ranking boundary rows | 20 |

## Campaign classification

| Classification | Rows |
|---|---:|
| autoresearch-supported | 20 |
| diagnostic-discard | 2000 |

## Interpretation

A serialized zero cannot recover a width from the CSV alone. Replay is eligible only
when all legacy lambda1 coordinates are finite and the historical `mh` provenance is
known. Re-run only zero-width or positive-width ranking-boundary rows; never a complete grid.
Aggregate source checksums and campaign classifications are in the JSON manifest.
