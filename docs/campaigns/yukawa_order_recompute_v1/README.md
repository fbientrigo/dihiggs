# Issue #62 bounded replay

This v1 replay is retained as historical engineering evidence only. It is
superseded for scientific recomputation by
`docs/campaigns/yukawa_order_recompute_v2/`, which uses fresh neighborhoods and
does not treat cross-schema historical comparisons as constraints.

This campaign replays the five canonical `dihiggs.point.v2` engineering-anchor
cases whose exact post-fix inputs are versioned in `input_grid.json`. It calls
the maintained `DihiggsPointV2Evaluator` and records a new output path and
checksums; the frozen pilot snapshots are not overwritten.

The run is a bounded pilot, not a claim about the unavailable historical scan
space. The Phase A inventory found 215 legacy scan payloads and recovery logs
without local Git-LFS payloads, per-row coordinates, or complete run manifests.
Those remain `UNRESOLVABLE_PROVENANCE` and must not be regenerated from file or
directory names.

Run after a clean build:

```bash
python3 scripts/run_yukawa_order_recompute_v1.py
```

The manifest records the exact producer commit, dirty state, Yukawa-fix
provenance, input/output checksums, schema, timing, row accounting, and
validation checks. `campaign_manifest.json` and `point_v2_anchors.csv` are the
bounded evidence products; per-case CSV/stdout/stderr files are retained for
replay diagnostics.
