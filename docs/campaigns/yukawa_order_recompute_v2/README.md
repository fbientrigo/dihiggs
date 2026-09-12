# Issue #62 fresh recomputation

This is a fresh, bounded recomputation using the maintained
`DihiggsPointV2Evaluator` and its current `dihiggs.point.v2` schema. Two known
historical points seed neighborhoods; two independent families test behavior
away from those seeds. Seed values are metadata only and do not constrain
acceptance.

The recorded 102-point run contains accepted observations in both off-seed
families; rejected neighbors and construction failures remain in the output as
evidence rather than being filtered out.

Run after building the evaluator:

```bash
python3 scripts/run_yukawa_order_recompute_v2.py
```

The manifest records evaluator/Yukawa-fix provenance, the exact input grid,
new-artifact integrity hashes, row accounting, and physical validation. The
campaign deliberately does not produce an old-vs-new checksum or numeric
comparison: historical artifacts have different column layouts and are not a
valid same-schema scientific baseline. Broad historical scans without exact
coordinates remain quarantined as `UNRESOLVABLE_PROVENANCE`.
