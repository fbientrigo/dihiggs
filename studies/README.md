# Studies

Code of the reports that were generated at the workspace root (`reports/<name>/`). The scripts read `$DIHIGGS_ROOT` (set by `dihiggs_hep_cross/scripts/env.sh`; defaults to the workspace root three levels above the script) and write their outputs into `$DIHIGGS_ROOT/reports/<name>/`, which is output-only and stored in the lake (restore: `python lake/lake.py restore 20261028T024759Z-050cde reports`).

| Study | Order | Notes |
|---|---|---|
| `lambda7_tables_20261014/` | 1 | `build_tables.py`, `recompute_decays.py`, MG5 runners (`run_mg5.sh` needs `PY311`) |
| `llp_families_20261019/` | 2 | `core.py` is imported by the cτ study |
| `ctau_control_20261020/` | 3 | `finalize_ctau.py` reproduces `selected_points.csv`, `raw/ranges.json`, `tex/table_ctau_mixed250.tex` bit-for-bit from `raw/*.jsonl` |
