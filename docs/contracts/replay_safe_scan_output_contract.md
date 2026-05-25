# Replay-safe scan output contract (additive)

Scope: additive metadata/output hardening for C++ scan engines. This contract does not change physics formulas, scan semantics, model construction order, or historical CSVs.

## 1) Legacy m12 status

- `m12` is a legacy compatibility column.
- `m12` is replay-unsafe unless explicit replay semantics are provided with it.
- Consumers must not assume that legacy `m12` alone is sufficient for CalcPhys parity.

## 2) Replay-safe core fields

Required fields for replay-safe interpretation:

1. `m12_2_used`
   - Canonical value used in `set_param_phys(...)` for replay parity.
   - This is the value that must be used for replay.

2. `m12_2_gen_after_set`
   - Generic-basis `m12^2` read back from the canonical model after `set_param_phys(...)`.

3. `delta_m12_2_gen_minus_used`
   - Defined as `m12_2_gen_after_set - m12_2_used`.

4. `replay_semantics_version`
   - Semantic tag for interpretation of replay-safe fields.
   - Current helper default: `pswf_m12_used_v1`.

## 3) Provenance/context fields

Required metadata fields for traceability:

- `yukawa_type`
- `higgs_state`
- `model_api_path`
- `calc_engine`
- `git_commit` (if available)
- `git_dirty` (if available)

Recommended default conventions for PhysScanWithFixings path:

- `yukawa_type = 1`
- `higgs_state = h2/DECAY35` (or project-equivalent explicit convention)
- `model_api_path = THDM::set_param_phys_lam1->THDM::get_param_gen->THDM::set_param_phys`
- `calc_engine = 2HDMC::DecayTable`

## 4) Precision requirements

- Floating-point output must use scientific notation with 17 digits of precision where possible.
- Replay-safe metadata fields (`m12_2_used`, `m12_2_gen_after_set`, `delta_m12_2_gen_minus_used`) must be emitted at this precision.
- Rationale: lower precision can reintroduce replay drift/catastrophic parity mismatch.

## 5) Additive-only compatibility rules

- Do not remove legacy columns (including `m12`).
- Add replay-safe columns without deleting historical fields.
- Do not recompute or overwrite historical database artifacts as part of this refactor branch.
- Do not claim global historical database validation from this contract alone.
