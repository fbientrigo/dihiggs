# Learnings — set_param_phys_lam1

## [2026-04-03T13:41] Session started
- Plan: Add THDM::set_param_phys_lam1(λ1 input instead of m12²)
- Execution strategy: Sequential (Wave 1 → Wave 2 → Final Verification)
- Critical: Preserve source conventions exactly (α=β−asin(sba), v2=v²)

## [2026-04-03 TASK 1] set_param_phys_lam1 Implementation Complete

- **Declaration**: Added to src/THDM.h line 209-212 with Doxygen documentation block (lines 189-208)
- **Implementation**: Added to src/THDM.cpp lines 355-389 immediately after set_param_phys (line 352)
- **Build Status**: SUCCESS (exit code 0, no errors)
- **Build Evidence**: `.sisyphus/evidence/task-1-build.txt` (11 lines, all compilation succeeded)

### Implementation Details
- Properly replicated both guard-check blocks from set_param_phys (lines 360-370)
- Computed local variables: beta_loc, cb, cb2, tb, alpha_loc, sa, sa2, ca, ca2 (lines 372-381)
- Reconstructed m12_2 using exact inversion formula from lambda[1] derivation (lines 383-385)
- Delegated to set_param_phys with reconstructed m12_2 (line 388)
- Used indentation style consistent with existing code (tabs)

### Formatting Observations
- Existing codebase uses tab indentation
- Lambda[1] formula inverts cleanly: `m12_2 = (m_H²ca²+m_h²sa² - v²cb²(λ1+1.5λ6·tb-0.5λ7·tb³))/tb`
- Convention verified: `alpha = -asin(sba) + beta` (not alternative branches)
- v2 is a member variable, accessible within method scope


## Task 2: Round-Trip Validator (CalcRoundTrip) — COMPLETED

### Deliverable
- **File**: `src/CalcRoundTrip.cpp` (3.1 KB)
- **Executable**: `CalcRoundTrip` (235 KB, dynamically linked)
- **Evidence**: `.sisyphus/evidence/task-2-roundtrip-stdout.txt`

### Test Results
All 3 test cases **PASS** with exit code 0:
1. **Generic**: Standard physics point (m_h=125, m_H=300, sin(β−α)=0.80, tan_β=2.0)
2. **Near-alignment**: Extreme mixing (sin(β−α)=0.9999, all other params same)
3. **Invalid-input guard**: Rejects tan_β=-1.0 as expected (returns false)

### Build
- Compiled with: `g++ src/CalcRoundTrip.cpp -Isrc -std=c++11 -Wall -O3 -lgsl -lgslcblas -Llib -l2HDMC -lgsl -lgslcblas -lm -o CalcRoundTrip`
- No warnings or errors
- Linked successfully against libTHDM

### Validation Logic
- Uses `nearly_equal(a, b, rtol=1e-10, atol=1e-12)` for floating-point comparisons
- Compares after round-trip:
  - λ1..λ7 (7 couplings)
  - m12² (soft mass parameter)
  - tan_β (vev ratio)
  - Higgs masses (via `get_hmass(1..4)`)
- Tests confirm `set_param_phys_lam1` reproduces identical internal state to `set_param_phys`

### Key Insight
The validator confirms that inverting the λ1 expression from `set_param_phys` and delegating back to `set_param_phys` produces bit-level consistency (within numeric tolerance), validating the inversion formula and delegation pattern.

## [2026-04-03 TASK F4] Scope Fidelity Check

- `set_param_phys` implementation remains unchanged; no refactor was applied to its body.
- New method `set_param_phys_lam1` is additive and delegates to `set_param_phys` after reconstructing `m12_2`.
- No new physics-constraint framework was introduced; only existing guard pattern was reused.
