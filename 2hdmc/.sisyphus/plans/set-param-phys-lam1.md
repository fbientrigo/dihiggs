# THDM: add `set_param_phys_lam1` (λ1 input instead of m12²)

## TL;DR
> **Summary**: Add `THDM::set_param_phys_lam1(...)` that reconstructs `m12_2` by *inverting the existing `set_param_phys` λ1-expression (source conventions)*, then delegates to `set_param_phys` to preserve identical internal state.
> **Deliverables**: new public setter + minimal round-trip validation program.
> **Effort**: Short
> **Parallel**: NO
> **Critical Path**: Implement setter → add round-trip validator → build + run validator

## Context
### Original Request
Implement:
```cpp
bool set_param_phys_lam1(double mh, double mH, double mA, double mHp,
                         double sin_ba, double lambda1,
                         double lambda6, double lambda7, double tan_beta);
```
Same physical point as `set_param_phys(...)`, but with **λ1** as input instead of **m12²**.

Critical rule: do **not** trust external formulas blindly; verify against 2HDMC source conventions.

### Interview / Assumptions
- No unit-test framework detected/required; validation will follow the repo’s existing **example-program** pattern (Makefile-built executable). (Default applied.)
- No additional physics constraints will be introduced; setters remain purely parameter-initialization like existing `set_param_phys`.

### Metis Review (gaps addressed)
- Add explicit, agent-executable acceptance criteria (build + run + exit codes).
- Use `get_param_gen(...)` / `get_param_phys(...)` for validation (no new public getters; `get_m12_2()` is private).
- Add explicit tolerance policy for floating comparisons.
- Avoid changing conventions: compute α,β exactly like `set_param_phys` does.

## A. Source audit (ground truth)
### Relevant files
- `src/THDM.h`
- `src/THDM.cpp`
- `src/SM.h`, `src/SM.cpp` (for `v2` meaning)
- `Makefile` (builds example programs)

### Relevant methods / APIs
- `THDM::set_param_phys(...)` — **physical basis** setter (masses, sba, λ6, λ7, m12², tanβ)
- `THDM::set_param_gen(...)` — **generic basis** setter (λ1..λ7, m12², tanβ)
- `THDM::get_param_gen(...)` — returns λ1..λ7, m12², tanβ (by reference)
- `THDM::get_param_phys(...)` — returns masses, sba, λ6, λ7, m12², tanβ (by reference)
- `THDM::get_hmass(int)` — returns physical masses by index
- `THDM::get_m12_2()` — exists but is **private**; used internally by getters

### Exact call path of `set_param_phys`
`set_param_phys` does **not** call helpers: it directly computes and stores internal state.

Diagram:
```
set_param_phys
  -> sets lambda[6], lambda[7]
  -> beta = atan(tan_beta)
  -> alpha = beta - asin(sba)
  -> computes lambda[1..5] (explicit formulas)
  -> computes m22_2 (explicit formula)
  -> sinba = sba
  -> params_set = true/false
  -> return params_set
```

### Where `m12²` enters internally
In `set_param_phys`, the input `m12_2` enters **linearly** in the formulas for:
- `lambda[1]..lambda[5]`
- `m22_2`

## B. Convention audit (source-consistent)
### v / v² convention
- `THDM::v2` is initialized from `SM::get_v2()`.
- `SM::get_v2()` is documented/implemented as:
  - \(v^2 = 1/(\sqrt{2} G_F)\) in GeV².
- Therefore **code `v2` equals physics `v²`** (not `v²/2`).

### α/β and sin(β−α) convention
In `set_param_phys`:
- `beta = atan(tan_beta)` (with required `tan_beta>0`, so β∈(0,π/2))
- `alpha = -asin(sba) + beta` (i.e. α = β − asin(sin(β−α)))
- `cba = sqrt(1 - sba*sba)` → forces **cos(β−α) ≥ 0** (no negative-branch)

**Branch/sign statement**: The source fixes the (β−α) branch via `asin(sba)` and enforces `cos(β−α)≥0` by construction (`sqrt`). There is no alternate sign/branch path.

### Candidate λ1 formula: confirmed vs corrected
**Confirmed (source match)**.

From `THDM::set_param_phys`:
\[
\lambda_1 = \frac{m_H^2\cos^2\alpha + m_h^2\sin^2\alpha - m_{12}^2\tan\beta}{v^2\cos^2\beta}
            - \frac{3}{2}\lambda_6\tan\beta + \frac{1}{2}\lambda_7\tan^3\beta
\]
with the identification **`v²` = `THDM::v2`**.

**Source-consistent inversion** (solve the code expression for `m12_2`):
\[
m_{12}^2 = \frac{m_H^2\cos^2\alpha + m_h^2\sin^2\alpha - v^2\cos^2\beta\,\big(\lambda_1 + \tfrac{3}{2}\lambda_6\tan\beta - \tfrac{1}{2}\lambda_7\tan^3\beta\big)}{\tan\beta}
\]

## C. Implementation plan (minimal diff)
### Files to modify
1) `src/THDM.h`
   - Add the declaration (plus Doxygen block) **immediately after** `set_param_phys(...)`.

2) `src/THDM.cpp`
   - Implement `THDM::set_param_phys_lam1(...)` near the existing `set_param_phys(...)`.
   - Implementation steps:
     - replicate the same basic “problematic parameter choices” checks needed to safely compute `asin(sba)` and divide by `tan_beta` (match `set_param_phys` behavior);
     - compute local `beta`, `alpha` using the same formulas as `set_param_phys`;
     - compute local `m12_2` using the **source-consistent inversion** above;
     - delegate: `return set_param_phys(..., m12_2, tan_beta);`

3) Validation harness (preferred minimal)
   - Add a new example-style program: `src/CalcRoundTrip.cpp`.
   - Build with `make CalcRoundTrip` (and optionally add it to `PROG` in `Makefile` so `make` builds it by default).

## D. Patch (planned diff — executor should apply exactly)
> Note: This is the intended minimal code change; executor should preserve existing formatting/style (tabs, line breaks) used in nearby code.

### 1) `src/THDM.h`
Add after the existing `set_param_phys(...)` declaration:

```diff
diff --git a/src/THDM.h b/src/THDM.h
--- a/src/THDM.h
+++ b/src/THDM.h
@@
 bool set_param_phys(double m_h,double m_H, double m_A, double m_Hp,
                     double sba, double lambda6, double lambda7,
                     double m12_2, double tan_beta);

+/**
+ * @brief Specifies 2HDM in the physical basis using \f$\lambda_1\f$
+ *
+ * Same physical point as set_param_phys(...), but with \f$\lambda_1\f$ as input
+ * instead of \f$m_{12}^2\f$. The method reconstructs \f$m_{12}^2\f$ using the
+ * internal 2HDMC conventions and delegates to set_param_phys(...).
+ *
+ * @param m_h  Mass of lightest CP-even Higgs \f$ h \f$
+ * @param m_H  Mass of heavier CP-even Higgs \f$ H \f$
+ * @param m_A  Mass of CP-odd Higgs \f$ A \f$
+ * @param m_Hp Mass of charged Higgs
+ * @param sba  Mixing parameter \f$ \sin(\beta-\alpha) \f$. NB: Correct sign on
+ *             \f$ \sin(\beta-\alpha) \f$ must be determined from the condition
+ *             \f$ \cos(\beta-\alpha) \geq 0 \f$
+ * @param lambda1 Value of \f$ \lambda_1 \f$ in generic potential
+ * @param lambda6 Value of \f$ \lambda_6 \f$ in generic potential
+ * @param lambda7 Value of \f$ \lambda_7 \f$ in generic potential
+ * @param tan_beta Ratio of vevs, \f$ \tan\beta=v_2/v_1 \f$
+ *
+ * @returns Boolean \a true if all parameters were set correctly, \a false otherwise
+ */
+bool set_param_phys_lam1(double m_h,double m_H, double m_A, double m_Hp,
+                         double sba, double lambda1,
+                         double lambda6, double lambda7,
+                         double tan_beta);
```

### 2) `src/THDM.cpp`
Add near `set_param_phys(...)`:

```diff
diff --git a/src/THDM.cpp b/src/THDM.cpp
--- a/src/THDM.cpp
+++ b/src/THDM.cpp
@@
 bool THDM::set_param_phys(double m_h,double m_H, double m_A, double m_Hp,
 			  double sba, double lambda6, double lambda7,
 			  double m12_2,double tan_beta) {
@@
 }

+bool THDM::set_param_phys_lam1(double m_h,double m_H, double m_A, double m_Hp,
+                               double sba, double lambda1,
+                               double lambda6, double lambda7,
+                               double tan_beta) {
+
+  if (m_h>m_H) {
+    cerr << "WARNING: Cannot set physical masses such that m_H < m_h\n";
+    params_set = false;
+    return params_set;
+  }
+
+  // Problematic parameter choices (match set_param_phys guards)
+  if ((tan_beta<=0)||(abs(sba)>1)||(m_h<0)||(m_H<0)||(m_A<0)||(m_Hp<0)) {
+    params_set=false;
+    return params_set;
+  }
+
+  // Reconstruct m12_2 by inverting the lambda[1] expression used in set_param_phys
+  double beta_loc = atan(tan_beta);
+  double cb = cos(beta_loc);
+  double cb2 = cb*cb;
+  double tb = tan(beta_loc);
+
+  double alpha_loc = -asin(sba)+beta_loc;
+  double sa  = sin(alpha_loc);
+  double sa2 = sa*sa;
+  double ca  = cos(alpha_loc);
+  double ca2 = ca*ca;
+
+  double m12_2 = (m_H*m_H*ca2+m_h*m_h*sa2
+                 - v2*cb2*(lambda1+1.5*lambda6*tb-0.5*lambda7*tb*tb*tb))/tb;
+
+  return set_param_phys(m_h,m_H,m_A,m_Hp,sba,lambda6,lambda7,m12_2,tan_beta);
+}
```

### 3) `src/CalcRoundTrip.cpp` (new)
Create a new example-style executable that hardcodes 2 test points and exits nonzero on mismatch.

```cpp
// src/CalcRoundTrip.cpp

#include "THDM.h"

#include <algorithm>
#include <cmath>
#include <iostream>

static bool nearly_equal(double a, double b, double rtol, double atol) {
  const double diff = std::fabs(a - b);
  const double scale = std::max(std::fabs(a), std::fabs(b));
  return diff <= (atol + rtol * scale);
}

static int run_case(const char* name,
                    double m_h, double m_H, double m_A, double m_Hp,
                    double sba, double lambda6, double lambda7,
                    double m12_2, double tan_beta,
                    double rtol, double atol) {

  THDM m1;
  if (!m1.set_param_phys(m_h, m_H, m_A, m_Hp, sba, lambda6, lambda7, m12_2, tan_beta)) {
    std::cout << "Case " << name << ": FAIL (set_param_phys returned false)\n";
    return 1;
  }

  double l1,l2,l3,l4,l5,l6,l7,m12_ref,tb_ref;
  m1.get_param_gen(l1,l2,l3,l4,l5,l6,l7,m12_ref,tb_ref);

  THDM m2;
  if (!m2.set_param_phys_lam1(m_h, m_H, m_A, m_Hp, sba, l1, lambda6, lambda7, tan_beta)) {
    std::cout << "Case " << name << ": FAIL (set_param_phys_lam1 returned false)\n";
    return 1;
  }

  double L1,L2,L3,L4,L5,L6,L7,m12_new,tb_new;
  m2.get_param_gen(L1,L2,L3,L4,L5,L6,L7,m12_new,tb_new);

  bool ok = true;
  ok &= nearly_equal(m12_new, m12_ref, rtol, atol);
  ok &= nearly_equal(tb_new, tb_ref, rtol, atol);
  ok &= nearly_equal(L1, l1, rtol, atol);
  ok &= nearly_equal(L2, l2, rtol, atol);
  ok &= nearly_equal(L3, l3, rtol, atol);
  ok &= nearly_equal(L4, l4, rtol, atol);
  ok &= nearly_equal(L5, l5, rtol, atol);
  ok &= nearly_equal(L6, l6, rtol, atol);
  ok &= nearly_equal(L7, l7, rtol, atol);

  for (int h = 1; h <= 4; ++h) {
    ok &= nearly_equal(m2.get_hmass(h), m1.get_hmass(h), rtol, atol);
  }

  if (ok) {
    std::cout << "Case " << name << ": PASS\n";
    return 0;
  }

  std::cout << "Case " << name << ": FAIL\n";
  std::cout << "  m12_2: ref=" << m12_ref << " new=" << m12_new << "\n";
  std::cout << "  tanb : ref=" << tb_ref  << " new=" << tb_new  << "\n";
  return 1;
}

int main() {
  const double rtol = 1e-10;
  const double atol = 1e-12;

  int rc = 0;

  // Generic point
  rc |= run_case("generic",
                 125.0, 300.0, 350.0, 360.0,
                 0.80,
                 0.10, -0.20,
                 20000.0,
                 2.0,
                 rtol, atol);

  // Near-alignment point
  rc |= run_case("near-alignment",
                 125.0, 300.0, 350.0, 360.0,
                 0.9999,
                 0.10, -0.20,
                 20000.0,
                 2.0,
                 rtol, atol);

  // Invalid-input guard (expected failure)
  {
    THDM m_bad;
    const bool ok = !m_bad.set_param_phys_lam1(125.0, 300.0, 350.0, 360.0,
                                               0.80,
                                               1.0,
                                               0.10, -0.20,
                                               -1.0);
    if (ok) {
      std::cout << "Case invalid-input: PASS\n";
    } else {
      std::cout << "Case invalid-input: FAIL\n";
      rc |= 1;
    }
  }

  if (rc == 0) {
    std::cout << "All tests passed\n";
  }
  return rc;
}
```

### 4) `Makefile` (optional but recommended)
- If the repo uses a `PROG = ...` list for the `all` target, append `CalcRoundTrip` so `make` builds it by default.

## E. Validation (agent-executed)
### Build + run
Commands (from repo root):
```bash
mkdir -p .sisyphus/evidence
make clean
make
make CalcRoundTrip
./CalcRoundTrip > .sisyphus/evidence/task-2-roundtrip-stdout.txt 2>&1
```

### Quantities compared (must match within tolerance)
- `m12_2` from `get_param_gen(...)`
- `lambda1..lambda7` from `get_param_gen(...)`
- `tan_beta` from `get_param_gen(...)`
- masses from `get_hmass(1..4)`

### Expected outcome
- `./CalcRoundTrip` prints (captured to `.sisyphus/evidence/task-2-roundtrip-stdout.txt`):
  - `Case generic: PASS`
  - `Case near-alignment: PASS`
  - `Case invalid-input: PASS`
  - `All tests passed`
- Exit code is **0**.

## Verification Strategy
> ZERO HUMAN INTERVENTION — all verification is agent-executed.
- Test decision: tests-after (example program)
- Evidence artifacts:
  - `.sisyphus/evidence/task-2-roundtrip-stdout.txt` (captured program output)

## Execution Strategy
### Parallel Execution Waves
Wave 1:
- T1 Add `set_param_phys_lam1` declaration + implementation

Wave 2:
- T2 Add round-trip validator program and run it

### Dependency Matrix
- T1 blocks T2 (validator calls the new method)

## TODOs

- [x] 1. Add `THDM::set_param_phys_lam1` public API

  **What to do**:
  - In `src/THDM.h`, declare `set_param_phys_lam1` immediately after `set_param_phys`, matching Doxygen style.
  - In `src/THDM.cpp`, implement `THDM::set_param_phys_lam1` near `set_param_phys`.
  - Implement reconstruction of `m12_2` by inverting the *exact* `lambda[1]` formula used in `set_param_phys`.
  - Delegate to `set_param_phys(...)` and return its boolean.
  - Ensure guard checks prevent invalid `asin(sba)` domain and division by nonpositive `tan_beta`, matching the existing `set_param_phys` checks.

  **Must NOT do**:
  - Do not refactor `set_param_phys`.
  - Do not add new physics constraints (unitarity, stability, etc.).
  - Do not change α/β conventions or introduce alternative sign branches.

  **Recommended Agent Profile**:
  - Category: `quick` — small, localized C++ change.
  - Skills: (none)

  **Parallelization**: Can Parallel: NO | Wave 1 | Blocks: T2 | Blocked By: —

  **References**:
  - Implementation pattern: `src/THDM.cpp:set_param_phys` (computes α=β−asin(sba), uses `v2`, defines `lambda[1]`)
  - API placement: `src/THDM.h:set_param_phys` block (declare new method directly after)
  - v² convention: `src/SM.h:get_v2()` docstring (v²=1/(sqrt(2)G_F))

  **Acceptance Criteria** (agent-executable):
  - [ ] `make` succeeds.
  - [ ] A minimal C++ snippet (or CalcRoundTrip once added) can compile/link calling `set_param_phys_lam1`.

  **QA Scenarios**:
  ```
  Scenario: Build succeeds with new API
    Tool: Bash
    Steps:
      1) mkdir -p .sisyphus/evidence
      2) make > .sisyphus/evidence/task-1-build.txt 2>&1
    Expected: Build succeeds (exit 0)
    Evidence: .sisyphus/evidence/task-1-build.txt
  ```

  **Commit**: NO (user did not request)

- [x] 2. Add round-trip consistency validator (`CalcRoundTrip`)

  **What to do**:
  - Add `src/CalcRoundTrip.cpp` (see section D.3) that:
    - initializes one model with `set_param_phys` using a chosen `m12_2`;
    - extracts `lambda1` via `get_param_gen`;
    - initializes a second model with `set_param_phys_lam1` using that `lambda1`;
    - compares (`m12_2`, λ1..λ7, tanβ, masses) with tolerance; exits nonzero on mismatch.
  - Ensure it runs two cases: generic + near-alignment.
  - Build and run it; capture stdout to evidence file.
  - Optionally append `CalcRoundTrip` to the Makefile `PROG` list if needed for `make` default builds.

  **Must NOT do**:
  - Do not add a new testing framework.
  - Do not add new public getters solely for testing.

  **Recommended Agent Profile**:
  - Category: `quick`
  - Skills: (none)

  **Parallelization**: Can Parallel: NO | Wave 2 | Blocks: Final verification | Blocked By: T1

  **References**:
  - Getter signatures: `src/THDM.h:get_param_gen`, `src/THDM.h:get_hmass`
  - Existing example style: `src/Demo.cpp`, `src/CalcPhys.cpp` (use includes/printing conventions)

  **Acceptance Criteria** (agent-executable):
  - [ ] `make CalcRoundTrip` succeeds.
  - [ ] `./CalcRoundTrip` exits 0 and prints PASS for both physics cases + invalid-input guard.

  **QA Scenarios**:
  ```
  Scenario: Round-trip equivalence (generic + near-alignment)
    Tool: Bash
    Steps:
      1) mkdir -p .sisyphus/evidence
      2) make CalcRoundTrip
      3) ./CalcRoundTrip > .sisyphus/evidence/task-2-roundtrip-stdout.txt 2>&1
    Expected:
      - Output contains:
        - "Case generic: PASS"
        - "Case near-alignment: PASS"
        - "Case invalid-input: PASS"
        - "All tests passed"
      - Exit code 0
    Evidence: .sisyphus/evidence/task-2-roundtrip-stdout.txt
  ```

  **Commit**: NO (user did not request)

## Final Verification Wave (MANDATORY)
> Run after ALL implementation tasks. Do NOT auto-proceed; wait for user approval.

- [x] F1. Plan Compliance Audit — oracle

  **Tool**: Bash + Read

  **Steps**:
  1) `mkdir -p .sisyphus/evidence`
  2) Capture changed-file list:
     - `git diff --name-only > .sisyphus/evidence/final-f1-files.txt`
  3) Confirm the new API exists in header + implementation:
     - `grep -n "set_param_phys_lam1" src/THDM.h src/THDM.cpp > .sisyphus/evidence/final-f1-grep.txt`
  4) Inspect `src/THDM.cpp` to confirm `m12_2` is reconstructed by inverting the *existing* `lambda[1]` expression in `set_param_phys` (same coefficients/signs; same `alpha=beta-asin(sba)` convention).

  **Expected**:
  - Only intended files changed (typically `src/THDM.h`, `src/THDM.cpp`, `src/CalcRoundTrip.cpp`, optional `Makefile`, plus `.sisyphus/evidence/*`).
  - `set_param_phys_lam1` delegates to `set_param_phys` and preserves conventions.

  **Evidence**:
  - `.sisyphus/evidence/final-f1-files.txt`
  - `.sisyphus/evidence/final-f1-grep.txt`

- [x] F2. Code Quality Review — unspecified-high

  **Tool**: Bash + Read

  **Steps**:
  1) `mkdir -p .sisyphus/evidence`
  2) Build from a clean state and capture output:
     - `make clean && make > .sisyphus/evidence/final-f2-build.txt 2>&1`

  **Expected**:
  - Exit code 0.
  - No new warnings/errors attributable to the change (review the build log).

  **Evidence**:
  - `.sisyphus/evidence/final-f2-build.txt`

- [x] F3. Real Manual QA — unspecified-high

  **Tool**: Bash

  **Steps**:
  1) `mkdir -p .sisyphus/evidence`
  2) Run the round-trip validator and capture output:
     - `./CalcRoundTrip > .sisyphus/evidence/final-f3-roundtrip.txt 2>&1`

  **Expected**:
  - Exit code 0.
  - Output includes the PASS lines for generic + near-alignment + invalid-input.

  **Evidence**:
  - `.sisyphus/evidence/final-f3-roundtrip.txt`

- [x] F4. Scope Fidelity Check — deep

  **Tool**: Bash + Read

  **Steps**:
  1) Confirm `set_param_phys` was not refactored (beyond inserting the new method nearby).
  2) Confirm no new constraint logic was introduced (no new physics checks beyond the existing `set_param_phys`-style guards).

  **Expected**:
  - Diff is minimal and limited to the requested feature + validator.

  **Evidence**:
  - `.sisyphus/evidence/final-f4-notes.txt` (brief bullet summary)

## Commit Strategy
- No commits unless the user explicitly requests.

## Success Criteria
- New API `set_param_phys_lam1(...)` exists and compiles.
- Round-trip validation (`CalcRoundTrip`) passes for both a generic and near-alignment point.
- A model initialized with `set_param_phys_lam1(...)` is numerically equivalent (within tolerance) to one initialized with `set_param_phys(...)` using the reconstructed `m12_2`.
