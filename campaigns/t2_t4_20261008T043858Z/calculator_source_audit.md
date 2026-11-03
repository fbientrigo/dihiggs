# Calculator source audit (read-only)

This campaign is exploratory and does not alter scientific authorities.

- `2hdmc/src/THDM.cpp:236`: `set_param_gen` accepts all seven stable quartics directly. It constructs sin/cos using tan beta, but stores beta through atan; downstream rotations require a measured precision audit. It reconstructs m22 from stationarity. Mass ordering and nearly degenerate CP-even masses require explicit handling.
- `2hdmc/src/THDM.cpp:321`: `set_param_phys` uses atan, cos, tan and a cancellation-prone lambda1 inverse. Physical masses are cached, so getters are not an independent spectrum check.
- `2hdmc/src/THDM.cpp:383`: `set_param_phys_lam1` derives m12, delegates to physical construction, saves the pre-correction residual, then forces requested lambda1. Agreement of the final lambda1 is insufficient evidence.
- `2hdmc/src/SM.cpp:136`: v = 1/sqrt(sqrt(2) GF). Analytic-mode matching can use the existing `SM::set_GF`, with GF = 1/(sqrt(2) v_requested^2). Preserve legacy GF and record requested/backend v separately.
- `2hdmc/src/THDM.cpp:1765`: scattering blocks follow Ginzburg–Ivanov. Existing backend checks the maximum absolute eigenvalue, not the complete exported spectrum.
- `2hdmc/src/THDM.cpp:1827`: legacy threshold is 16 pi and permits equality. Campaign criterion is strict |Lambda| < 8 pi, corresponding to a0 = Lambda/(16 pi) and |Re a0| < 1/2. Do not silently change the legacy threshold.
- `appendix_2hdm_with_gao.tex:417`: independent diphoton amplitude uses fermion modifier -cot beta and charged-scalar coefficient -v^2 Z7/(2 mHp^2); widths cannot be trusted without signed coupling and loop-function checks against DecayTable.

References verified online: [Ginzburg–Ivanov](https://arxiv.org/abs/hep-ph/0508020), [2026 correction](https://link.springer.com/article/10.1140/epjc/s10052-026-15990-9). The correction changes the marginal case: positivity of the quadratic part only on quartic-zero directions is not sufficient to certify arbitrary paths to infinity. Marginal cases remain INCONCLUSIVE without the corrected path analysis.

Existing calculator tests: 56 passed. This verifies characterized legacy behavior, not general BFB, global-minimum certification, or analytic backend agreement. Heavy execution is BLOCKED because lab-health and safe-run were not found. Existing tmux sessions were inspected and left untouched.
