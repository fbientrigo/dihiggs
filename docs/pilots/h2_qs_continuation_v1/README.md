# Task 2 pilot — deterministic Q,S continuation (150 → 175 → 200 GeV)

Full evidence: [`continuation_report.json`](continuation_report.json).

## Question

Starting from the frozen, theory-valid `H2scan_mH150_tb300000` anchor
(m_H2=150 GeV, m_A=m_Hp=450 GeV, tan_beta=3e5, sin(β−α)=1, λ6=1e-10, λ7=0,
Type I), does holding

```
Q = (m_H2^2 - M2) * tan_beta^2
S = m_A^2 - m_H2^2
```

fixed while raising m_H2 keep the point on the theory-valid 2HDM branch? All
other protected coordinates (tan_beta, λ6, λ7, sin(β−α), Yukawa type) are
held fixed throughout; only m_A=m_Hp and M2 are re-derived from Q,S.

## Method

`dihiggs_boundary/python/dhb/continuation.py` computes each child proposal
as pure arithmetic (no 2HDMC, no theory evaluation, no search). Each
proposal was submitted **one at a time** through the frozen `dihiggs`
proposal-batch socket (`docs/contracts/adaptive_proposal_batch_v1.yaml`,
`dihiggs/app/orchestrator/proposal_batch.py`) via
`dihiggs_boundary/scripts/qs_continuation_step.py`, which contains no mass
loop, no retry, and no rescue. Only `DihiggsPointV2Evaluator` decided
validity.

## Result

| m_H2 (GeV) | role | m_A=m_Hp (GeV) | theory_ok_v1 | rejection |
|---|---|---|---|---|
| 150 | anchor | 450.000 | 1 | accepted |
| 175 | direct Q,S continuation | 458.939 | 1 | accepted |
| 200 | direct Q,S continuation | 469.042 | 1 | accepted |

Both 175 and 200 GeV children were theory-valid on the first and only
proposal submitted (no repair, no neighbor search). Q and S were preserved
within their documented Float64 ULP-scale rounding budgets. λ2–λ5 were
unchanged to machine precision across all three points; λ1 (which absorbs
most of Q in the alignment limit) shifted by only ~1e-5 relative per 25 GeV
step.

**Verdict: `DIRECT_QS_CONTINUATION_VALID_TO_200_GEV`.**

## Next task justified by this evidence

Per the Task 2 handoff's decision rule for a Gate‑4 PASS: **not** Task 3
(structured local Q,S rescue) — direct continuation needed no rescue. The
justified next step is a bounded extension to the next deliberate mass
milestone, still without stochastic search (Monte Carlo remains a later
task).
