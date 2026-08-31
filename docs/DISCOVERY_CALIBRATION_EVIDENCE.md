# Discovery calibration evidence — READY_FOR_DISCOVERY_SEARCH

Status: master-side deterministic calibration, completed before the multi-agent
campaign. Every number below was read from the frozen
`dihiggs/app/DihiggsPointV2Evaluator` binary, invoked with the same argv the
frozen substrate uses. No physics code was modified.

Repository: `/home/fabian/atlas_dihiggs/dihiggs`
Branch: `feat/h2-highmass-continuation-mc`
HEAD: `ea2681069d2edbe9b059e7b2e17f3169d5a4653d`
Evaluator sha256: `a6153efd4ddc0712094163e85306214e5e6ec4cb32b0872c5fc89afef529cce2`

## 1. Readiness model

The previous mission required `READY_FOR_MULTIAGENT_SEARCH` and stopped when it
was absent. That gate has been replaced, by explicit scientific-owner
instruction, with **`READY_FOR_DISCOVERY_SEARCH`**.

The distinction is deliberate:

| gate | requires | status |
|---|---|---|
| `READY_FOR_MULTIAGENT_SEARCH` | authoritative mixed/photonic thresholds, active-125.20 fixtures, justified lifetime envelope | still NOT met |
| `READY_FOR_DISCOVERY_SEARCH` | frozen evaluator intact, gates intact, provenance intact, non-empty search space demonstrated | **MET** |

The unresolved classification question is treated as a *search question*, not a
blocker. Family labels produced during discovery are explicitly PROVISIONAL
(see `search_discovery/families.py`) and are never written into the frozen
contract or the frozen fail-closed archive.

## 2. What remains immutable

Untouched, and imported rather than reimplemented by the discovery layer:
evaluator semantics and the physics binary; Yukawa initialization; the
construction / numerical / positivity / unitarity / perturbativity / theory /
width gates (`search_substrate.evaluator._validity`); the ctau calculation and
`HBAR_C_GEV_MM = 1.973269804e-13`; parameter normalization semantics and
candidate identity hashing; provenance stamping; the append-only ledger;
deterministic execution (`OMP_NUM_THREADS=1`, `LC_ALL=C`).

`tests/test_discovery_layer.py` pins this: candidate identity, the canonical
physical JSON, the derived parameters and the **argv handed to the physics
binary** are asserted byte-identical between the frozen path and the discovery
path for a point inside both envelopes. The 8 immutable substrate tests still
pass unchanged (21 passed total).

## 3. Why the previous probe found 0/200 valid

The frozen v1 contract requires the derived quantity `X = lambda6 * tan_beta` to
lie in `[200, 100000]`. Direct measurement across three decades of `tan_beta`
gives, exactly:

```
lambda1_reconstructed = m_h^2/v^2 - 1.5 * X        (coefficient 1.5000)
```

| tan_beta | X | lambda1 | positivity |
|---|---|---|---|
| 1e4 | 0.15 | +0.033559 | 1 |
| 1e4 | 0.20 | -0.041441 | 0 |
| 1e5 | 0.15 | +0.033558 | 1 |
| 1e6 | 0.15 | +0.033516 | 1 |
| 1e6 | 1.00 | -1.241484 | 0 |

Positivity requires `lambda1 > 0`, i.e. `X < m_h^2 / (1.5 v^2) ~ 0.172`. The
frozen envelope's floor of `X = 200` therefore forces `lambda1 <= -300` at every
point it permits. **The 0/200 result was an envelope artifact, not evidence that
the model has no valid region.** The earlier report was right to label it
evidence rather than proof; it is now superseded.

## 4. Q, not M2, is the physical dial

At `tan_beta = 1e6`, changing `M2` by 1e-4 GeV^2 out of 40000 flips unitarity.
That is not physics; it is the coordinate. The gates respond to the invariant

```
Q = (mH^2 - M2) * tan_beta^2
```

which is exactly the quantity the frozen family continuation already preserves.
Scanning Q directly at fixed `X = 0.10`, `mA = 500`, `mH = 200`:

| Q | ctau_mm | br_bb | br_gammagamma | theory_ok_v1 |
|---|---|---|---|---|
| 0 | 33.97 | 0.6754 | 0.00051 | 1 |
| 1e4 | 33.97 | 0.6755 | 0.00051 | 1 |
| 1e5 | 33.98 | 0.6755 | 0.00048 | 1 |
| 3e5 | 33.98 | 0.6755 | 0.00041 | 0 (perturbativity) |
| 1e6 | 33.99 | 0.6757 | 0.00020 | 0 (unitarity) |

Results depend on Q, not on `tan_beta` separately. `ctau` scales as
`tan_beta^2` at fixed Q. The discovery layer therefore searches in
`(mH, mA, tan_beta, X, Q)` and derives `M2` and `lambda6`.

## 5. Demonstrated non-empty target region

A theory-valid point **inside the required 500-1000 mm anchor interval**:

```
mH = 200 GeV, mA = mHp = 500 GeV, M2 = 40000 GeV^2 (Q = 0),
tan_beta = 4.54e6, lambda6 = 2.2026e-8  (X = 0.10)
-> ctau = 699.99 mm, theory_ok_v1 = 1
   br_bb 0.6753, br_gg 0.2221, br_tautau 0.0713, br_gammagamma 0.00051
```

`tan_beta` tunes the lifetime across the window: 4.0e6 -> 543.7 mm,
4.54e6 -> 700.0 mm, 5.0e6 -> 849.4 mm, 6.0e6 -> 1223.5 mm.

A 16-point Sobol' batch in envelope E0 returned **7/16 theory-valid**, including
a valid point at `ctau = 808.9 mm`. Compare 0/200 under the frozen v1 envelope.

This retires frozen-contract scientific blocker #3 (no justified lifetime
envelope): the 500-1000 mm interval is now demonstrably reachable.

## 6. The two basins, and the open question

The evaluator produces a sharply bimodal branching structure with essentially
nothing between the modes:

- fermion-dominated: `br_bb + br_gg + br_tautau ~ 0.97`
- photon-dominated: `br_gammagamma + br_Zgamma ~ 0.685 + 0.315 ~ 1.00`

The photon modes are driven by the charged-Higgs loop, whose coupling grows with
`|mH^2 - M2|` and is *not* cot(beta)-suppressed, whereas every fermionic mode is.
This is why the branching pattern is invariant under `tan_beta` at fixed Q.

**Open question for the campaign:** in master calibration the photon-dominated
pattern appeared only where unitarity and perturbativity FAIL — photon dominance
needs large `mH^2 - M2`, while the theory gates cap `|Q| = |mH^2 - M2| tan_beta^2`.
Whether a theory-valid photonic point exists anywhere is *not settled* and is
assigned to Strategy B and to worker C2. A well-evidenced negative is an
acceptable and valuable campaign outcome; an unevidenced positive is not.

## 7. Three-tier bounds model

1. **`GLOBAL_PHYSICS_BOUNDS`** (`search_discovery/bounds.py`) — immutable
   physics/safety limits, deliberately permissive. They span the entire frozen
   v1 contract range, so the discovery layer can reach any point the old
   contract allowed and can re-measure the positivity and unitarity walls rather
   than inherit them. These bounds never decide validity; the frozen gates do.
2. **`Envelope`** (`search_discovery/envelopes.py`) — master-controlled adaptive
   search region, validated on construction to lie inside tier 1. `E0` is the
   calibration-wave envelope and is deliberately WIDER than the known-valid
   region, spanning both theory boundaries from both sides so the campaign
   measures them. `E1_mixed_exploit` and `E2_photonic_hunt` are exploitation
   envelopes.
3. **worker-local neighbourhoods** — mutation radii and CMA-ES step sizes,
   always clipped to the active tier-2 envelope.

## 8. Caveat carried forward

The whole substrate, and now the discovery layer, remain **untracked in git**.
`repository_identity()` stamps `repository_dirty: true` on every ledger record,
and the frozen contract hash references an uncommitted file. This does not block
discovery, but it must be resolved before any result here is treated as a
provenance-bearing scientific artifact.

## 9. Scope

This campaign selects physics benchmark candidates. It does NOT establish
MadGraph production rates, does NOT calculate detector acceptance, does NOT
validate a recast, and does NOT produce event-yield tables.
