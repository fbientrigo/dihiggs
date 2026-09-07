# Project physics authority

Status: **active**  
Machine contract: [`conventions/physics_conventions.yaml`](../conventions/physics_conventions.yaml)  
Decision: [issue #83](https://github.com/fbientrigo/dihiggs/issues/83)  
Adopted: 2026-09-07

## 1. Authority and scope

This document and `conventions/physics_conventions.yaml` freeze the physics
language shared by the five active repositories. The YAML is the normative
machine contract; this document gives the exact human interpretation. A
consumer must reject a mismatch between them.

The authority applies to **new production and new cross-repository handoffs**.
It does not alter the numerical content, hashes, or interpretation of frozen
historical artifacts. It also does not select a benchmark, run an evaluator,
generate events, recast events, or calculate a limit.

`fbientrigo/dihiggs` is the only owner of the project-wide contract. A copy in
another repository is a versioned cache and cannot become authoritative by
being newer, locally edited, or convenient.

## 2. Physics model and fixed basis

The project model is a **CP-conserving general 2HDM with a Type-I Yukawa
assignment**. “Type-I” fixes the Yukawa basis: `Phi2` couples to up-type quarks,
down-type quarks, and charged leptons. In that basis,

```text
tan(beta) = v2/v1,
<Phi1> = v1/sqrt(2),
<Phi2> = v2/sqrt(2).
```

`tan(beta)` is basis-dependent. It is meaningful here only because the Type-I
Yukawa basis has been fixed. A model point that declares `tan(beta)` without
the basis is incomplete.

The generic-basis quadratic convention relevant to this contract is

```text
V contains
  + m11_sq Phi1^dagger Phi1
  + m22_sq Phi2^dagger Phi2
  - (m12_sq Phi1^dagger Phi2 + h.c.).
```

The quartic convention contains

```text
V contains
  + 1/2 lambda1 (Phi1^dagger Phi1)^2
  + 1/2 lambda2 (Phi2^dagger Phi2)^2
  + ...
  + [lambda6 (Phi1^dagger Phi1)(Phi1^dagger Phi2)
     + lambda7 (Phi2^dagger Phi2)(Phi1^dagger Phi2) + h.c.].
```

Consequently, `m12_sq` is the soft `Z2`-breaking coefficient. Nonzero
`lambda6` or `lambda7` is hard `Z2` breaking. The active model must not be
described as “softly broken through lambda6/lambda7.”

## 3. State identities

| State | Canonical new label | Accepted compatibility label | PDG | New mass field |
|---|---|---|---:|---|
| SM-like CP-even scalar | `h` | `h_SM` | 25 | `m_h_GeV` |
| Additional CP-even scalar | `phi` | `H2` | 35 | `m_phi_GeV` |
| CP-odd scalar | `A` | — | 36 | `mA_GeV` |
| Charged scalar | `H+` (`H-`) | — | 37 (-37) | `mHp_GeV` |

Existing fields such as `mH`, `mH_input_GeV`, and `mH2_GeV` remain valid only
as declared compatibility fields for `phi`. New interfaces must not use the
bare label `H`: it is ambiguous between the SM-like state, the additional
CP-even state, and a mediator in existing UFO process syntax. A consumer must
fail closed rather than infer the state from a mass ordering.

## 4. SM-like Higgs mass conventions

The active project convention for new production is the decimal string

```text
m_h = "125.13" GeV
```

This adopts the PDG 2026 Higgs world-average central value, reported as
`125.13 +/- 0.11 GeV`, as the project convention. `125.20 GeV` was a prior
project convention and remains a named, replay-only convention while its
downstream consumers are migrated. It is not the PDG 2026 central value.

The decimal string is authoritative. Producers may convert it for numerical
calculation, but card writers must not reconstruct the original decimal text
from a binary float.

| Convention ID | Value [GeV] | Status | Permitted use |
|---|---:|---|---|
| `mh_125p13_pdg_2026` | `"125.13"` | active | New production and new handoffs |
| `mh_125p20_project_pre_v3` | `"125.20"` | superseded / pending migration | Exact replay of already-provenanced work only |
| `mh_125p13_high_mass_v1` | `"125.13"` | historical | High-mass v1, 150 GeV benchmark, point-v2 replay |
| `mh_125p09_lambda1_boundary_legacy` | `"125.09"` | historical | Lambda1 characterization and legacy boundary replay |
| `mh_125p0_ufo_autoresearch_legacy` | `"125.0"` | historical | Frozen UFO defaults and replay-only legacy paths |

Historical artifacts retain their original value and checksum. They are not
silently recalculated, renamed as active, or mixed across convention IDs. A
recalculation at another `m_h` is a new point with new provenance. The exact
locations still using `125.20` are enumerated in
`migration.pending_125p20_uses` in the machine contract and in the discrepancy
ledger; they are a required migration queue, not an authorization to produce.

For `dihiggs.point.v2` new production, the orchestrator reads the active mass
from the YAML and emits `convention_id`, decimal mass text, schema version,
source repository, commit, path, and SHA-256 in `mass_convention`. A requested
mass different from the active value is rejected by this new-production
interface; a historical replay must use an explicit historical interface rather
than silently override the contract. The existing point-v2 spelling `mh_GeV`
is retained as a compatibility duplicate of canonical `m_h_GeV` until the
point schema itself is versioned.

## 5. Campaign choices are not model identities

The following are common project choices, but the general model does not imply
them:

| Choice | Expression | Contract rule |
|---|---|---|
| Exact alignment | `sin(beta-alpha) = 1` | Explicit in every new manifest |
| Vanishing generic-basis quartic | `lambda7 = 0` | Explicit in every new manifest |
| Heavy-state degeneracy | `mA = mHp` | Explicit in every new manifest |

The authority supplies no default for these choices. Missing values are a
contract error, not permission to assume the values used by an older campaign.

## 6. Three distinct mass-squared quantities

All three quantities have units of `GeV^2`, but they are not aliases:

| Code name | Symbol | Definition and role |
|---|---|---|
| `M2` | `M^2` | Derived coordinate `m12_sq/(sin(beta) cos(beta))` in the fixed basis |
| `m12_sq` | `m12^2` | Off-diagonal quadratic coefficient; `M2 sin(beta) cos(beta)`; soft `Z2` breaking |
| `m22_sq` | `m22^2` | Diagonal coefficient multiplying `Phi2^dagger Phi2` in the declared generic-basis potential |

`m22_sq` is a calculator-produced model quantity and carries the calculator's
fixed-basis convention. It must travel with basis and producer provenance.
It must never be synthesized from `M2` or `m12_sq` by name similarity. The
definition of `M2` is invalid when `sin(beta) cos(beta) = 0`; such input must be
rejected.

## 7. Trilinear prescription registry

The prescription label is mandatory metadata. These three prescriptions answer
different questions and must never share an unlabeled `coupling` column.

### `physical_2HDM`

The model-derived trilinear is obtained from the validated 2HDMC model point
through `get_coupling_hhh(1,2,2,c)`. The current compatibility field
`g_hH2H2_GeV` stores an absolute magnitude. It does not preserve the signed
potential coefficient or the complete Feynman-rule convention.

Signed export is deliberately `unresolved` until issue #85 audits the
potential coefficient, `L_int = -V`, the complex 2HDMC return value, symmetry
factors, and the UFO vertex. Any sign-sensitive use must fail closed in the
meantime.

### `effective_mass_only`

The project-defined effective prescription is

```text
g_effective_hphiphi_GeV = 8 m_phi^2 / v.
```

It requires explicit `m_phi_GeV` and `v_GeV`. It is not a general 2HDM identity.

### `effective_m22_shifted`

The second project-defined effective prescription is

```text
g_effective_hphiphi_GeV = 8 (m_phi^2 - m22_sq) / v.
```

The dependence of `m22_sq` on the model point comes from the calculator; the
formula combining it with `m_phi` is project-defined. It requires
`m_phi_GeV`, `m22_sq_GeV2`, `v_GeV`, and the declared generic basis. It is not
a general 2HDM identity.

This authority intentionally does not assign a numerical default to `v_GeV`
and does not infer a UFO sign/phase mapping for either effective prescription.
Those values must be explicit and provenance-bearing at the implementing
stage; until then, sign-sensitive vertex use fails closed.

## 8. Repository ownership

| Repository | Owns | Must not claim ownership of |
|---|---|---|
| `dihiggs` | Model points, theory predicates, masses, potential parameters, widths, BRs, proper lifetime | UFO implementation, MadGraph cross sections, recast response, final limits |
| `dihiggs_ufo` | UFO implementation, parameter/vertex mapping, model-pack validation | Model-point authority, production result, recast response, final limits |
| `dihiggs_hep_cross` | MadGraph execution, cross sections, cards, production provenance | Model-point authority, recast response, final limits |
| `dihiggs_llp_recast` | Cutflows, acceptance, efficiency, response provenance | Production normalization, model-point authority, final composition |
| `dihiggs_boundary` | Composition of upstream quantities, signal yields, statistical limits, boundary products | Recalculation or silent replacement of upstream model, production, or response quantities |

The boundary between stages is scientific: cross section, branching ratio,
and acceptance-efficiency remain separate quantities until `dihiggs_boundary`
performs an explicitly provenance-linked composition.

## 9. Migration protocol

After this contract is merged, every downstream migration follows the same
sequence:

1. Read the canonical file from `fbientrigo/dihiggs` at a full 40-character
   commit SHA. Never use a floating branch name as provenance.
2. Calculate SHA-256 over the exact bytes of the canonical YAML.
3. If a repository keeps a copy, copy the YAML byte-for-byte and create
   `conventions/physics_conventions.lock.yaml` containing:

   ```yaml
   source_repository: fbientrigo/dihiggs
   source_commit: <40-character commit SHA>
   source_path: conventions/physics_conventions.yaml
   source_sha256: <64 lowercase hexadecimal characters>
   schema_version: physics_conventions_v3
   copied_utc: <UTC timestamp>
   ```

4. Verify the checksum before loading the cache. A local edit requires a new
   upstream canonical revision; it must not create a downstream authority.
5. Every new model-point handoff records `physics_conventions_schema`,
   `physics_conventions_commit`, and `physics_conventions_sha256`.
6. Unknown schema, missing provenance, checksum mismatch, or missing required
   physics input is a hard rejection.

The canonical commit cannot be embedded in the canonical file itself without
creating a self-reference. It is recorded by each downstream sidecar and each
handoff after the canonical commit exists.

## 10. Machine-to-human field crosswalk

| YAML path | Human section |
|---|---|
| `schema_version`, `authority` | Section 1 |
| `constants` | Machine constants retained from the previous contract |
| `mass_conventions` | Section 4 |
| `states` | Section 3 |
| `model` | Section 2 |
| `campaign_choices` | Section 5 |
| `parameters` | Section 6 |
| `trilinear_prescriptions` | Section 7 |
| `stage_ownership` | Section 8 |
| `migration` | Section 9 |
| `failure_policy` | Section 11 |

## 11. Fail-closed rule

Reject rather than guess when any of the following occurs:

- a bare `H` is used to identify a state;
- the mass convention or Type-I basis is absent or unknown;
- a campaign choice is defaulted rather than declared;
- `M2`, `m12_sq`, and `m22_sq` are conflated;
- the trilinear prescription or one of its required inputs is absent;
- an unresolved sign/phase mapping is used in a sign-sensitive calculation;
- a downstream copy lacks canonical commit/checksum provenance or mismatches it.

## 12. Scientific references

- F. Takahashi et al. (Particle Data Group), *Review of Particle Physics*
  (2026), Higgs-boson listing: the cited average is `125.13 +/- 0.11 GeV`.
- D. Eriksson, J. Rathsman, and O. Stal, *2HDMC — Two-Higgs-Doublet Model
  Calculator*, arXiv:0902.0851.
- G. C. Branco et al., *Theory and phenomenology of two-Higgs-doublet
  models*, arXiv:1106.0034.

The PDG listing supplies the adopted `"125.13"` mass convention. The two
effective trilinear formulas remain project decisions.
