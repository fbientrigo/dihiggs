# Artifacts and schema semantics

The repository has two canonical row-preserving schemas:

- `dihiggs.lambda1.v2`;
- `dihiggs.point.v2`.

Both preserve attempted rows rather than dropping construction failures or
theory-rejected points.

## Core status fields

Interpret these fields separately:

- `construction_ok`: 2HDMC parameter construction succeeded;
- `positivity`, `unitarity`, `perturbativity`: theory predicates;
- `theory_ok_v1`: current theory-only conjunction of those predicates;
- `width_ok`: lifetime-producing width is numerically usable;
- `ctau_mm`: proper decay length in millimetres when `width_ok` is true.

A theory-valid row is not an experimental acceptance decision.

## Mass-squared quantities

The canonical project language is:

```text
M2      = m12_sq / (sin(beta) cos(beta))
m12_sq  = M2 * sin(beta) cos(beta)
m22_sq  = diagonal generic-basis quadratic coefficient
```

All have units of GeV^2, but they are not aliases.

## Higgs mass provenance

New point-v2 production resolves the active Higgs-mass convention from the
project physics authority and records provenance in `run_manifest.json`.
The active convention is `mh_125p13_pdg_2026`, decimal text `"125.13"`.

Historical artifacts keep their original convention ID and value.

## Trilinear fields

For the physical `h phi phi` coupling, the signed observable is
`g_physical_hphiphi_GeV`. The compatibility field `g_hH2H2_GeV` is the
absolute magnitude and must not be used when the sign matters.

See [Project physics authority](PROJECT_PHYSICS_AUTHORITY.md) for the complete
serialization and sign convention.
