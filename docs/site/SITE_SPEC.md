# DiHiggs scientific site specification

## Product boundary

This is a static scientific catalog for `fbientrigo/dihiggs`, not a generic
reference manual and not a result database. Its smallest useful unit is an
evidence chain:

```text
physics convention → implementation → contract/test → benchmark data → figure → reproduction command
```

The MVP keeps Sphinx and MyST. It uses Sphinx Design for cards/tabs, small
custom templates/CSS for catalog patterns, and no client application beyond
the theme's search UI. GitHub Pages hosts the generated output. React, a
database, and an API add no requirement that static generation does not meet.

The source order for physics claims is: the canonical theoretical appendix;
`docs/PROJECT_PHYSICS_AUTHORITY.md`; then the maintained v2 evaluator
contracts. Campaign currency is separately governed by the cross-repository
current-data authority. The site must never turn a frozen or replay-only
artifact into current evidence by presenting it beside a canonical entity.

## Information architecture

The primary catalog units are scientific entities. Every detail page has a
stable ID, a status, provenance, explicit relations, source/code links, and a
reproduction affordance. Collections are views over those records, not copies
of their facts.

| Collection | Papers-with-Code analogue | DiHiggs content |
|---|---|---|
| Scientific questions | Tasks | answerable questions with an evidence path |
| Physics concepts | Methods | conventions, bases, decay/lifetime and coupling concepts |
| Evaluators | Code | maintained executable producers and orchestration entry points |
| Campaigns and data | Datasets | manifests, pilots, bounded datasets and their status |
| Results and figures | Results | observables, verification records and generated plots |
| Benchmarks | Benchmarks | frozen regression anchors with exact provenance |

The MVP home page has a short scientific purpose statement, the two task
cards, a status legend, a compact "canonical core" row, and a separate
non-canonical/replay row. It must not show a headline numerical result.

### First vertical slice

1. `/concepts/lifetime/` traces **model point → `DecayTable` → partial widths
   → total width → `ctau_mm`**. Its observable table is generated from the
   evaluator contract/selected records, and makes `width_ok` visible.
2. `/concepts/physical-h-phi-phi-coupling/` traces **appendix → Higgs-basis
   expression → 2HDMC `c` → v2 serialization → benchmark**. It shows the
   signed field and compatibility alias as different semantics.
3. `/benchmarks/h2scan-mh150-tb300000/` is a replay-only anchor page. It may
   display its frozen evidence but never labels it a current physics result.
4. `/figures/lifetime-width-closure-v2/` is the required generated v2 figure.
   Its bounded-pilot status is visible in the card, caption, gallery tile,
   alt text, and search result.

## URL taxonomy

Routes contain semantic IDs rather than filenames or source commits. A source
revision is provenance, not identity.

```text
/
/tasks/<question-id>/
/concepts/<concept-id>/
/evaluators/<evaluator-id>/
/benchmarks/<benchmark-id>/
/figures/<figure-id>/
/contracts/<contract-id>/
/verifications/<verification-id>/
/campaigns/<campaign-id>/
/status/<scientific-status>/
```

Collection routes (`/concepts/`, `/benchmarks/`, and so on) accept only
static query/filter links emitted at build time, for example
`/benchmarks/?status=replay-only`. Do not make filter state a server feature.
Deprecated IDs get Sphinx redirects only when a record supplies a stable
successor; never redirect a historical identifier to a current one without a
visible supersession note.

## Page types and interaction model

All detail pages share the following order. The schema controls their data;
templates control rendering.

1. **Identity block** — title, entity kind, one sentence, scientific-status
   badge, convention ID, and an explicit scope/disclaimer where non-canonical.
2. **Evidence chain** — a compact, accessible sequence of linked entities;
   it is a list on mobile, not a canvas graph.
3. **Scientific content** — page-type sections from `CONTENT_MAP.yaml`.
4. **Evidence and provenance** — source artifacts, checksums where supplied,
   repository/commit, producer and verification scope.
5. **Code, data, reproduce** — GitHub source links, raw artifact links, and a
   copyable command. Commands are not executed in the website.
6. **Related entities** — relationship-derived cards, not manually curated
   duplicate link lists.

Page variations are deliberately small:

| Type | Primary body | Required evidence |
|---|---|---|
| Physics concept | definition, convention, derivation/mapping | primary authority and implementation links |
| Evaluator | input/output contract and acceptance semantics | source file, contract, verification command |
| Benchmark | declared inputs, checks, replay result | frozen manifest, source/checksum, replay command |
| Figure | visual, caption, selection, limits | full figure sidecar required by `FIGURE_SCHEMA.yaml` |
| Contract | normative fields and failure rules | authoritative source, revision/checksum |
| Verification | tested scope, commands, pass/fail artifacts | report and outputs, known limits |
| Campaign | manifest and data status | campaign currency decision and data links |

Tabs are permitted only for mutually exclusive views of the same entity:
`Overview`, `Convention`, `Implementation`, `Evidence`, `Reproduce`. The
default tab is always `Overview`, and tab panels remain in the server-rendered
HTML so they are indexed and usable without JavaScript.

## Relationship model

`ENTITY_SCHEMA.yaml` is the registry. Relations use typed edges such as
`defines`, `implements`, `serializes`, `verifies`, `benchmarks`, `derives`,
`visualizes`, and `reproduces`. Rendering follows those edges in both
directions: a lifetime page links to its evaluator; the evaluator page lists
the lifetime entity it implements.

The first two paths must render exactly as follows:

```text
Lifetime:
model point → DecayTable → selected partial widths → total_width_GeV → ctau_mm

Physical coupling:
canonical appendix → Higgs-basis convention → get_coupling_hhh(1,2,2,c)
→ two_hdmc_hphiphi_* / g_physical_hphiphi_GeV → benchmark evidence
```

`g_physical_hphiphi_GeV` and `g_hH2H2_GeV` cannot share a generic "coupling"
label: the former is signed; the latter is a legacy absolute-magnitude alias.
The figure/data template likewise keeps physical lifetime separate from a
forced-response lifetime when both exist.

## Navigation and search

The desktop header contains the project wordmark, `Concepts`, `Code`,
`Benchmarks`, `Figures`, `Campaigns & data`, and `About/status`, followed by
the existing Sphinx search control. The mobile menu keeps those destinations,
then Search; it does not hide status documentation behind an icon-only control.

Use Sphinx's static search index for the MVP. Build an additional small,
generated catalog index that contributes aliases, entity IDs, kind, status,
convention ID, and source paths to search text. Search results show title,
kind, status badge, and one-line summary. They must not rank a historical or
replay-only entry above a canonical exact match merely because it has more
keywords. Status filters are static collection links, not JavaScript queries.

## Status semantics

Status is scientific scope, not a color preference. It appears on every
entity card, page hero, table row, figure caption, data download, and search
result. Text accompanies color and iconography.

| Status | Meaning | Default behavior |
|---|---|---|
| `CANONICAL` | Current authority, maintained contract, or maintained implementation for its stated scope | included in default discovery |
| `EXPERIMENTAL` | Bounded pilot, validation, or exploratory evidence with explicit limits | visible but scope callout required |
| `HISTORICAL` | Superseded material retained to explain past work | excluded from default result listings |
| `REPLAY_ONLY` | Frozen artifact permitted solely for exact replay/regression | excluded from default listings; reproduction label says replay |

The renderer computes the page's effective status as the most restrictive of
the entity and its directly displayed data/figure sources. It fails the build
when a page declared `CANONICAL` uses a non-canonical source. A current
MadGraph result must additionally have the external campaign authority's
complete approved manifest; none may be implied from the 150 GeV benchmark.

## Provenance, source links, and reproducibility

Every entity binds to an authoritative path and revision. Pages display:

- a GitHub permalink when a full commit SHA is available;
- the repository-relative source path and artifact checksum when supplied;
- the evaluator/producer identity and declared convention ID;
- a data link to the raw CSV/JSON/manifest;
- a copyable, noninteractive reproduction command and expected outputs.

The build resolves `source_commit: null` to the clean checkout SHA and records
that resolved SHA in the generated page sidecar. A dirty checkout fails a
release build. Values are extracted once from declared JSON/YAML/CSV selectors
and written to generated fragments; MyST pages only reference the fragment.
This prevents drift from hand-copying the same width, coupling, BR, or
lifetime across concept, benchmark, and figure pages.

The `Reproduce` panel has three compact rows: **build environment**, **run
command**, and **expected artifacts/checksums**. It always says whether the
command replays a frozen result, runs a bounded verification, or produces a
new canonical evaluation. A command requiring another repository links to its
pinned handoff contract; it is not presented as runnable from this repository.

## Static-generation strategy

At build time, a single deterministic metadata step shall:

1. validate `ENTITY_SCHEMA.yaml` and `FIGURE_SCHEMA.yaml`;
2. resolve source paths, source hashes, clean checkout SHA, and requested
   selectors from repository artifacts;
3. validate status propagation, convention compatibility, figure provenance,
   and referenced relation IDs;
4. generate MyST fragments, collection tables, related-entity blocks, search
   catalog text, and figure sidecars under `docs/site/_generated/`;
5. run Sphinx with MyST, Sphinx Design, the existing search mechanism, custom
   templates, and CSS; and
6. publish only the static `_build/html` tree to GitHub Pages.

The metadata step must be deterministic and reviewable: sort IDs and tables,
record input hashes, and fail rather than guess an absent convention, status,
or source. It does not calculate physics, run MadGraph, or rewrite benchmark
artifacts. A release workflow runs the metadata validation, the figure
reproduction command(s), and `sphinx-build -W` from a clean checkout.

## Explicit non-goals

- No client-side database, live job runner, accounts, comments, or API.
- No automatic interpretation of legacy notebooks or plots as scientific
  authority.
- No cross-section, acceptance, or limit result page until the owning
  repository's provenance and campaign-currency contract allow it.
- No global knowledge graph visualization; the evidence chain list is the
  accessible static MVP.

## Bounded implementation backlog

1. Add a Sphinx/MyST site root under `docs/site/` that consumes the registry
   and builds the collection and first-slice routes.
2. Implement one deterministic registry validator/generator for relations,
   provenance, status propagation, selector fragments, and catalog search text.
3. Add the site base template, status badges, evidence-chain component, and
   responsive custom CSS using Sphinx Design primitives.
4. Implement the lifetime concept and `DihiggsPointV2Evaluator` pages from
   generated fragments and the v2 contract.
5. Implement the physical h-phi-phi coupling page with signed/alias semantics
   and source permalinks.
6. Implement the replay-only `H2scan_mH150_tb300000` page and its exact replay
   panel; add a test that it cannot render as canonical.
7. Add `build_lifetime_width_closure.py`, its JSON figure sidecar, and the
   figure page; validate `FIGURE_SCHEMA.yaml` in CI.
8. Add a GitHub Pages workflow that requires a clean checkout, metadata/figure
   validation, and `sphinx-build -W` before deploy.
