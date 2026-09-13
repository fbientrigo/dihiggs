# DiHiggs design specification

## Visual direction

The site borrows Papers with Code's research-catalog rhythm: a dense but calm
header, discoverable topic pages, compact result cards, clear code/data links,
and tables that reward comparison. It is not a clone: DiHiggs puts convention,
provenance, and scientific status ahead of popularity metrics.

Use the existing Sphinx theme as the foundation, Sphinx Design cards/tabs for
semantic components, and a small `docs/site/_static/site.css`. Custom
templates should be limited to the header, entity hero, status badge, evidence
chain, card, and provenance/reproduce panels. No component framework is
needed.

## Layout

Desktop pages use a centered content rail (maximum 1200px) with a primary
reading column (about 720px) and an optional 280px evidence rail. The evidence
rail is sticky only when it fits in the viewport; it contains status,
convention, source/data/code links, and the reproduction button. A page must
remain complete when that rail moves below content.

```text
┌──────────────────────────────────────────────────────────────────────┐
│ wordmark  Concepts  Code  Benchmarks  Figures  Campaigns  [Search]  │
├──────────────────────────────────────────────────────────────────────┤
│ breadcrumb / entity kind                         [STATUS]             │
│ Title and one-sentence scientific scope                                │
│ convention • producer • commit                                         │
├───────────────────────────────────────┬──────────────────────────────┤
│ Evidence chain / scientific content   │ Status + scope                │
│ tables, figure, tabs                  │ Code • Data • Source           │
│ related entities                      │ Reproduce                      │
└───────────────────────────────────────┴──────────────────────────────┘
```

The evidence chain uses short linked pills and arrows on desktop. At narrow
width it becomes a numbered list, preserving order and target links.

## Header and navigation

Header height is 56–64px, with a visible wordmark and text labels for all
primary destinations. Search is a text control with a keyboard-focus style;
it opens the Sphinx search page or overlay only if the selected theme already
supports it. The header is sticky after a small scroll offset and has a solid
background, not translucent text over scientific figures.

Breadcrumbs show collection → entity. The current location is textual, not
only distinguished by color. The mobile menu has the same route order and a
separate visible Search action.

## Cards, tables, and figures

Cards use a 1px neutral border, 8px radius, white/off-white surface, 16–20px
padding, and a modest shadow only on hover/focus. Cards always contain: kind,
title, status, one-line scope, and at least one useful destination (source,
data, code, or reproduce). Never use a card grid to hide scientific caveats.

Benchmark/result tables are compact, left-align labels, right-align numerical
columns, keep units in headers, and include a dedicated `Status` column.
Tables use row striping only lightly, horizontal scrolling on small screens,
and a source/manifest link directly beneath the table. No value appears in a
table unless the registry binds it to a source selector.

Figures have a clear title, `figcaption`, meaningful alt text, visible status,
convention ID, producer, data link, and reproduce link. Gallery tiles use a
fixed visual frame but never crop axes, legends, status banners, or captions.
The gallery defaults to canonical figures; a status filter reveals other
material with its warning intact.

## Badges and status callouts

Badges have icon, text, and a tooltip/expanded definition; color is redundant.
Use one badge per entity identity, not a rainbow of granular states.

| Status | Badge treatment | Callout copy pattern |
|---|---|---|
| CANONICAL | deep teal, check icon | "Canonical for the stated scope." |
| EXPERIMENTAL | indigo, flask icon | "Bounded evidence; see scope and limits." |
| HISTORICAL | muted amber, archive icon | "Historical; not current evidence." |
| REPLAY_ONLY | charcoal/amber, replay icon | "Replay-only; not for new-production claims." |

Every non-canonical detail page starts with a full-width callout immediately
below the hero. Status colors meet WCAG AA against the surface; status text is
always present in exports/print styles.

## Typography and spacing

Use a system sans-serif UI stack for interface text and a readable serif or
the theme's body face for prose only if it is already available without a
network font. Mathematical expressions use MyST/MathJax conventions; source
code uses the theme monospace face. Keep body text at 16–18px with 1.5–1.65
line height. Code and provenance strings may reduce to 14px but remain
copyable and wrap safely.

Use a four-point spacing scale: 4, 8, 16, 24, 32, 48px. Typical entity hero
spacing is 32px top/bottom; cards use 16px; sections use 32px; dense table
cells use 8px vertical/12px horizontal. Titles use a restrained type scale
rather than oversized marketing typography.

## Tabs and interaction

Use Sphinx Design tab sets only for the same page's alternate view; URL
anchors must target each panel. Tabs have visible active and keyboard-focus
states. The first panel is meaningful without JavaScript, and all panels are
in the static DOM for search and printing. The only permitted custom vanilla
JavaScript is a short copy-to-clipboard enhancement for reproduction commands;
the command remains selectable text and links still work without it.

## Responsive and accessible behavior

At ≤960px, move the evidence rail below the hero and collapse card grids to
two columns. At ≤640px, use one column, convert the evidence chain to the
ordered list, retain a 44px minimum target size, and allow tables to scroll in
a labelled wrapper. Do not hide source, data, or status information on mobile.

All controls work with keyboard navigation. Tables retain header associations;
figures need informative alt text plus captions; badge contrast is tested;
focus rings are never removed. The print stylesheet expands tab content,
prints status and provenance, and places URLs beside source/reproduction links.
