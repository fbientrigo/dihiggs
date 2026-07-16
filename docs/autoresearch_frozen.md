# Autoresearch frozen status

`autoresearch/` is a preserved historical experiment, not part of the canonical
DiHiggs core. Its fixed/15 CSV path loses long-lived widths and its `width > 0`
lifetime filter can discard affected points. Existing files, campaigns, and
entry points remain available for historical replay, but receive no repairs or
new adapters here.

New model evaluation must use `dihiggs/app/Lambda1EvaluatorV2`, whose source is
`dihiggs/src/Lambda1EvaluatorV2.cpp`. Canonical code under `dihiggs/` must not
import `autoresearch`; the test suite enforces that one-way boundary.
