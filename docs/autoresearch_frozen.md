# Autoresearch frozen status

`autoresearch/` is a preserved historical experiment, not part of the canonical
DiHiggs core. Its fixed/15 CSV path loses long-lived widths and its `width > 0`
lifetime filter can discard affected points. Existing files, campaigns, and
entry points remain available for historical replay, but receive no repairs or
new adapters here.

The canonical producers are `dihiggs/app/Lambda1EvaluatorV2` from
`dihiggs/src/Lambda1EvaluatorV2.cpp` and
`dihiggs/app/DihiggsPointV2Evaluator` from
`dihiggs/src/DihiggsPointV2Evaluator.cpp`. Autoresearch is frozen; no canonical
code imports it. New LLP lifetime evaluation must use one of these v2
producers. The test suite enforces that one-way boundary.
