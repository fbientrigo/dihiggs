# Autoresearch results plotting

Manual refresh command (intentionally manual):

`/home/fabi/higgs_env_py312/bin/python /home/fabi/dihiggs/autoresearch/results/build_ctau_evolution.py`

This regenerates:
- `logs/ctau_discoveries.csv`
- `logs/ctau_frontier_events.csv`
- `ctau_evolution.png`

Notes:
- Bubble size in the plot reflects stage attempted quota used at discovery.
- Data source is run artifacts under `scripts/out/gemini_harness/boulder_explore_ladder/` plus validated baseline from `boulder.json`.
