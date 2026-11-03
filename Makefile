# autoresearch: tests + real end-to-end smoke (needs ../dihiggs_boundary built: run `make setup` in the workspace root)
PY ?= $(if $(wildcard ../.venv/bin/python),../.venv/bin/python,python3)

.PHONY: autoresearch-smoke
autoresearch-smoke:
	$(PY) -m pytest -q autoresearch/tests
	$(PY) -m autoresearch.smoke
