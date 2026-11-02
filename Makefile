# autoresearch: tests + real end-to-end smoke (needs ../dihiggs_boundary/build and ../.venv)
PY ?= ../.venv/bin/python

.PHONY: autoresearch-smoke
autoresearch-smoke:
	$(PY) -m pytest -q autoresearch/tests
	$(PY) -m autoresearch.smoke
