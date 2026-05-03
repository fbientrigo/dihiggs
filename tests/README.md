# Test layout after scripts cleanup

Active Colab analyzer tests (primary)
- `tests/test_validate_colab.py`
- `tests/test_show_width_rel_stats.py`
- `tests/test_group_width_errors_by_config.py`

Legacy drift/lambda1 tests (optional)
- Moved under `tests/legacy/` and renamed with `legacy_test_*.py` to avoid confusion with active analyzer workflow.
- They validate legacy code now located in `legacy/colab_compare/`.
- Run manually only when auditing legacy drift helpers/reports.

Manual legacy test commands
```bash
python -m pytest tests/legacy/legacy_test_stat_drift_helpers.py
python -m pytest tests/legacy/legacy_test_export_lambda1_compact_report.py
python -m pytest tests/legacy/legacy_test_adaptive_gate_and_export.py
```
