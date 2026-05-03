"""
Unified TUI for autoresearch physics exploration.

Run with:
    python -m autoresearch.tui
    python -m autoresearch.tui --config path/to/config.json
    python -m autoresearch.tui --campaign path/to/campaign/
"""

from .app import AutoresearchApp, run_tui, main

__all__ = ["AutoresearchApp", "run_tui", "main"]
