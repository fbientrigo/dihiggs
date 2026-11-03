"""#19: the TUI shows budget remaining, the action menu and coverage from journal.db (Textual pilot, headless)."""
import asyncio
import json

from autoresearch.harness.journal import Journal
from autoresearch.tui.app import AutoresearchApp
from autoresearch.tui.services.monitor import journal_text, journal_view

A, B = "survivor_v2|tb=100.0|r", "survivor_v2|tb=10000000.0|r"


def campaign(tmp_path):
    (tmp_path / "action_menu.json").write_text(json.dumps({"digest": "x", "actions": [{"key": A}, {"key": B}]}))
    j = Journal(tmp_path / "journal.db")
    j.append("limits:x", "limits", {"cpu_core_s": 3600, "disk_mb": 1024})
    j.reserve("a0000", {"cpu_core_s": 100, "disk_mb": 10}, {}, {"action": A})
    j.close("a0000", {"cpu_core_s": 8, "disk_mb": 1, "wall_s": 4}, {"action": A, "status": "success", "fresh_rows": 32})
    return tmp_path


def test_journal_view_remaining_actions_coverage(tmp_path):
    v = journal_view(campaign(tmp_path))
    assert v["remaining"] == {"cpu_core_s": 3592, "disk_mb": 1023}
    assert v["covered"] == 1 and v["menu_size"] == 2 and v["actions"][B]["runs"] == 0 and v["actions"][A]["fresh_rows"] == 32


def test_overview_renders_journal(tmp_path):
    out = campaign(tmp_path)

    async def run():
        app = AutoresearchApp()
        async with app.run_test(size=(140, 50)) as pilot:
            app.load_campaign(out)
            await pilot.pause()
            return str(app.query_one("#journal-info").render())

    text = asyncio.run(run())
    assert text == journal_text(journal_view(out))
    assert "cpu_core_s: 3592/3600" in text and "Coverage: 1/2" in text and f"{A}  runs=1 fresh_rows=32" in text
