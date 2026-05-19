from __future__ import annotations

from pathlib import Path


def test_gemini_md_exists_and_has_guardrails() -> None:
    p = Path("/home/fabi/dihiggs/GEMINI.md")
    assert p.exists()
    text = p.read_text(encoding="utf-8")
    assert "Never run broad scans" in text
    assert "triple_ok-only" in text
    assert "2x" in text
