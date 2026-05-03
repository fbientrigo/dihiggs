import json
import tempfile
from pathlib import Path

from autoresearch.harness.status_dashboard import write_status_json, write_status_html, print_cli_status

def test_write_status_json():
    with tempfile.TemporaryDirectory() as temp_dir:
        outdir = Path(temp_dir)
        snapshot = {
            "campaign_state": "RUNNING",
            "round_index": 5,
            "last_successful_round": 4,
            "rolling_metrics": {
                "yield": 0.4,
                "coverage": 0.5,
                "diversity": 1.2,
                "composite": 0.71
            },
            "active_alerts": [],
            "last_scaling_action": {"action": "scale_up", "threads": 8},
            "manifest_compatibility": "Compatible",
            "file_timestamps": {
                "scan_meta.json": "2026-04-06T12:00:00Z"
            }
        }
        
        path = write_status_json(outdir, snapshot)
        assert path.exists()
        assert path.name == "campaign_status.json"
        
        data = json.loads(path.read_text(encoding="utf-8"))
        assert data["campaign_state"] == "RUNNING"
        assert data["round_index"] == 5
        assert "rolling_metrics" in data
        assert data["rolling_metrics"]["yield"] == 0.4
        assert data["last_successful_round"] == 4

def test_write_status_html():
    with tempfile.TemporaryDirectory() as temp_dir:
        outdir = Path(temp_dir)
        snapshot = {
            "campaign_state": "FAILED",
            "active_alerts": [
                {"alert_type": "round_failed", "severity": "critical", "message": "Round 10 failed."}
            ]
        }
        
        path = write_status_html(outdir, snapshot)
        assert path.exists()
        assert path.name == "campaign_status.html"
        
        html_content = path.read_text(encoding="utf-8")
        assert "FAILED" in html_content
        assert "round_failed" in html_content
        assert "critical" in html_content
        assert "Round 10 failed." in html_content

def test_print_cli_status(capsys):
    with tempfile.TemporaryDirectory() as temp_dir:
        outdir = Path(temp_dir)
        snapshot = {
            "campaign_state": "CONVERGED",
            "round_index": 12,
            "rolling_metrics": {
                "yield": 0.8,
                "coverage": 1.0,
                "diversity": 2.0,
                "composite": 1.5
            },
            "active_alerts": [
                {"alert_type": "low_yield", "severity": "warning", "message": "Yield dropping"}
            ]
        }
        
        write_status_json(outdir, snapshot)
        
        # Capture stdout
        print_cli_status(outdir)
        captured = capsys.readouterr()
        
        assert "CONVERGED" in captured.out
        assert "12" in captured.out
        assert "Yield=0.8000" in captured.out
        assert "low_yield" in captured.out
        assert "warning" in captured.out.lower()
