from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict

def write_status_json(outdir: str | Path, snapshot: dict[str, Any]) -> Path:
    """Write campaign status to campaign_status.json."""
    target = Path(outdir) / "campaign_status.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(snapshot, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    return target


def write_status_html(outdir: str | Path, snapshot: dict[str, Any]) -> Path:
    """Write a standalone HTML dashboard to campaign_status.html."""
    target = Path(outdir) / "campaign_status.html"
    target.parent.mkdir(parents=True, exist_ok=True)

    # Extract required fields with safe defaults
    state = snapshot.get("campaign_state", "UNKNOWN")
    active_round = snapshot.get("round_index", "N/A")
    last_success = snapshot.get("last_successful_round", "N/A")
    
    metrics = snapshot.get("rolling_metrics", {})
    m_yield = metrics.get("yield", 0.0)
    m_coverage = metrics.get("coverage", 0.0)
    m_diversity = metrics.get("diversity", 0.0)
    m_composite = metrics.get("composite", 0.0)
    
    alerts = snapshot.get("active_alerts", [])
    last_scaling = snapshot.get("last_scaling_action", "None")
    manifest_compat = snapshot.get("manifest_compatibility", "Unknown")
    file_timestamps = snapshot.get("file_timestamps", {})

    # Build HTML
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Campaign Status Dashboard</title>
    <style>
        body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; padding: 20px; line-height: 1.6; color: #333; background: #f9f9f9; }}
        .container {{ max-width: 1000px; margin: 0 auto; background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        h1, h2 {{ color: #2c3e50; border-bottom: 1px solid #eee; padding-bottom: 10px; }}
        .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 20px; margin-bottom: 20px; }}
        .card {{ background: #f8f9fa; padding: 15px; border-radius: 6px; border: 1px solid #e9ecef; }}
        .card h3 {{ margin-top: 0; font-size: 0.9em; color: #6c757d; text-transform: uppercase; letter-spacing: 1px; }}
        .card .value {{ font-size: 1.5em; font-weight: bold; color: #2c3e50; }}
        .status-RUNNING {{ color: #28a745; }}
        .status-FAILED {{ color: #dc3545; }}
        .status-COMPLETED {{ color: #007bff; }}
        .alerts {{ background: #fff3cd; border: 1px solid #ffeeba; color: #856404; padding: 15px; border-radius: 6px; }}
        .alert-item {{ margin-bottom: 10px; padding: 10px; background: white; border-left: 4px solid #dc3545; }}
        table {{ width: 100%; border-collapse: collapse; margin-top: 10px; }}
        th, td {{ padding: 8px; text-align: left; border-bottom: 1px solid #ddd; }}
        th {{ background-color: #f2f2f2; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Campaign Dashboard</h1>
        
        <div class="grid">
            <div class="card">
                <h3>Campaign State</h3>
                <div class="value status-{state}">{state}</div>
            </div>
            <div class="card">
                <h3>Active Round</h3>
                <div class="value">{active_round}</div>
            </div>
            <div class="card">
                <h3>Last Successful Round</h3>
                <div class="value">{last_success}</div>
            </div>
            <div class="card">
                <h3>Manifest Compatibility</h3>
                <div class="value">{manifest_compat}</div>
            </div>
        </div>

        <h2>Rolling Metrics</h2>
        <div class="grid">
            <div class="card">
                <h3>Yield</h3>
                <div class="value">{m_yield:.4f}</div>
            </div>
            <div class="card">
                <h3>Coverage</h3>
                <div class="value">{m_coverage:.4f}</div>
            </div>
            <div class="card">
                <h3>Diversity</h3>
                <div class="value">{m_diversity:.4f}</div>
            </div>
            <div class="card">
                <h3>Composite</h3>
                <div class="value">{m_composite:.4f}</div>
            </div>
        </div>

        <h2>Active Alerts</h2>
        <div class="alerts">
            {_render_alerts(alerts)}
        </div>

        <div class="grid" style="margin-top: 20px;">
            <div class="card">
                <h3>Last Scaling Action</h3>
                <div><pre style="margin:0; white-space: pre-wrap;">{json.dumps(last_scaling, indent=2)}</pre></div>
            </div>
        </div>

        <h2>File Timestamps (Latest Raw Artifacts)</h2>
        <table>
            <tr><th>File</th><th>Timestamp</th></tr>
            {_render_timestamps(file_timestamps)}
        </table>
    </div>
</body>
</html>"""
    target.write_text(html, encoding="utf-8")
    return target


def _render_alerts(alerts: list[Any]) -> str:
    if not alerts:
        return "<p>No active alerts.</p>"
    items = []
    for alert in alerts:
        alert_dict = alert if isinstance(alert, dict) else {"raw": str(alert)}
        type_val = alert_dict.get("alert_type", "Unknown")
        sev_val = alert_dict.get("severity", "info")
        msg_val = alert_dict.get("message", json.dumps(alert_dict))
        items.append(f'<div class="alert-item"><strong>{type_val} ({sev_val})</strong>: {msg_val}</div>')
    return "".join(items)


def _render_timestamps(timestamps: dict[str, Any]) -> str:
    if not timestamps:
        return "<tr><td colspan='2'>No timestamp data available.</td></tr>"
    rows = []
    for fname, ts in timestamps.items():
        rows.append(f"<tr><td>{fname}</td><td>{ts}</td></tr>")
    return "".join(rows)


def print_cli_status(outdir: str | Path) -> None:
    """Read campaign_status.json and print a compact summary to stdout."""
    status_file = Path(outdir) / "campaign_status.json"
    if not status_file.exists():
        print(f"Error: No campaign_status.json found in {outdir}", file=sys.stderr)
        sys.exit(1)
        
    try:
        data = json.loads(status_file.read_text(encoding="utf-8"))
    except Exception as e:
        print(f"Error reading {status_file}: {e}", file=sys.stderr)
        sys.exit(1)

    state = data.get("campaign_state", "UNKNOWN")
    round_idx = data.get("round_index", "N/A")
    metrics = data.get("rolling_metrics", {})
    alerts = data.get("active_alerts", [])
    
    print(f"Campaign Status: {state}")
    print(f"Active Round:    {round_idx}")
    print(f"Metrics:         Yield={metrics.get('yield', 0.0):.4f}, Coverage={metrics.get('coverage', 0.0):.4f}, "
          f"Diversity={metrics.get('diversity', 0.0):.4f}, Composite={metrics.get('composite', 0.0):.4f}")
    if alerts:
        print(f"Alerts:          {len(alerts)} active alerts")
        for a in alerts:
            a_dict = a if isinstance(a, dict) else {}
            print(f"  - [{a_dict.get('severity', 'info').upper()}] {a_dict.get('alert_type', 'Unknown')}: {a_dict.get('message', '')}")
    else:
        print("Alerts:          None")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Autoresearch Campaign Status Dashboard CLI")
    parser.add_argument("outdir", help="Path to the campaign output directory")
    args = parser.parse_args()
    print_cli_status(args.outdir)
