"""
Monitoring service - adapts live_monitor.py for TUI use.
"""
import json
from pathlib import Path
from typing import Dict, Optional, Any
from dataclasses import dataclass
from datetime import datetime


@dataclass
class CampaignStatus:
    """Campaign status snapshot for TUI display."""
    
    state: str = "UNKNOWN"
    round_index: int = -1
    total_events: int = 0
    total_run_dirs: int = 0
    elapsed_hours: float = 0.0
    max_duration_hours: Optional[float] = None
    
    # Metrics
    yield_val: float = 0.0
    coverage: float = 0.0
    diversity: float = 0.0
    composite: float = 0.0
    
    # Health
    is_stale: bool = False
    last_activity: Optional[datetime] = None
    consecutive_empty_rounds: int = 0
    
    # Errors
    error: Optional[str] = None
    
    @property
    def progress_pct(self) -> Optional[float]:
        """Calculate progress percentage if time budget is set."""
        if self.max_duration_hours and self.max_duration_hours > 0:
            return min(100.0, (self.elapsed_hours / self.max_duration_hours) * 100)
        return None
    
    @property
    def state_emoji(self) -> str:
        """Get emoji for current state."""
        return {
            "INIT": "⏳",
            "RUNNING": "🔄",
            "COMPLETED": "✅",
            "CONVERGED": "🎯",
            "FAILED": "❌",
            "PREFLIGHT_BLOCKED": "🚫",
            "UNKNOWN": "❓",
        }.get(self.state, "❓")


class MonitoringService:
    """Service for monitoring campaign state via file polling."""
    
    def __init__(self, outdir: Optional[Path] = None):
        self.outdir = outdir
    
    def set_outdir(self, outdir: Path):
        """Set the campaign output directory to monitor."""
        self.outdir = Path(outdir)
    
    def get_status(self) -> CampaignStatus:
        """Get current campaign status by reading state files."""
        status = CampaignStatus()
        
        if self.outdir is None or not self.outdir.exists():
            status.error = "No campaign directory set"
            return status
        
        # Load campaign_state.json (authoritative)
        state_file = self.outdir / "campaign_state.json"
        if state_file.exists():
            try:
                state = json.loads(state_file.read_text())
                status.state = state.get("campaign_state", "UNKNOWN")
                status.round_index = state.get("round_index", -1)
                status.elapsed_hours = state.get("elapsed_hours", 0.0)
                status.consecutive_empty_rounds = state.get("consecutive_empty_rounds", 0)
                
                # Check for time budget
                if "campaign_start_time" in state:
                    import time
                    elapsed = time.time() - state["campaign_start_time"]
                    status.elapsed_hours = elapsed / 3600.0
                
            except Exception as e:
                status.error = f"Failed to read state: {e}"
        
        # Load campaign_status.json for metrics
        status_file = self.outdir / "campaign_status.json"
        if status_file.exists():
            try:
                snapshot = json.loads(status_file.read_text())
                metrics = snapshot.get("rolling_metrics", {})
                status.yield_val = metrics.get("yield", 0.0)
                status.coverage = metrics.get("coverage", 0.0)
                status.diversity = metrics.get("diversity", 0.0)
                status.composite = metrics.get("composite", 0.0)
            except Exception:
                pass
        
        # Count events
        events_file = self.outdir / "events.jsonl"
        if events_file.exists():
            try:
                with open(events_file) as f:
                    status.total_events = sum(1 for _ in f)
            except Exception:
                pass
        
        # Count run directories
        checkpoint_dir = self.outdir / "checkpoints"
        if checkpoint_dir.exists():
            try:
                status.total_run_dirs = sum(
                    1 for d in checkpoint_dir.glob("**/run_*")
                    if d.is_dir()
                )
            except Exception:
                pass
        
        # Check staleness (files not modified in 5+ minutes)
        try:
            import time
            mtime = state_file.stat().st_mtime if state_file.exists() else 0
            status.is_stale = (time.time() - mtime) > 300  # 5 minutes
            if mtime > 0:
                status.last_activity = datetime.fromtimestamp(mtime)
        except Exception:
            pass
        
        return status
    
    def list_campaigns(self, runs_dir: Path) -> list[Dict[str, Any]]:
        """List all campaign directories with basic info."""
        campaigns = []
        
        if not runs_dir.exists():
            return campaigns
        
        for d in sorted(runs_dir.iterdir()):
            if not d.is_dir():
                continue
            
            state_file = d / "campaign_state.json"
            if not state_file.exists():
                continue
            
            try:
                state = json.loads(state_file.read_text())
                campaigns.append({
                    "name": d.name,
                    "path": str(d),
                    "state": state.get("campaign_state", "UNKNOWN"),
                    "round_index": state.get("round_index", -1),
                })
            except Exception:
                campaigns.append({
                    "name": d.name,
                    "path": str(d),
                    "state": "ERROR",
                    "round_index": -1,
                })
        
        return campaigns
