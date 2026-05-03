#!/usr/bin/env python3
"""
Live Campaign Monitor

Real-time campaign monitoring with automatic stale detection.
Replaces manual ps/jq/tail checks with a single command.

Usage:
    python -m autoresearch.harness.live_monitor /path/to/campaign/outdir [--watch]
    
Options:
    --watch         Refresh display every N seconds (default: 10)
    --interval N    Set refresh interval in seconds
    --compact       Minimal output (state + progress only)
"""

import json
import os
import sys
import time
import argparse
from pathlib import Path
from datetime import datetime, timezone
from typing import Dict, Optional, Tuple
import subprocess


class CampaignMonitor:
    """Monitor campaign state and detect staleness."""
    
    def __init__(self, outdir: Path):
        self.outdir = Path(outdir)
        self.state_file = self.outdir / "campaign_state.json"
        self.status_file = self.outdir / "campaign_status.json"
        self.events_file = self.outdir / "events.jsonl"
        self.supervisor_events_file = self.outdir / "supervisor_events.jsonl"
        self.log_file = self.outdir / "supervisor.log"
        
    def load_state(self) -> Optional[Dict]:
        """Load campaign_state.json (authoritative)."""
        if not self.state_file.exists():
            return None
        try:
            with open(self.state_file, 'r') as f:
                return json.load(f)
        except Exception as e:
            return {"error": str(e)}
    
    def load_status(self) -> Optional[Dict]:
        """Load campaign_status.json (last snapshot)."""
        if not self.status_file.exists():
            return None
        try:
            with open(self.status_file, 'r') as f:
                return json.load(f)
        except Exception as e:
            return {"error": str(e)}
    
    def count_events(self) -> int:
        """Count total ATTEMPT_EVALUATED events."""
        if not self.events_file.exists():
            return 0
        try:
            with open(self.events_file, 'r') as f:
                return sum(1 for _ in f)
        except Exception:
            return 0
    
    def count_rounds(self) -> int:
        """Count completed rounds from supervisor_events.jsonl."""
        if not self.supervisor_events_file.exists():
            return 0
        try:
            with open(self.supervisor_events_file, 'r') as f:
                return sum(1 for line in f if '"event_type": "ROUND_COMPLETED"' in line)
        except Exception:
            return 0
    
    def count_run_dirs(self) -> int:
        """Count run_* directories created."""
        checkpoints_dir = self.outdir / "checkpoints"
        if not checkpoints_dir.exists():
            return 0
        try:
            # Find all run_* directories
            count = 0
            for arm_dir in checkpoints_dir.iterdir():
                if arm_dir.is_dir():
                    for iter_dir in arm_dir.iterdir():
                        if iter_dir.is_dir() and iter_dir.name.startswith("iter_"):
                            count += sum(1 for d in iter_dir.iterdir() if d.is_dir() and d.name.startswith("run_"))
            return count
        except Exception:
            return 0
    
    def get_last_activity(self) -> Optional[Tuple[datetime, str]]:
        """Get timestamp of most recently modified run directory."""
        checkpoints_dir = self.outdir / "checkpoints"
        if not checkpoints_dir.exists():
            return None
        
        try:
            latest_time = None
            latest_path = None
            
            for arm_dir in checkpoints_dir.iterdir():
                if arm_dir.is_dir():
                    for iter_dir in arm_dir.iterdir():
                        if iter_dir.is_dir() and iter_dir.name.startswith("iter_"):
                            for run_dir in iter_dir.iterdir():
                                if run_dir.is_dir() and run_dir.name.startswith("run_"):
                                    mtime = run_dir.stat().st_mtime
                                    if latest_time is None or mtime > latest_time:
                                        latest_time = mtime
                                        latest_path = str(run_dir.relative_to(self.outdir))
            
            if latest_time:
                dt = datetime.fromtimestamp(latest_time)
                return (dt, latest_path)
            return None
        except Exception:
            return None
    
    def check_process(self) -> Optional[Dict]:
        """Check if supervisor process is running."""
        try:
            # Search for run_supervisor process with this campaign's outdir
            result = subprocess.run(
                ["ps", "aux"],
                capture_output=True,
                text=True,
                timeout=5
            )
            
            for line in result.stdout.splitlines():
                if "run_supervisor" in line and str(self.outdir) in line and "grep" not in line:
                    parts = line.split()
                    return {
                        "pid": parts[1],
                        "cpu": parts[2],
                        "mem": parts[3],
                        "time": parts[9],
                        "running": True
                    }
            
            return {"running": False}
        except Exception:
            return {"running": False, "error": "Could not check process"}
    
    def detect_stale(self, state: Dict, process_info: Dict, last_activity: Optional[Tuple]) -> Tuple[bool, str]:
        """
        Detect if campaign is stale.
        
        Returns:
            (is_stale, reason)
        """
        campaign_state = state.get("campaign_state", "UNKNOWN")
        
        # Terminal states are not stale
        if campaign_state in ["COMPLETED", "CONVERGED", "FAILED"]:
            return (False, f"Campaign finished: {campaign_state}")
        
        # If state says RUNNING but no process → stale
        if campaign_state == "RUNNING" and not process_info.get("running"):
            return (True, "State says RUNNING but no supervisor process found")
        
        # If process running, check last activity
        if process_info.get("running") and last_activity:
            last_time, _ = last_activity
            age_seconds = (datetime.now() - last_time).total_seconds()
            
            # If no activity in 5 minutes → likely stale/stuck
            if age_seconds > 300:
                return (True, f"No file activity in {int(age_seconds)}s (>5min)")
        
        # Process running and recent activity → active
        if process_info.get("running"):
            return (False, "Active (process running, recent activity)")
        
        # Unknown state
        return (False, "State unclear")
    
    def get_metrics_summary(self, status: Optional[Dict]) -> str:
        """Extract metrics from status snapshot."""
        if not status or "rolling_metrics" not in status:
            return "No metrics available"
        
        metrics = status["rolling_metrics"]
        return (
            f"Yield: {metrics.get('yield', 0):.4f} | "
            f"Coverage: {metrics.get('coverage', 0):.2f} | "
            f"Diversity: {metrics.get('diversity', 0):.2f} | "
            f"Composite: {metrics.get('composite', 0):.2f}"
        )
    
    def format_report(self, compact: bool = False) -> str:
        """Generate monitoring report."""
        state = self.load_state()
        status = self.load_status()
        process_info = self.check_process()
        last_activity = self.get_last_activity()
        
        if not state:
            return f"❌ Campaign not found: {self.outdir}"
        
        campaign_state = state.get("campaign_state", "UNKNOWN")
        round_number = state.get("round_index", 0)
        
        # Staleness detection
        is_stale, stale_reason = self.detect_stale(state, process_info, last_activity)
        
        # Event/round counts
        events_count = self.count_events()
        rounds_count = self.count_rounds()
        run_dirs_count = self.count_run_dirs()
        
        # Process status
        process_running = process_info.get("running", False)
        pid = process_info.get("pid", "N/A")
        
        # Last activity
        if last_activity:
            last_time, last_path = last_activity
            age = datetime.now() - last_time
            age_str = f"{int(age.total_seconds())}s ago"
        else:
            age_str = "No activity"
            last_path = "N/A"
        
        # Build report
        lines = []
        
        if compact:
            # Compact mode: single line status
            status_icon = "✅" if campaign_state in ["COMPLETED", "CONVERGED"] else "🔄" if process_running else "⚠️"
            stale_icon = "🚨" if is_stale else ""
            lines.append(
                f"{status_icon}{stale_icon} {campaign_state} | Round {round_number} | "
                f"Events: {events_count} | Runs: {run_dirs_count} | PID: {pid} | Last: {age_str}"
            )
        else:
            # Full report
            lines.append("=" * 80)
            lines.append(f"Campaign Monitor: {self.outdir.name}")
            lines.append("=" * 80)
            lines.append("")
            
            # State
            state_icon = "✅" if campaign_state in ["COMPLETED", "CONVERGED"] else "🔄" if process_running else "⚠️"
            lines.append(f"State:       {state_icon} {campaign_state}")
            lines.append(f"Round:       {round_number}")
            
            # Staleness
            if is_stale:
                lines.append(f"Status:      🚨 STALE - {stale_reason}")
            else:
                lines.append(f"Status:      {stale_reason}")
            
            lines.append("")
            
            # Progress
            lines.append("Progress:")
            lines.append(f"  Events:         {events_count}")
            lines.append(f"  Rounds:         {rounds_count}")
            lines.append(f"  Run Dirs:       {run_dirs_count}")
            
            lines.append("")
            
            # Process
            lines.append("Supervisor Process:")
            if process_running:
                lines.append(f"  PID:            {pid}")
                lines.append(f"  CPU:            {process_info.get('cpu', 'N/A')}%")
                lines.append(f"  Memory:         {process_info.get('mem', 'N/A')}%")
                lines.append(f"  Time:           {process_info.get('time', 'N/A')}")
            else:
                lines.append(f"  Status:         ❌ Not running")
            
            lines.append("")
            
            # Last activity
            lines.append("Last Activity:")
            lines.append(f"  Time:           {age_str}")
            if last_activity:
                lines.append(f"  Path:           {last_path}")
            
            lines.append("")
            
            # Metrics
            if status:
                lines.append("Metrics:")
                lines.append(f"  {self.get_metrics_summary(status)}")
                lines.append("")
            
            # Convergence info
            if state.get("convergence_state"):
                lines.append("Convergence:")
                lines.append(f"  State:          {state.get('convergence_state')}")
                if state.get("convergence_reason"):
                    lines.append(f"  Reason:         {state.get('convergence_reason')}")
                lines.append("")
            
            lines.append("=" * 80)
        
        return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Live campaign monitoring with automatic stale detection",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single status check
  python -m autoresearch.harness.live_monitor /path/to/campaign
  
  # Watch mode (refresh every 10s)
  python -m autoresearch.harness.live_monitor /path/to/campaign --watch
  
  # Compact watch (minimal output)
  python -m autoresearch.harness.live_monitor /path/to/campaign --watch --compact
  
  # Custom refresh interval
  python -m autoresearch.harness.live_monitor /path/to/campaign --watch --interval 5
        """
    )
    
    parser.add_argument(
        "outdir",
        type=Path,
        help="Campaign output directory"
    )
    parser.add_argument(
        "--watch",
        action="store_true",
        help="Refresh display periodically"
    )
    parser.add_argument(
        "--interval",
        type=int,
        default=10,
        help="Refresh interval in seconds (default: 10)"
    )
    parser.add_argument(
        "--compact",
        action="store_true",
        help="Minimal output (single line)"
    )
    
    args = parser.parse_args()
    
    if not args.outdir.exists():
        print(f"❌ Campaign directory not found: {args.outdir}", file=sys.stderr)
        sys.exit(1)
    
    monitor = CampaignMonitor(args.outdir)
    
    if args.watch:
        # Watch mode
        try:
            iteration = 0
            while True:
                if not args.compact and iteration > 0:
                    # Clear screen for full report (skip first iteration)
                    print("\033[2J\033[H", end="")
                
                print(monitor.format_report(compact=args.compact))
                
                if args.compact:
                    sys.stdout.flush()
                
                iteration += 1
                time.sleep(args.interval)
        except KeyboardInterrupt:
            print("\n\n👋 Monitoring stopped.")
            sys.exit(0)
    else:
        # Single check
        print(monitor.format_report(compact=args.compact))


if __name__ == "__main__":
    main()
