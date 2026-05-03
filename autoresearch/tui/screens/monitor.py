"""
Monitor screen - live campaign monitoring.
"""
from textual.app import ComposeResult
from textual.containers import Container, Horizontal, Vertical
from textual.widgets import Static, DataTable, Log, ProgressBar, TabbedContent, TabPane
from textual.reactive import reactive

from ..services.monitor import CampaignStatus


class MonitorScreen(Container):
    """Live campaign monitoring screen."""
    
    DEFAULT_CSS = """
    MonitorScreen {
        layout: vertical;
        padding: 1;
    }
    
    .section-title {
        text-style: bold;
        margin-bottom: 1;
        color: $primary;
    }
    
    .status-bar {
        height: 5;
        border: solid $surface;
        padding: 1;
        margin-bottom: 1;
    }
    
    .status-running {
        background: $success 20%;
        border: solid $success;
    }
    
    .status-completed {
        background: $primary 20%;
        border: solid $primary;
    }
    
    .status-failed {
        background: $error 20%;
        border: solid $error;
    }
    
    .metrics-row {
        layout: horizontal;
        height: 4;
        margin-bottom: 1;
    }
    
    .metric-box {
        width: 1fr;
        height: 100%;
        border: solid $surface;
        padding: 0 1;
        margin-right: 1;
        content-align: center middle;
    }
    
    .progress-section {
        height: 3;
        margin-bottom: 1;
    }
    
    .tabs-section {
        height: 1fr;
    }
    
    #events-log {
        height: 100%;
    }
    
    #metrics-table {
        height: 100%;
    }
    """
    
    status: reactive[CampaignStatus] = reactive(CampaignStatus)
    
    def compose(self) -> ComposeResult:
        yield Static("[b]📡 Live Monitor[/]", classes="section-title")
        
        # Status bar
        with Horizontal(classes="status-bar", id="monitor-status-bar"):
            yield Static("⏳ Waiting for campaign...", id="status-text")
        
        # Progress bar (for time budget)
        with Container(classes="progress-section"):
            yield Static("[dim]Time Budget Progress[/]", id="progress-label")
            yield ProgressBar(id="time-progress", total=100, show_eta=False)
        
        # Metrics row
        with Horizontal(classes="metrics-row"):
            yield Static("[b]Yield[/]\n—", id="metric-yield", classes="metric-box")
            yield Static("[b]Coverage[/]\n—", id="metric-coverage", classes="metric-box")
            yield Static("[b]Diversity[/]\n—", id="metric-diversity", classes="metric-box")
            yield Static("[b]Composite[/]\n—", id="metric-composite", classes="metric-box")
        
        # Tabbed content for logs/metrics/alerts
        with TabbedContent(classes="tabs-section"):
            with TabPane("Events", id="tab-events"):
                yield Log(id="events-log", highlight=True, auto_scroll=True)
            
            with TabPane("Metrics History", id="tab-metrics"):
                yield DataTable(id="metrics-table")
            
            with TabPane("Alerts", id="tab-alerts"):
                yield Log(id="alerts-log", highlight=True, auto_scroll=True)
    
    def on_mount(self) -> None:
        """Initialize the metrics table."""
        try:
            table = self.query_one("#metrics-table", DataTable)
            table.add_columns("Round", "Yield", "Coverage", "Diversity", "Composite", "Events")
        except Exception:
            pass
    
    def watch_status(self, status: CampaignStatus) -> None:
        """React to status updates."""
        self._update_display(status)
    
    def _update_display(self, status: CampaignStatus):
        """Update all display elements."""
        try:
            # Status bar
            status_bar = self.query_one("#monitor-status-bar", Horizontal)
            status_text = self.query_one("#status-text", Static)
            
            status_bar.remove_class("status-running", "status-completed", "status-failed")
            
            if status.state == "RUNNING":
                status_bar.add_class("status-running")
                elapsed_str = f"{status.elapsed_hours:.1f}h" if status.elapsed_hours else "0h"
                status_text.update(
                    f"🔄 [b]{status.state}[/] | "
                    f"Round {status.round_index + 1} | "
                    f"Events: {status.total_events:,} | "
                    f"Elapsed: {elapsed_str}"
                )
            elif status.state in ("COMPLETED", "CONVERGED"):
                status_bar.add_class("status-completed")
                status_text.update(
                    f"✅ [b]{status.state}[/] | "
                    f"Rounds: {status.round_index + 1} | "
                    f"Events: {status.total_events:,}"
                )
            elif status.state in ("FAILED", "PREFLIGHT_BLOCKED"):
                status_bar.add_class("status-failed")
                status_text.update(f"❌ [b]{status.state}[/] | {status.error or 'Unknown error'}")
            else:
                status_text.update(f"⏳ {status.state} | Waiting for data...")
            
            # Progress bar
            progress = self.query_one("#time-progress", ProgressBar)
            label = self.query_one("#progress-label", Static)
            
            if status.progress_pct is not None:
                progress.update(progress=status.progress_pct)
                label.update(
                    f"Time Budget: {status.elapsed_hours:.1f}h / {status.max_duration_hours:.1f}h "
                    f"({status.progress_pct:.0f}%)"
                )
            else:
                progress.update(progress=0)
                label.update("[dim]No time budget set[/]")
            
            # Metrics
            self.query_one("#metric-yield", Static).update(f"[b]Yield[/]\n{status.yield_val:.2f}")
            self.query_one("#metric-coverage", Static).update(f"[b]Coverage[/]\n{status.coverage:.1f}%")
            self.query_one("#metric-diversity", Static).update(f"[b]Diversity[/]\n{status.diversity:.2f}")
            self.query_one("#metric-composite", Static).update(f"[b]Composite[/]\n{status.composite:.2f}")
            
        except Exception:
            pass
    
    def add_event(self, message: str):
        """Add an event to the events log."""
        try:
            self.query_one("#events-log", Log).write_line(message)
        except Exception:
            pass
    
    def add_alert(self, message: str):
        """Add an alert to the alerts log."""
        try:
            self.query_one("#alerts-log", Log).write_line(message)
        except Exception:
            pass
    
    def add_metrics_row(self, round_idx: int, metrics: dict):
        """Add a row to the metrics history table."""
        try:
            table = self.query_one("#metrics-table", DataTable)
            table.add_row(
                str(round_idx + 1),
                f"{metrics.get('yield', 0):.2f}",
                f"{metrics.get('coverage', 0):.1f}%",
                f"{metrics.get('diversity', 0):.2f}",
                f"{metrics.get('composite', 0):.2f}",
                str(metrics.get("events", 0)),
            )
        except Exception:
            pass
    
    def clear_logs(self):
        """Clear all logs."""
        try:
            self.query_one("#events-log", Log).clear()
            self.query_one("#alerts-log", Log).clear()
        except Exception:
            pass
