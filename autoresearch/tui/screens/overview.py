"""
Overview screen - campaign summary and quick status.
"""
from textual.app import ComposeResult
from textual.containers import Container, Horizontal, Vertical
from textual.widgets import Static, Label, Button
from textual.reactive import reactive

from ..services.monitor import CampaignStatus


class MetricCard(Static):
    """A card displaying a single metric."""
    
    def __init__(self, title: str, value: str = "—", **kwargs):
        super().__init__(**kwargs)
        self.title = title
        self._value = value
    
    def compose(self) -> ComposeResult:
        yield Static(f"[dim]{self.title}[/]", classes="metric-title")
        yield Static(self._value, classes="metric-value", id=f"value-{self.id}")
    
    def update_value(self, value: str):
        """Update the displayed value."""
        self._value = value
        value_widget = self.query_one(f"#value-{self.id}", Static)
        value_widget.update(value)


class OverviewScreen(Container):
    """Overview screen showing campaign summary."""
    
    DEFAULT_CSS = """
    OverviewScreen {
        layout: vertical;
        padding: 1;
    }
    
    .section-title {
        text-style: bold;
        margin-bottom: 1;
        color: $primary;
    }
    
    .metric-row {
        layout: horizontal;
        height: auto;
        margin-bottom: 1;
    }
    
    .metric-card {
        width: 1fr;
        height: 5;
        border: solid $surface;
        padding: 0 1;
        margin-right: 1;
    }
    
    .metric-title {
        height: 1;
    }
    
    .metric-value {
        height: 2;
        text-style: bold;
        content-align: center middle;
    }
    
    .status-running {
        background: $success 20%;
        border: solid $success;
    }
    
    .status-failed {
        background: $error 20%;
        border: solid $error;
    }
    
    .status-completed {
        background: $primary 20%;
        border: solid $primary;
    }
    
    .info-section {
        margin-top: 1;
        padding: 1;
        border: solid $surface;
        height: auto;
    }
    
    .actions-row {
        layout: horizontal;
        height: 3;
        margin-top: 1;
    }
    
    .actions-row Button {
        margin-right: 1;
    }
    """
    
    status: reactive[CampaignStatus] = reactive(CampaignStatus)
    
    def compose(self) -> ComposeResult:
        yield Static("[b]📊 Campaign Overview[/]", classes="section-title")
        
        # Status cards row
        with Horizontal(classes="metric-row"):
            yield MetricCard("State", "—", id="state-card", classes="metric-card")
            yield MetricCard("Round", "—", id="round-card", classes="metric-card")
            yield MetricCard("Events", "—", id="events-card", classes="metric-card")
            yield MetricCard("Run Dirs", "—", id="runs-card", classes="metric-card")
        
        # Metrics row
        with Horizontal(classes="metric-row"):
            yield MetricCard("Yield", "—", id="yield-card", classes="metric-card")
            yield MetricCard("Coverage", "—", id="coverage-card", classes="metric-card")
            yield MetricCard("Diversity", "—", id="diversity-card", classes="metric-card")
            yield MetricCard("Composite", "—", id="composite-card", classes="metric-card")
        
        # Info section
        with Container(classes="info-section"):
            yield Static("[b]Campaign Info[/]", classes="section-title")
            yield Static("No campaign selected. Go to [b]Run[/] to start a campaign.", id="campaign-info")
        
        # Quick actions
        with Horizontal(classes="actions-row"):
            yield Button("📁 Select Campaign", id="btn-select", variant="default")
            yield Button("▶️ Go to Run", id="btn-run", variant="primary")
            yield Button("📊 View Results", id="btn-results", variant="success")
    
    def watch_status(self, status: CampaignStatus) -> None:
        """React to status changes."""
        self._update_cards(status)
    
    def _update_cards(self, status: CampaignStatus):
        """Update all metric cards with new status."""
        try:
            # State card with styling
            state_card = self.query_one("#state-card", MetricCard)
            state_card.update_value(f"{status.state_emoji} {status.state}")
            state_card.remove_class("status-running", "status-failed", "status-completed")
            if status.state == "RUNNING":
                state_card.add_class("status-running")
            elif status.state in ("FAILED", "PREFLIGHT_BLOCKED"):
                state_card.add_class("status-failed")
            elif status.state in ("COMPLETED", "CONVERGED"):
                state_card.add_class("status-completed")
            
            # Other cards
            self.query_one("#round-card", MetricCard).update_value(str(status.round_index + 1))
            self.query_one("#events-card", MetricCard).update_value(f"{status.total_events:,}")
            self.query_one("#runs-card", MetricCard).update_value(f"{status.total_run_dirs:,}")
            
            self.query_one("#yield-card", MetricCard).update_value(f"{status.yield_val:.2f}")
            self.query_one("#coverage-card", MetricCard).update_value(f"{status.coverage:.1f}%")
            self.query_one("#diversity-card", MetricCard).update_value(f"{status.diversity:.2f}")
            self.query_one("#composite-card", MetricCard).update_value(f"{status.composite:.2f}")
            
        except Exception:
            pass  # Widget not yet mounted
    
    def set_campaign_info(self, info: str):
        """Update the campaign info text."""
        try:
            self.query_one("#campaign-info", Static).update(info)
        except Exception:
            pass
