"""
Autoresearch TUI - Unified Terminal User Interface for campaign management.

This is the main application entry point combining all screens into a
cohesive dashboard similar to htop/btop.
"""
from pathlib import Path
from typing import Optional, Iterable
import json

from textual.app import App, ComposeResult, SystemCommand
from textual.binding import Binding
from textual.containers import Container, Horizontal, Vertical
from textual.screen import Screen
from textual.widgets import (
    Static, Footer, Header, OptionList, ContentSwitcher,
    Button, Label
)
from textual.widgets.option_list import Option
from textual.timer import Timer

from .screens.overview import OverviewScreen
from .screens.config import ConfigScreen
from .screens.run import RunScreen
from .screens.monitor import MonitorScreen
from .screens.results import ResultsScreen
from .services.campaign import CampaignService, CampaignProcess
from .services.monitor import MonitoringService, CampaignStatus


class Sidebar(Container):
    """Navigation sidebar with screen selection."""
    
    DEFAULT_CSS = """
    Sidebar {
        width: 24;
        height: 100%;
        dock: left;
        border-right: solid $surface;
        padding: 1;
    }
    
    Sidebar Static.title {
        text-style: bold;
        color: $primary;
        margin-bottom: 1;
        text-align: center;
    }
    
    Sidebar OptionList {
        height: 1fr;
        border: none;
        padding: 0;
    }
    
    Sidebar .status-section {
        height: auto;
        margin-top: 1;
        padding: 1;
        border-top: solid $surface;
    }
    
    Sidebar .status-line {
        height: 1;
    }
    """
    
    def compose(self) -> ComposeResult:
        yield Static("🔬 Autoresearch", classes="title")
        yield OptionList(
            Option("📊 Overview", id="overview"),
            Option("⚙️  Config", id="config"),
            Option("▶️  Run", id="run"),
            Option("📡 Monitor", id="monitor"),
            Option("📈 Results", id="results"),
            id="nav-list"
        )
        with Container(classes="status-section"):
            yield Static("State: [dim]Idle[/]", id="sidebar-state", classes="status-line")
            yield Static("Round: [dim]—[/]", id="sidebar-round", classes="status-line")
    
    def update_status(self, status: CampaignStatus):
        """Update sidebar status display."""
        try:
            self.query_one("#sidebar-state", Static).update(
                f"State: {status.state_emoji} {status.state}"
            )
            self.query_one("#sidebar-round", Static).update(
                f"Round: [b]{status.round_index + 1}[/]"
            )
        except Exception:
            pass


class MainContent(Container):
    """Main content area with screen switcher."""
    
    DEFAULT_CSS = """
    MainContent {
        width: 1fr;
        height: 100%;
    }
    
    MainContent ContentSwitcher {
        width: 100%;
        height: 100%;
    }
    """
    
    def compose(self) -> ComposeResult:
        with ContentSwitcher(initial="overview", id="content-switcher"):
            yield OverviewScreen(id="overview")
            yield ConfigScreen(id="config")
            yield RunScreen(id="run")
            yield MonitorScreen(id="monitor")
            yield ResultsScreen(id="results")


class StatusBar(Container):
    """Bottom status bar showing quick info."""
    
    DEFAULT_CSS = """
    StatusBar {
        dock: bottom;
        height: 1;
        background: $surface;
        padding: 0 1;
    }
    
    StatusBar Static {
        width: 1fr;
    }
    
    StatusBar .left {
        text-align: left;
    }
    
    StatusBar .right {
        text-align: right;
    }
    """
    
    def compose(self) -> ComposeResult:
        yield Static("[dim]Campaign: None[/]", id="status-campaign", classes="left")
        yield Static("[dim]Press ? for help[/]", classes="right")
    
    def set_campaign(self, name: str):
        """Update campaign name in status bar."""
        try:
            self.query_one("#status-campaign", Static).update(f"Campaign: [b]{name}[/]")
        except Exception:
            pass


class AutoresearchApp(App):
    """
    Autoresearch TUI Application.
    
    A unified dashboard for managing DiHiggs autoresearch campaigns.
    """
    
    CSS = """
    Screen {
        layout: vertical;
    }
    
    .main-layout {
        layout: horizontal;
        height: 1fr;
    }
    """
    
    BINDINGS = [
        Binding("q", "quit", "Quit"),
        Binding("?", "help", "Help"),
        Binding("1", "switch_screen('overview')", "Overview", show=False),
        Binding("2", "switch_screen('config')", "Config", show=False),
        Binding("3", "switch_screen('run')", "Run", show=False),
        Binding("4", "switch_screen('monitor')", "Monitor", show=False),
        Binding("5", "switch_screen('results')", "Results", show=False),
        Binding("r", "refresh", "Refresh"),
        Binding("escape", "focus_sidebar", "Sidebar", show=False),
    ]
    
    TITLE = "Autoresearch TUI"
    SUB_TITLE = "DiHiggs Campaign Manager"
    
    # Services
    campaign_service: Optional[CampaignService] = None
    monitoring_service: Optional[MonitoringService] = None
    
    # State
    current_campaign_dir: Optional[Path] = None
    current_config: Optional[dict] = None
    
    # Polling timer
    _monitor_timer: Optional[Timer] = None
    
    def compose(self) -> ComposeResult:
        yield Header()
        with Horizontal(classes="main-layout"):
            yield Sidebar(id="sidebar")
            yield MainContent(id="main-content")
        yield StatusBar(id="status-bar")
        yield Footer()
    
    def on_mount(self) -> None:
        """Initialize app on mount."""
        # Select overview by default
        nav_list = self.query_one("#nav-list", OptionList)
        nav_list.highlighted = 0
        
        # Start monitoring timer (polls every 2 seconds when campaign is running)
        self._monitor_timer = self.set_interval(2.0, self._poll_campaign_status, pause=True)
    
    def on_option_list_option_selected(self, event: OptionList.OptionSelected) -> None:
        """Handle navigation selection."""
        screen_id = str(event.option.id)
        self._switch_to_screen(screen_id)
    
    def _switch_to_screen(self, screen_id: str) -> None:
        """Switch to a specific screen."""
        switcher = self.query_one("#content-switcher", ContentSwitcher)
        switcher.current = screen_id
        
        # Update nav highlight
        nav_list = self.query_one("#nav-list", OptionList)
        screen_indices = {"overview": 0, "config": 1, "run": 2, "monitor": 3, "results": 4}
        if screen_id in screen_indices:
            nav_list.highlighted = screen_indices[screen_id]
    
    def action_switch_screen(self, screen_id: str) -> None:
        """Action to switch screens (for keybindings)."""
        self._switch_to_screen(screen_id)
    
    def action_focus_sidebar(self) -> None:
        """Focus the sidebar navigation."""
        self.query_one("#nav-list", OptionList).focus()
    
    def action_refresh(self) -> None:
        """Refresh current view."""
        if self.current_campaign_dir:
            self._poll_campaign_status()
    
    def action_help(self) -> None:
        """Show help dialog."""
        self.notify(
            "Keys: 1-5 switch screens, r refresh, q quit, Ctrl+P command palette",
            title="Help",
            timeout=5
        )
    
    def get_system_commands(self, screen: Screen) -> Iterable[SystemCommand]:
        """Return system commands for command palette (Ctrl+P)."""
        # Always include default commands
        yield from super().get_system_commands(screen)
        
        # Navigation commands
        yield SystemCommand(
            "Go to Overview",
            "Switch to Overview screen",
            lambda: self._switch_to_screen("overview"),
        )
        yield SystemCommand(
            "Go to Config",
            "Switch to Config screen",
            lambda: self._switch_to_screen("config"),
        )
        yield SystemCommand(
            "Go to Run",
            "Switch to Run screen",
            lambda: self._switch_to_screen("run"),
        )
        yield SystemCommand(
            "Go to Monitor",
            "Switch to Monitor screen",
            lambda: self._switch_to_screen("monitor"),
        )
        yield SystemCommand(
            "Go to Results",
            "Switch to Results screen",
            lambda: self._switch_to_screen("results"),
        )
        
        # Config screen commands
        yield SystemCommand(
            "Save Config",
            "Save current configuration to file",
            self._command_save_config,
        )
        yield SystemCommand(
            "Load Config",
            "Load configuration from file",
            self._command_load_config,
        )
        yield SystemCommand(
            "Reset Config",
            "Reset configuration to defaults",
            self._command_reset_config,
        )
        
        # Run screen commands
        yield SystemCommand(
            "Browse Config File",
            "Select a config file for the campaign",
            self._command_browse_config,
        )
        yield SystemCommand(
            "Start Campaign",
            "Start the campaign with current config",
            self._command_start_campaign,
        )
        yield SystemCommand(
            "Stop Campaign",
            "Stop the running campaign",
            self.stop_campaign,
        )
        
        # Results screen commands
        yield SystemCommand(
            "Load Campaign Results",
            "Load results from a campaign directory",
            self._command_load_campaign,
        )
        yield SystemCommand(
            "Export Summary CSV",
            "Export campaign summary to CSV file",
            self._command_export_csv,
        )
        yield SystemCommand(
            "Generate Plots",
            "Generate result plots for current campaign",
            self._command_generate_plots,
        )
    
    # Command palette action methods
    
    def _command_save_config(self) -> None:
        """Trigger save config button."""
        config_screen = self.query_one("#config", ConfigScreen)
        save_btn = config_screen.query_one("#btn-save", Button)
        save_btn.press()
    
    def _command_load_config(self) -> None:
        """Trigger load config button."""
        config_screen = self.query_one("#config", ConfigScreen)
        load_btn = config_screen.query_one("#btn-load", Button)
        load_btn.press()
    
    def _command_reset_config(self) -> None:
        """Trigger reset config button."""
        config_screen = self.query_one("#config", ConfigScreen)
        reset_btn = config_screen.query_one("#btn-reset", Button)
        reset_btn.press()
    
    def _command_browse_config(self) -> None:
        """Trigger browse config button."""
        run_screen = self.query_one("#run", RunScreen)
        browse_btn = run_screen.query_one("#btn-browse-config", Button)
        browse_btn.press()
    
    def _command_start_campaign(self) -> None:
        """Trigger start campaign button."""
        run_screen = self.query_one("#run", RunScreen)
        start_btn = run_screen.query_one("#btn-start", Button)
        start_btn.press()
    
    def _command_load_campaign(self) -> None:
        """Trigger load campaign button."""
        results_screen = self.query_one("#results", ResultsScreen)
        load_btn = results_screen.query_one("#btn-load", Button)
        load_btn.press()
    
    def _command_export_csv(self) -> None:
        """Trigger export CSV button."""
        results_screen = self.query_one("#results", ResultsScreen)
        export_btn = results_screen.query_one("#btn-export-csv", Button)
        export_btn.press()
    
    def _command_generate_plots(self) -> None:
        """Trigger generate plots button."""
        results_screen = self.query_one("#results", ResultsScreen)
        plots_btn = results_screen.query_one("#btn-gen-plots", Button)
        plots_btn.press()
    
    # === Campaign Management ===
    
    def set_config(self, config: dict) -> None:
        """Set the current configuration."""
        self.current_config = config
        campaign_id = config.get("campaign_id", "unknown")
        self.query_one("#status-bar", StatusBar).set_campaign(campaign_id)
        
        # Update config screen
        config_screen = self.query_one("#config", ConfigScreen)
        config_screen.load_config(config)
    
    def start_campaign(self, config: dict, outdir: Path) -> bool:
        """Start a new campaign."""
        self.current_config = config
        self.current_campaign_dir = outdir
        
        # Initialize services
        self.campaign_service = CampaignService()
        self.monitoring_service = MonitoringService(outdir)
        
        # Save config to temp file for subprocess
        import tempfile
        config_path = Path(tempfile.mktemp(suffix=".json", prefix="campaign_"))
        config_path.write_text(json.dumps(config, indent=2))
        
        # Create and start campaign subprocess
        def on_state_change(proc: CampaignProcess) -> None:
            run_screen = self.query_one("#run", RunScreen)
            run_screen.process_state = proc.state
        
        self.campaign_service.create_campaign(config_path, outdir, on_state_change=on_state_change)
        success = self.campaign_service.start_campaign()
        
        if success:
            # Update UI
            campaign_id = config.get("campaign_id", "unknown")
            self.query_one("#status-bar", StatusBar).set_campaign(campaign_id)
            
            # Update run screen
            run_screen = self.query_one("#run", RunScreen)
            run_screen.set_running(True, str(outdir))
            
            # Start polling
            if self._monitor_timer:
                self._monitor_timer.resume()
            
            # Switch to monitor screen
            self._switch_to_screen("monitor")
            
            self.notify(f"Campaign started: {campaign_id}", title="Started")
        else:
            self.notify("Failed to start campaign", title="Error", severity="error")
        
        return success
    
    def stop_campaign(self) -> None:
        """Stop the running campaign."""
        if self.campaign_service:
            self.campaign_service.stop_campaign()
            
            # Pause polling
            if self._monitor_timer:
                self._monitor_timer.pause()
            
            # Update run screen
            run_screen = self.query_one("#run", RunScreen)
            run_screen.set_running(False)
            
            self.notify("Campaign stopped", title="Stopped")
    
    def load_campaign(self, campaign_dir: Path) -> bool:
        """Load an existing campaign for viewing."""
        self.current_campaign_dir = campaign_dir
        
        # Initialize monitoring service
        self.monitoring_service = MonitoringService(campaign_dir)
        
        # Load config if exists
        config_path = campaign_dir / "campaign_manifest.json"
        if config_path.exists():
            try:
                with open(config_path) as f:
                    manifest = json.load(f)
                    self.current_config = manifest.get("config", {})
                    self.set_config(self.current_config)
            except Exception:
                pass
        
        # Update status bar
        campaign_id = self.current_config.get("campaign_id", campaign_dir.name) if self.current_config else campaign_dir.name
        self.query_one("#status-bar", StatusBar).set_campaign(campaign_id)
        
        # Load results
        results_screen = self.query_one("#results", ResultsScreen)
        results_screen.load_results(campaign_dir)
        
        # Poll initial status
        self._poll_campaign_status()
        
        return True
    
    def _poll_campaign_status(self) -> None:
        """Poll campaign status and update displays."""
        if not self.monitoring_service:
            return
        
        status = self.monitoring_service.get_status()
        
        # Update sidebar
        sidebar = self.query_one("#sidebar", Sidebar)
        sidebar.update_status(status)
        
        # Update overview screen
        overview = self.query_one("#overview", OverviewScreen)
        overview.status = status
        
        # Update monitor screen
        monitor = self.query_one("#monitor", MonitorScreen)
        monitor.status = status
        
        # Check if campaign completed
        if self.campaign_service and status.state in ("COMPLETED", "CONVERGED", "FAILED"):
            if self._monitor_timer:
                self._monitor_timer.pause()
            
            run_screen = self.query_one("#run", RunScreen)
            run_screen.set_running(False)
            
            self.notify(f"Campaign {status.state.lower()}", title="Campaign Ended")
    
    # === Event Handlers ===
    
    def on_button_pressed(self, event: Button.Pressed) -> None:
        """Handle button presses from any screen."""
        button_id = event.button.id
        
        # Overview screen buttons
        if button_id == "btn-run":
            self._switch_to_screen("run")
        elif button_id == "btn-results":
            self._switch_to_screen("results")
        elif button_id == "btn-select":
            # Open campaign directory picker
            async def _select_campaign():
                from autoresearch.tui.widgets import FilePickerModal
                path = await self.push_screen(
                    FilePickerModal(
                        title="Select Campaign Directory",
                        initial_path="runs",
                        select_dirs=True,
                    )
                )
                if path:
                    self.load_campaign(path)
                    self.notify(f"Campaign loaded: {path.name}", title="Campaign")
            
            self.call_later(_select_campaign)
        
        # Run screen buttons
        elif button_id == "btn-start":
            if self.current_config:
                outdir = Path(self.current_config.get("paths", {}).get("outdir", "runs/tui-campaign"))
                self.start_campaign(self.current_config, outdir)
            else:
                self.notify("Load a config first", title="Error", severity="error")
                self._switch_to_screen("config")
        
        elif button_id == "btn-stop":
            self.stop_campaign()
        
        # Results screen buttons
        elif button_id == "btn-load":
            if self.current_campaign_dir:
                results_screen = self.query_one("#results", ResultsScreen)
                results_screen.load_results(self.current_campaign_dir)
        
        elif button_id == "btn-refresh":
            self._poll_campaign_status()


def run_tui(config_path: Optional[Path] = None, campaign_dir: Optional[Path] = None):
    """
    Run the Autoresearch TUI application.
    
    Args:
        config_path: Optional path to configuration JSON file
        campaign_dir: Optional path to existing campaign directory to load
    """
    app = AutoresearchApp()
    
    # Load config if provided
    if config_path and config_path.exists():
        try:
            with open(config_path) as f:
                config = json.load(f)
            app.set_config(config)
        except Exception as e:
            print(f"Warning: Failed to load config: {e}")
    
    # Load campaign if provided
    if campaign_dir and campaign_dir.exists():
        app.current_campaign_dir = campaign_dir
        app.load_campaign(campaign_dir)
    
    app.run()


def main():
    """CLI entry point."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Autoresearch TUI - Unified campaign management interface",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Launch empty TUI
  python -m autoresearch.tui
  
  # Launch with config
  python -m autoresearch.tui --config configs/my_config.json
  
  # Load existing campaign
  python -m autoresearch.tui --campaign runs/my-campaign
  
  # Both config and campaign
  python -m autoresearch.tui --config configs/lam1.json --campaign runs/lam1-test
        """
    )
    
    parser.add_argument(
        "--config", "-c",
        type=Path,
        help="Path to campaign configuration JSON file"
    )
    parser.add_argument(
        "--campaign", "-d",
        type=Path,
        help="Path to existing campaign directory to load"
    )
    
    args = parser.parse_args()
    
    run_tui(config_path=args.config, campaign_dir=args.campaign)


if __name__ == "__main__":
    main()
