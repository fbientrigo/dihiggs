"""
Run screen - campaign execution control.
"""
from pathlib import Path
from typing import Optional
import json

from textual.app import ComposeResult
from textual.containers import Container, Horizontal, Vertical
from textual.widgets import Static, Button, Input, DirectoryTree, Log
from textual.reactive import reactive

from .config import ConfigScreen
from ..services.campaign import CampaignProcess, ProcessState


class RunScreen(Container):
    """Campaign run control screen."""
    
    DEFAULT_CSS = """
    RunScreen {
        layout: vertical;
        padding: 1;
    }
    
    .section-title {
        text-style: bold;
        margin-bottom: 1;
        color: $primary;
    }
    
    .config-section {
        height: auto;
        padding: 1;
        border: solid $surface;
        margin-bottom: 1;
    }
    
    .form-row {
        layout: horizontal;
        height: 3;
        margin-bottom: 1;
    }
    
    .form-label {
        width: 15;
        padding-top: 1;
    }
    
    .form-input {
        width: 1fr;
    }
    
    .actions-row {
        layout: horizontal;
        height: 3;
        margin-bottom: 1;
    }
    
    .actions-row Button {
        margin-right: 1;
    }
    
    .status-section {
        height: auto;
        padding: 1;
        border: solid $surface;
        margin-bottom: 1;
    }
    
    .status-running {
        background: $success 20%;
        border: solid $success;
    }
    
    .status-stopped {
        background: $surface;
    }
    
    .status-failed {
        background: $error 20%;
        border: solid $error;
    }
    
    .log-section {
        height: 1fr;
        border: solid $surface;
    }
    
    #run-log {
        height: 100%;
    }
    """
    
    process_state: reactive[ProcessState] = reactive(ProcessState.IDLE)
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.config_path: Optional[Path] = None
        self.outdir: Optional[Path] = None
    
    def compose(self) -> ComposeResult:
        yield Static("[b]▶️ Run Campaign[/]", classes="section-title")
        
        # Config selection
        with Container(classes="config-section"):
            yield Static("[b]Configuration[/]")
            
            with Horizontal(classes="form-row"):
                yield Static("Config File:", classes="form-label")
                yield Input(
                    placeholder="Path to config JSON...",
                    id="input-config-path",
                    classes="form-input"
                )
            
            with Horizontal(classes="form-row"):
                yield Static("Output Dir:", classes="form-label")
                yield Input(
                    placeholder="Path for campaign output (runs/my-campaign)",
                    id="input-outdir",
                    classes="form-input"
                )
            
            with Horizontal(classes="actions-row"):
                yield Button("📂 Browse Config", id="btn-browse-config", variant="default")
                yield Button("📋 Use Current Config", id="btn-use-current", variant="default")
        
        # Process control
        with Container(classes="status-section", id="status-container"):
            yield Static("[b]Process Control[/]")
            yield Static("Status: [dim]Idle[/]", id="process-status")
            
            with Horizontal(classes="actions-row"):
                yield Button("▶️ Start Campaign", id="btn-start", variant="success")
                yield Button("⏹️ Stop Campaign", id="btn-stop", variant="error", disabled=True)
                yield Button("🔄 Refresh", id="btn-refresh", variant="default")
        
        # Log output
        with Container(classes="log-section"):
            yield Static("[b]Process Log[/]")
            yield Log(id="run-log", highlight=True, auto_scroll=True)
    
    def watch_process_state(self, state: ProcessState) -> None:
        """React to process state changes."""
        self._update_status_display(state)
    
    def _update_status_display(self, state: ProcessState):
        """Update the status display based on process state."""
        try:
            status_widget = self.query_one("#process-status", Static)
            container = self.query_one("#status-container", Container)
            start_btn = self.query_one("#btn-start", Button)
            stop_btn = self.query_one("#btn-stop", Button)
            
            # Remove all status classes
            container.remove_class("status-running", "status-stopped", "status-failed")
            
            if state == ProcessState.IDLE:
                status_widget.update("Status: [dim]Idle - Ready to start[/]")
                container.add_class("status-stopped")
                start_btn.disabled = False
                stop_btn.disabled = True
                
            elif state == ProcessState.STARTING:
                status_widget.update("Status: [yellow]Starting...[/]")
                start_btn.disabled = True
                stop_btn.disabled = True
                
            elif state == ProcessState.RUNNING:
                status_widget.update("Status: [green]🔄 Running[/]")
                container.add_class("status-running")
                start_btn.disabled = True
                stop_btn.disabled = False
                
            elif state == ProcessState.STOPPING:
                status_widget.update("Status: [yellow]Stopping...[/]")
                start_btn.disabled = True
                stop_btn.disabled = True
                
            elif state == ProcessState.STOPPED:
                status_widget.update("Status: [blue]✅ Stopped[/]")
                container.add_class("status-stopped")
                start_btn.disabled = False
                stop_btn.disabled = True
                
            elif state == ProcessState.FAILED:
                status_widget.update("Status: [red]❌ Failed[/]")
                container.add_class("status-failed")
                start_btn.disabled = False
                stop_btn.disabled = True
                
        except Exception:
            pass
    
    def set_config_path(self, path: Path):
        """Set the config file path."""
        self.config_path = path
        try:
            self.query_one("#input-config-path", Input).value = str(path)
        except Exception:
            pass
    
    def set_outdir(self, path: Path):
        """Set the output directory."""
        self.outdir = path
        try:
            self.query_one("#input-outdir", Input).value = str(path)
        except Exception:
            pass
    
    def get_config_path(self) -> Optional[Path]:
        """Get config path from input."""
        try:
            value = self.query_one("#input-config-path", Input).value
            return Path(value) if value else None
        except Exception:
            return self.config_path
    
    def get_outdir(self) -> Optional[Path]:
        """Get output dir from input."""
        try:
            value = self.query_one("#input-outdir", Input).value
            return Path(value) if value else None
        except Exception:
            return self.outdir
    
    def write_log(self, message: str):
        """Add a message to the run log."""
        try:
            self.query_one("#run-log", Log).write_line(message)
        except Exception:
            pass
    
    def clear_log(self):
        """Clear the run log."""
        try:
            self.query_one("#run-log", Log).clear()
        except Exception:
            pass
    
    def on_button_pressed(self, event: Button.Pressed) -> None:
        """Handle button presses."""
        if event.button.id == "btn-browse-config":
            # Open file picker for config
            async def _browse():
                from autoresearch.tui.widgets import FilePickerModal
                path = await self.app.push_screen(
                    FilePickerModal(
                        title="Select Config File",
                        initial_path="autoresearch/configs",
                        filter_ext=[".json"],
                    )
                )
                if path:
                    self.set_config_path(path)
                    self.app.notify(f"Config selected: {path.name}", title="Config")
            
            self.app.call_later(_browse)
        
        elif event.button.id == "btn-use-current":
            # Use config from config screen
            try:
                config_screen = self.app.query_one("#config", ConfigScreen)
                config = config_screen.get_config()
                if config:
                    # Save to temp file and set path
                    import tempfile
                    from pathlib import Path
                    tmp = Path(tempfile.mktemp(suffix=".json"))
                    tmp.write_text(json.dumps(config, indent=2))
                    self.set_config_path(tmp)
                    self.app.notify("Using config from Config screen", title="Config")
                else:
                    self.app.notify("No config available in Config screen", title="Warning", severity="warning")
            except Exception as e:
                self.app.notify(f"Failed to load config: {e}", title="Error", severity="error")

    def set_running(self, running: bool, outdir: Optional[str] = None) -> None:
        """Update the running state of the campaign."""
        if running:
            self.process_state = ProcessState.RUNNING
            if outdir:
                self.set_outdir(Path(outdir))
        else:
            self.process_state = ProcessState.STOPPED
