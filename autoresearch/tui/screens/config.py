"""
Config screen - campaign configuration wizard.
"""
import json
from pathlib import Path
from typing import Dict, Any, Optional

from textual.app import ComposeResult
from textual.containers import Container, Horizontal, Vertical, ScrollableContainer
from textual.widgets import Static, Input, Button, Select, Checkbox, TextArea, Collapsible
from textual.reactive import reactive


# Default config template
DEFAULT_CONFIG = {
    "supervisor": {
        "max_rounds": 0,
        "max_duration_hours": 8.0,
        "enable_convergence": True,
        "checkpoint_interval": 1,
    },
    "runtime": {
        "timeout_sec": 300,
        "threads": 4,
    },
    "limits": {
        "max_new_run_dirs_per_round": 100,
        "max_new_run_dirs_per_arm_call": 100,
    },
    "convergence": {
        "plateau_threshold": 0.01,
        "plateau_window": 3,
        "max_empty_rounds": 2,
    },
    "autoscaling": {
        "enabled": True,
    },
    "orchestrator": {
        "arms": []
    }
}


class ConfigScreen(Container):
    """Configuration wizard screen."""
    
    DEFAULT_CSS = """
    ConfigScreen {
        layout: horizontal;
        padding: 1;
    }
    
    .config-form {
        width: 2fr;
        padding-right: 1;
    }
    
    .config-preview {
        width: 1fr;
        border: solid $surface;
        padding: 1;
    }
    
    .section-title {
        text-style: bold;
        margin-bottom: 1;
        color: $primary;
    }
    
    .form-section {
        margin-bottom: 1;
        padding: 1;
        border: solid $surface;
    }
    
    .form-row {
        layout: horizontal;
        height: 3;
        margin-bottom: 1;
    }
    
    .form-label {
        width: 20;
        padding-top: 1;
    }
    
    .form-input {
        width: 1fr;
    }
    
    .actions-row {
        layout: horizontal;
        height: 3;
        margin-top: 1;
        dock: bottom;
    }
    
    .actions-row Button {
        margin-right: 1;
    }
    
    #preview-area {
        height: 100%;
    }
    """
    
    config: reactive[Dict[str, Any]] = reactive(dict)
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._config = DEFAULT_CONFIG.copy()
    
    def compose(self) -> ComposeResult:
        # Left side: form
        with ScrollableContainer(classes="config-form"):
            yield Static("[b]⚙️ Campaign Configuration[/]", classes="section-title")
            
            # Campaign settings
            with Collapsible(title="Campaign Settings", collapsed=False):
                with Container(classes="form-section"):
                    with Horizontal(classes="form-row"):
                        yield Static("Campaign Name:", classes="form-label")
                        yield Input(placeholder="my-campaign", id="input-name", classes="form-input")
                    
                    with Horizontal(classes="form-row"):
                        yield Static("Max Rounds:", classes="form-label")
                        yield Input(value="0", id="input-max-rounds", classes="form-input", type="integer")
                    
                    with Horizontal(classes="form-row"):
                        yield Static("Time Budget (hrs):", classes="form-label")
                        yield Input(value="8.0", id="input-duration", classes="form-input", type="number")
                    
                    with Horizontal(classes="form-row"):
                        yield Static("Convergence:", classes="form-label")
                        yield Checkbox("Enable convergence detection", id="chk-convergence", value=True)
            
            # Runtime settings
            with Collapsible(title="Runtime Settings", collapsed=True):
                with Container(classes="form-section"):
                    with Horizontal(classes="form-row"):
                        yield Static("Threads:", classes="form-label")
                        yield Input(value="4", id="input-threads", classes="form-input", type="integer")
                    
                    with Horizontal(classes="form-row"):
                        yield Static("Timeout (sec):", classes="form-label")
                        yield Input(value="300", id="input-timeout", classes="form-input", type="integer")
                    
                    with Horizontal(classes="form-row"):
                        yield Static("Evals/Round:", classes="form-label")
                        yield Input(value="100", id="input-evals", classes="form-input", type="integer")
            
            # Actions
            with Horizontal(classes="actions-row"):
                yield Button("💾 Save Config", id="btn-save", variant="primary")
                yield Button("📂 Load Config", id="btn-load", variant="default")
                yield Button("🔄 Reset", id="btn-reset", variant="warning")
        
        # Right side: preview
        with Container(classes="config-preview"):
            yield Static("[b]JSON Preview[/]", classes="section-title")
            yield TextArea(id="preview-area", read_only=True, language="json")
    
    def on_mount(self) -> None:
        """Initialize preview on mount."""
        self._update_preview()
    
    def on_input_changed(self, event: Input.Changed) -> None:
        """Update config and preview when inputs change."""
        self._update_config_from_form()
        self._update_preview()
    
    def on_checkbox_changed(self, event: Checkbox.Changed) -> None:
        """Update config and preview when checkboxes change."""
        self._update_config_from_form()
        self._update_preview()
    
    def _update_config_from_form(self):
        """Build config dict from form inputs."""
        try:
            # Campaign settings
            max_rounds = int(self.query_one("#input-max-rounds", Input).value or "0")
            duration = float(self.query_one("#input-duration", Input).value or "8.0")
            convergence = self.query_one("#chk-convergence", Checkbox).value
            
            # Runtime
            threads = int(self.query_one("#input-threads", Input).value or "4")
            timeout = int(self.query_one("#input-timeout", Input).value or "300")
            evals = int(self.query_one("#input-evals", Input).value or "100")
            
            self._config = {
                "supervisor": {
                    "max_rounds": max_rounds,
                    "max_duration_hours": duration if duration > 0 else None,
                    "enable_convergence": convergence,
                    "checkpoint_interval": 1,
                },
                "runtime": {
                    "timeout_sec": timeout,
                    "threads": threads,
                },
                "limits": {
                    "max_new_run_dirs_per_round": evals,
                    "max_new_run_dirs_per_arm_call": evals,
                },
                "convergence": {
                    "plateau_threshold": 0.01,
                    "plateau_window": 3,
                    "max_empty_rounds": 2,
                },
                "autoscaling": {"enabled": True},
                "orchestrator": {"arms": []},  # add boundary action arms; see autoresearch/README.md
            }
            
        except (ValueError, TypeError):
            pass  # Ignore conversion errors during typing
    
    def _update_preview(self):
        """Update the JSON preview."""
        try:
            preview = self.query_one("#preview-area", TextArea)
            preview.load_text(json.dumps(self._config, indent=2))
        except Exception:
            pass
    
    def get_config(self) -> Dict[str, Any]:
        """Get the current config dict."""
        return self._config.copy()
    
    def get_campaign_name(self) -> str:
        """Get the campaign name from form."""
        try:
            return self.query_one("#input-name", Input).value or "my-campaign"
        except Exception:
            return "my-campaign"
    
    async def save_config(self, path: Path) -> bool:
        """Save config to file."""
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(self._config, indent=2))
            return True
        except Exception:
            return False
    
    def load_config(self, config_or_path: Path | dict[str, Any]) -> bool:
        """Load config from file path or dict."""
        try:
            if isinstance(config_or_path, dict):
                self._config = config_or_path
            else:
                self._config = json.loads(config_or_path.read_text())
            self._update_preview()
            # TODO: Update form inputs from loaded config
            return True
        except Exception:
            return False
    
    def on_button_pressed(self, event: Button.Pressed) -> None:
        """Handle button presses."""
        if event.button.id == "btn-save":
            # Save config to file
            async def _save():
                from autoresearch.tui.widgets import FilePickerModal
                path = await self.app.push_screen(
                    FilePickerModal(
                        title="Save Config As...",
                        initial_path="autoresearch/configs",
                        filter_ext=[".json"],
                    )
                )
                if path:
                    success = await self.save_config(path)
                    if success:
                        self.app.notify(f"Config saved to {path}", title="Success")
                    else:
                        self.app.notify("Failed to save config", title="Error", severity="error")
            
            self.app.call_later(_save)
        
        elif event.button.id == "btn-load":
            # Load config from file
            async def _load():
                from autoresearch.tui.widgets import FilePickerModal
                path = await self.app.push_screen(
                    FilePickerModal(
                        title="Load Config File",
                        initial_path="autoresearch/configs",
                        filter_ext=[".json"],
                    )
                )
                if path:
                    success = self.load_config(path)
                    if success:
                        self.app.notify(f"Config loaded from {path}", title="Success")
                    else:
                        self.app.notify("Failed to load config", title="Error", severity="error")
            
            self.app.call_later(_load)
        
        elif event.button.id == "btn-reset":
            # Reset form to defaults
            self.query_one("#input-max-rounds", Input).value = "0"
            self.query_one("#input-duration", Input).value = "8.0"
            self.query_one("#chk-convergence", Checkbox).value = True
            self.query_one("#input-threads", Input).value = "4"
            self.query_one("#input-timeout", Input).value = "300"
            self.query_one("#input-evals", Input).value = "100"
            self._update_config_from_form()
            self._update_preview()
            self.app.notify("Form reset to defaults", title="Info")
