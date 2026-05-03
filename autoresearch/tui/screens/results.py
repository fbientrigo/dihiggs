"""
Results screen - aggregated analysis of campaign results.
"""
from pathlib import Path
from typing import Optional
import json

from textual.app import ComposeResult
from textual.containers import Container, Horizontal, Vertical, ScrollableContainer
from textual.widgets import Static, Label, Button, DataTable, TabbedContent, TabPane
from textual.reactive import reactive
from dataclasses import dataclass


@dataclass
class ResultsSummary:
    """Summary of campaign results."""
    total_evaluations: int = 0
    unique_iterations: int = 0
    viable_points: int = 0
    positivity_pass_rate: float = 0.0
    unitarity_pass_rate: float = 0.0
    perturbativity_pass_rate: float = 0.0
    all_constraints_pass_rate: float = 0.0
    lam1_range: tuple[float, float] = (0.0, 0.0)
    mphi_range: tuple[float, float] = (0.0, 0.0)
    br_gaga_mean: float = 0.0
    br_gaga_max: float = 0.0
    loaded: bool = False
    error: Optional[str] = None


class StatCard(Static):
    """A compact stat display card."""
    
    def __init__(self, label: str, value: str = "—", **kwargs):
        super().__init__(**kwargs)
        self.label = label
        self._value = value
    
    def compose(self) -> ComposeResult:
        yield Static(f"[dim]{self.label}[/]", classes="stat-label")
        yield Static(self._value, classes="stat-value", id=f"stat-{self.id}")
    
    def update_value(self, value: str):
        self._value = value
        try:
            self.query_one(f"#stat-{self.id}", Static).update(value)
        except Exception:
            pass


class ResultsScreen(Container):
    """Results screen showing aggregated campaign analysis."""
    
    DEFAULT_CSS = """
    ResultsScreen {
        layout: vertical;
        padding: 1;
    }
    
    .section-title {
        text-style: bold;
        margin-bottom: 1;
        color: $primary;
    }
    
    .stats-grid {
        layout: grid;
        grid-size: 4;
        grid-gutter: 1;
        height: auto;
        margin-bottom: 1;
    }
    
    .stat-card {
        height: 4;
        border: solid $surface;
        padding: 0 1;
    }
    
    .stat-label {
        height: 1;
    }
    
    .stat-value {
        height: 2;
        text-style: bold;
        content-align: center middle;
    }
    
    .constraint-row {
        layout: horizontal;
        height: 3;
        margin-bottom: 1;
    }
    
    .constraint-item {
        width: 1fr;
        height: 3;
        border: solid $surface;
        padding: 0 1;
        margin-right: 1;
    }
    
    .pass-high {
        background: $success 30%;
        border: solid $success;
    }
    
    .pass-medium {
        background: $warning 30%;
        border: solid $warning;
    }
    
    .pass-low {
        background: $error 30%;
        border: solid $error;
    }
    
    .results-tabs {
        height: 1fr;
    }
    
    .results-table {
        height: 1fr;
    }
    
    .no-data {
        content-align: center middle;
        height: 100%;
        color: $text-muted;
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
    
    .plot-preview {
        height: auto;
        padding: 1;
        border: solid $surface;
    }
    
    .plot-list {
        height: auto;
        max-height: 10;
    }
    """
    
    summary: reactive[ResultsSummary] = reactive(ResultsSummary)
    campaign_dir: reactive[Optional[Path]] = reactive(None)
    
    def compose(self) -> ComposeResult:
        yield Static("[b]📈 Results Analysis[/]", classes="section-title")
        
        # Stats overview grid
        with Container(classes="stats-grid"):
            yield StatCard("Total Evaluations", "—", id="total-evals", classes="stat-card")
            yield StatCard("Iterations", "—", id="iterations", classes="stat-card")
            yield StatCard("Viable Points", "—", id="viable", classes="stat-card")
            yield StatCard("Viability Rate", "—", id="viability-rate", classes="stat-card")
        
        # Constraint pass rates
        yield Static("[b]Constraint Pass Rates[/]", classes="section-title")
        with Horizontal(classes="constraint-row"):
            yield Static("Positivity: —%", id="positivity-rate", classes="constraint-item")
            yield Static("Unitarity: —%", id="unitarity-rate", classes="constraint-item")
            yield Static("Perturbativity: —%", id="perturbativity-rate", classes="constraint-item")
            yield Static("All: —%", id="all-constraints-rate", classes="constraint-item")
        
        # Tabbed content for different views
        with TabbedContent(classes="results-tabs"):
            with TabPane("Parameters", id="tab-params"):
                yield DataTable(id="params-table", classes="results-table")
            
            with TabPane("Plots", id="tab-plots"):
                with ScrollableContainer(classes="plot-list"):
                    yield Static("📊 Available plots will appear here after loading results.", id="plots-info")
            
            with TabPane("Export", id="tab-export"):
                with Vertical():
                    yield Static("Export campaign results to various formats:", classes="section-title")
                    yield Button("📄 Export Summary CSV", id="btn-export-csv", variant="default")
                    yield Button("📊 Generate Plots", id="btn-gen-plots", variant="default")
                    yield Button("📋 Copy Summary", id="btn-copy-summary", variant="default")
        
        # Actions
        with Horizontal(classes="actions-row"):
            yield Button("🔄 Refresh", id="btn-refresh", variant="default")
            yield Button("📂 Load Campaign", id="btn-load", variant="primary")
    
    def on_mount(self) -> None:
        """Initialize tables on mount."""
        params_table = self.query_one("#params-table", DataTable)
        params_table.add_columns("Parameter", "Min", "Max", "Mean")
        params_table.add_row("λ₁", "—", "—", "—")
        params_table.add_row("m_φ", "—", "—", "—")
        params_table.add_row("tan(β)", "—", "—", "—")
        params_table.add_row("BR(h→γγ)", "—", "—", "—")
    
    def watch_summary(self, summary: ResultsSummary) -> None:
        """React to summary changes."""
        self._update_display(summary)
    
    def _update_display(self, summary: ResultsSummary):
        """Update all display elements with summary data."""
        if not summary.loaded:
            return
        
        try:
            # Update stat cards
            self.query_one("#total-evals", StatCard).update_value(f"{summary.total_evaluations:,}")
            self.query_one("#iterations", StatCard).update_value(str(summary.unique_iterations))
            self.query_one("#viable", StatCard).update_value(f"{summary.viable_points:,}")
            
            if summary.total_evaluations > 0:
                rate = summary.viable_points / summary.total_evaluations * 100
                self.query_one("#viability-rate", StatCard).update_value(f"{rate:.1f}%")
            
            # Update constraint rates with color coding
            self._update_constraint_item("#positivity-rate", "Positivity", summary.positivity_pass_rate)
            self._update_constraint_item("#unitarity-rate", "Unitarity", summary.unitarity_pass_rate)
            self._update_constraint_item("#perturbativity-rate", "Perturbativity", summary.perturbativity_pass_rate)
            self._update_constraint_item("#all-constraints-rate", "All", summary.all_constraints_pass_rate)
            
            # Update params table
            params_table = self.query_one("#params-table", DataTable)
            params_table.clear()
            params_table.add_row("λ₁", f"{summary.lam1_range[0]:.3f}", f"{summary.lam1_range[1]:.3f}", "—")
            params_table.add_row("m_φ", f"{summary.mphi_range[0]:.1f}", f"{summary.mphi_range[1]:.1f}", "—")
            params_table.add_row("BR(h→γγ)", "—", f"{summary.br_gaga_max:.6f}", f"{summary.br_gaga_mean:.6f}")
            
        except Exception:
            pass  # Widget not yet mounted
    
    def _update_constraint_item(self, widget_id: str, label: str, rate: float):
        """Update a constraint rate item with color coding."""
        try:
            widget = self.query_one(widget_id, Static)
            pct = rate * 100
            widget.update(f"{label}: {pct:.1f}%")
            
            widget.remove_class("pass-high", "pass-medium", "pass-low")
            if pct >= 80:
                widget.add_class("pass-high")
            elif pct >= 50:
                widget.add_class("pass-medium")
            else:
                widget.add_class("pass-low")
        except Exception:
            pass
    
    def load_results(self, campaign_dir: Path) -> ResultsSummary:
        """Load and aggregate results from a campaign directory."""
        self.campaign_dir = campaign_dir
        
        try:
            # Import aggregate_results functions
            from autoresearch.harness.aggregate_results import (
                find_results_files,
                aggregate_results,
                summarize_results
            )
            
            results_files = find_results_files(campaign_dir)
            if not results_files:
                return ResultsSummary(loaded=True, error="No results.csv files found")
            
            df = aggregate_results(results_files)
            if df.empty:
                return ResultsSummary(loaded=True, error="Failed to load results")
            
            raw_summary = summarize_results(df)
            
            # Convert to ResultsSummary dataclass
            summary = ResultsSummary(
                total_evaluations=raw_summary.get('total_evaluations', 0),
                unique_iterations=raw_summary.get('unique_iterations', 0),
                viable_points=raw_summary.get('viable_points', 0),
                positivity_pass_rate=raw_summary.get('positivity_pass_rate', 0.0),
                unitarity_pass_rate=raw_summary.get('unitarity_pass_rate', 0.0),
                perturbativity_pass_rate=raw_summary.get('perturbativity_pass_rate', 0.0),
                all_constraints_pass_rate=raw_summary.get('all_constraints_pass_rate', 0.0),
                lam1_range=(
                    raw_summary.get('lam1_min', 0.0),
                    raw_summary.get('lam1_max', 0.0)
                ),
                mphi_range=(
                    raw_summary.get('m_phi_min', 0.0),
                    raw_summary.get('m_phi_max', 0.0)
                ),
                br_gaga_mean=raw_summary.get('br_gaga_mean', 0.0),
                br_gaga_max=raw_summary.get('br_gaga_max', 0.0),
                loaded=True
            )
            
            self.summary = summary
            self._update_plots_list(campaign_dir)
            
            return summary
            
        except Exception as e:
            return ResultsSummary(loaded=True, error=str(e))
    
    def _update_plots_list(self, campaign_dir: Path):
        """Update the available plots list."""
        try:
            plots_dir = campaign_dir / "analysis_plots"
            plots_info = self.query_one("#plots-info", Static)
            
            if plots_dir.exists():
                plot_files = list(plots_dir.glob("*.png"))
                if plot_files:
                    lines = ["📊 Generated plots:"]
                    for pf in sorted(plot_files):
                        lines.append(f"  • {pf.name}")
                    lines.append(f"\n📁 Location: {plots_dir}")
                    plots_info.update("\n".join(lines))
                else:
                    plots_info.update("No plots generated yet. Click 'Generate Plots' to create them.")
            else:
                plots_info.update("No plots generated yet. Click 'Generate Plots' to create them.")
        except Exception:
            pass
    
    async def generate_plots(self) -> bool:
        """Generate analysis plots for the loaded campaign."""
        if self.campaign_dir is None:
            return False
        
        try:
            from autoresearch.harness.aggregate_results import (
                find_results_files,
                aggregate_results,
                plot_results
            )
            
            results_files = find_results_files(self.campaign_dir)
            df = aggregate_results(results_files)
            
            if not df.empty:
                plot_results(df, self.campaign_dir)
                self._update_plots_list(self.campaign_dir)
                return True
            return False
        except Exception:
            return False
    
    def on_button_pressed(self, event: Button.Pressed) -> None:
        """Handle button presses."""
        if event.button.id == "btn-gen-plots":
            if self.campaign_dir:
                self.run_worker(self.generate_plots())
            else:
                self.app.notify("No campaign loaded", title="Error", severity="error")
        
        elif event.button.id == "btn-load":
            # Open directory picker for campaign
            async def _browse():
                from autoresearch.tui.widgets import FilePickerModal
                path = await self.app.push_screen(
                    FilePickerModal(
                        title="Select Campaign Directory",
                        initial_path="runs",
                        select_dirs=True,
                    )
                )
                if path:
                    self.load_results(path)
                    self.app.notify(f"Campaign loaded: {path.name}", title="Results")
            
            self.app.call_later(_browse)
        
        elif event.button.id == "btn-refresh":
            # Reload results from current campaign
            if self.campaign_dir:
                self.load_results(self.campaign_dir)
                self.app.notify("Results refreshed", title="Info")
            else:
                self.app.notify("No campaign loaded", title="Warning", severity="warning")
        
        elif event.button.id == "btn-export-csv":
            # Export summary to CSV
            if not self.campaign_dir:
                self.app.notify("No campaign loaded", title="Error", severity="error")
                return
            
            async def _export():
                from autoresearch.tui.widgets import FilePickerModal
                from pathlib import Path
                import csv
                
                # Ask for save location
                path = await self.app.push_screen(
                    FilePickerModal(
                        title="Export Summary As...",
                        initial_path=str(self.campaign_dir),
                        filter_ext=[".csv"],
                    )
                )
                
                if path:
                    try:
                        # Write summary data to CSV
                        with open(path, 'w', newline='') as f:
                            writer = csv.writer(f)
                            writer.writerow(["Metric", "Value"])
                            writer.writerow(["Total Evaluations", self.summary.total_evaluations])
                            writer.writerow(["Unique Iterations", self.summary.unique_iterations])
                            writer.writerow(["Viable Points", self.summary.viable_points])
                            writer.writerow(["All Constraints Pass", f"{self.summary.all_constraints_pass_rate:.1f}%"])
                            writer.writerow(["Positivity Pass", f"{self.summary.positivity_pass_rate:.1f}%"])
                            writer.writerow(["Unitarity Pass", f"{self.summary.unitarity_pass_rate:.1f}%"])
                            writer.writerow(["Perturbativity Pass", f"{self.summary.perturbativity_pass_rate:.1f}%"])
                            writer.writerow(["Lambda1 Range", f"{self.summary.lam1_range[0]:.2f} - {self.summary.lam1_range[1]:.2f}"])
                            writer.writerow(["BR(h->gaga) Mean", f"{self.summary.br_gaga_mean:.6f}"])
                            writer.writerow(["BR(h->gaga) Max", f"{self.summary.br_gaga_max:.6f}"])
                        self.app.notify(f"Summary exported to {path.name}", title="Success")
                    except Exception as e:
                        self.app.notify(f"Export failed: {e}", title="Error", severity="error")
            
            self.app.call_later(_export)
        
        elif event.button.id == "btn-copy-summary":
            # Copy summary to clipboard (fallback: show in notification)
            if not self.campaign_dir:
                self.app.notify("No campaign loaded", title="Error", severity="error")
                return
            
            summary = f"""Campaign Summary:
Total Evaluations: {self.summary.total_evaluations:,}
Viable Points: {self.summary.viable_points:,}
All Constraints Pass: {self.summary.all_constraints_pass_rate:.1f}%
BR(h->gaga) Max: {self.summary.br_gaga_max:.6f}"""
            
            # Try to copy to clipboard (requires pyperclip or similar)
            # For now, just show in a longer notification
            self.app.notify(
                f"Summary copied:\n{summary}",
                title="Summary",
                timeout=10,
            )
