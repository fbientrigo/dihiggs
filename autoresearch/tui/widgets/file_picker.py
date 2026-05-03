"""File picker modal widget for selecting files and directories."""

from pathlib import Path
from typing import Callable, Iterable

from textual.app import ComposeResult
from textual.containers import Container, Horizontal, Vertical
from textual.screen import ModalScreen
from textual.widgets import Button, DirectoryTree, Input, Static


class FilePickerModal(ModalScreen[Path | None]):
    """Modal file picker dialog using DirectoryTree."""
    
    DEFAULT_CSS = """
    FilePickerModal {
        align: center middle;
        background: $primary 30%;
    }
    
    #picker-container {
        width: 70;
        height: 80%;
        border: thick $accent;
        background: $surface;
        padding: 1 2;
    }
    
    #picker-title {
        text-align: center;
        text-style: bold;
        margin-bottom: 1;
    }
    
    #path-input {
        margin-bottom: 1;
    }
    
    FilePickerModal DirectoryTree {
        height: 1fr;
        border: solid $primary;
        margin-bottom: 1;
    }
    
    #picker-buttons {
        height: auto;
        align: right middle;
    }
    
    #picker-buttons Button {
        margin-left: 1;
    }
    """
    
    def __init__(
        self,
        title: str = "Select File",
        initial_path: str = ".",
        filter_ext: Iterable[str] | None = None,
        select_dirs: bool = False,
    ) -> None:
        """Initialize file picker.
        
        Args:
            title: Dialog title
            initial_path: Starting directory path
            filter_ext: File extensions to show (e.g., [".json", ".yaml"])
            select_dirs: If True, allow directory selection; if False, only files
        """
        super().__init__()
        self.title = title
        self.initial_path = Path(initial_path).expanduser().resolve()
        self.filter_ext = set(filter_ext) if filter_ext else None
        self.select_dirs = select_dirs
        self._selected_path: Path | None = None
        self._updating_from_tree: bool = False  # Guard against circular updates
    
    def compose(self) -> ComposeResult:
        """Compose the file picker UI."""
        with Vertical(id="picker-container"):
            yield Static(f"📂 {self.title}", id="picker-title")
            yield Input(
                str(self.initial_path),
                id="path-input",
                placeholder="Enter path or select from tree...",
            )
            yield DirectoryTree(str(self.initial_path), id="tree")
            with Horizontal(id="picker-buttons"):
                yield Button("Cancel", id="btn-cancel")
                yield Button("Select", id="btn-select", variant="primary")
    
    def on_mount(self) -> None:
        """Focus tree on mount."""
        self.query_one("#tree", DirectoryTree).focus()
    
    def on_directory_tree_file_selected(self, event: DirectoryTree.FileSelected) -> None:
        """Handle file selection from tree."""
        path = event.path
        
        # Apply extension filter if specified
        if self.filter_ext and not any(str(path).endswith(ext) for ext in self.filter_ext):
            self.notify(
                f"Please select a file with extension: {', '.join(self.filter_ext)}",
                severity="warning",
            )
            return
        
        self._selected_path = path
        self.query_one("#path-input", Input).value = str(path)
    
    def on_directory_tree_directory_selected(self, event: DirectoryTree.DirectorySelected) -> None:
        """Handle directory selection from tree."""
        # Update the input to show current directory path (with guard)
        self._updating_from_tree = True
        self.query_one("#path-input", Input).value = str(event.path)
        self._updating_from_tree = False
        # For directory selection mode, track it as a potential selection
        if self.select_dirs:
            self._selected_path = event.path
    
    def on_input_changed(self, event: Input.Changed) -> None:
        """Handle manual path input."""
        # Skip if this change came from tree selection (avoid circular update)
        if self._updating_from_tree:
            return
        if event.input.id == "path-input":
            path = Path(event.value).expanduser()
            if path.exists() and path.is_dir():
                tree = self.query_one("#tree", DirectoryTree)
                tree.path = str(path)
    
    def on_button_pressed(self, event: Button.Pressed) -> None:
        """Handle button presses."""
        if event.button.id == "btn-select":
            # Get the path from input field (most reliable source)
            input_val = self.query_one("#path-input", Input).value.strip()
            if not input_val:
                self.notify("No path selected", severity="warning")
                return
            
            path = Path(input_val).expanduser().resolve()
            
            if not path.exists():
                self.notify("Path does not exist", severity="error")
                return
            
            # Validate selection type
            if self.select_dirs:
                if not path.is_dir():
                    self.notify("Please select a directory", severity="warning")
                    return
                self.dismiss(path)
            else:
                if not path.is_file():
                    self.notify("Please select a file", severity="warning")
                    return
                
                # Apply extension filter for files
                if self.filter_ext:
                    if not any(str(path).endswith(ext) for ext in self.filter_ext):
                        self.notify(
                            f"Please select a file with extension: {', '.join(self.filter_ext)}",
                            severity="warning",
                        )
                        return
                
                self.dismiss(path)
        
        elif event.button.id == "btn-cancel":
            self.dismiss(None)
