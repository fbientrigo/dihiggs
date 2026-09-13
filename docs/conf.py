from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

project = "DiHiggs 2HDMC Evaluation Core"
author = "fbientrigo"
html_title = "DiHiggs Documentation"
root_doc = "index"

extensions = [
    "myst_parser",
    "sphinx.ext.intersphinx",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

myst_heading_anchors = 3
myst_enable_extensions = ["dollarmath", "amsmath"]
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "**/.ipynb_checkpoints",
]
suppress_warnings = [
    "myst.xref_missing",
    "toc.not_readable",
]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
}

html_theme = "sphinx_rtd_theme"
html_static_path = []
