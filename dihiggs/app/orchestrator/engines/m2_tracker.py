"""
engines/m2_tracker.py — Engine adapter for Phys_M2BandTracker.

CLI contract:
    ./Phys_M2BandTracker \
        --mphi-min=... \
        --mphi-max=... \
        --mphi-step=... \
        --m2-min=... \
        --m2-max=... \
        --ma=... \
        --sin-ba=... \
        --tan-beta=... \
        --lam6=... \
        --lam7=... \
        --out-points=output.csv \
        --out-intervals=... \
        --out-summary=... \
        --out-meta=...

Unlike legacy engines, it uses long-form named flags.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List

from dihiggs.app.orchestrator.engines.base import EngineAdapter, ScanAxis
from dihiggs.app.orchestrator.grid import ScanGrid
from dihiggs.app.orchestrator.models import FixedParams


class M2TrackerEngine:
    """
    Engine adapter for ``Phys_M2BandTracker``.

    Second scan axis = M^2 = m12_sq/(sin(beta)cos(beta)) (GeV^2).
    """

    @property
    def engine_name(self) -> str:
        return "m2_tracker"

    @property
    def scan_axis(self) -> ScanAxis:
        return ScanAxis.M2

    @property
    def executable_basename(self) -> str:
        return "Phys_M2BandTracker"

    def build_command(
        self,
        executable: Path,
        grid: ScanGrid,
        fixed: FixedParams,
        output_csv: Path,
    ) -> List[str]:
        """
        Build the command for Phys_M2BandTracker using named flags.
        """
        # Determine mphi_step
        # If n_mphi > 1, step = (max - min) / (n_mphi - 1), else just fallback to a default like 1.0
        if grid.n_mphi > 1:
            step_mphi = (grid.mphi_max - grid.mphi_min) / (grid.n_mphi - 1)
        else:
            step_mphi = 1.0

        # Define aux output paths next to the main points CSV
        out_dir = output_csv.parent
        base_name = output_csv.stem
        
        out_intervals = out_dir / f"{base_name}_intervals.csv"
        out_summary = out_dir / f"{base_name}_summary.jsonl"
        out_meta = out_dir / f"{base_name}_meta.json"

        # The seed center is historically taken near the grid minimum or some safe point.
        # For this engine, we allow the grid.axis_min and axis_max to act as the global fallback bounds,
        # but we must pick a seed center. We can default it to the middle of the global bounds.
        seed_m2_center = (grid.axis_min + grid.axis_max) / 2.0
        seed_m2_halfwidth = (grid.axis_max - grid.axis_min) / 2.0

        return [
            str(executable),
            f"--mphi-min={grid.mphi_min:.6g}",
            f"--mphi-max={grid.mphi_max:.6g}",
            f"--mphi-step={step_mphi:.6g}",
            f"--m2-min={grid.axis_min:.6g}",
            f"--m2-max={grid.axis_max:.6g}",
            f"--ma={fixed.mA:.6g}",
            f"--mh={(125.20 if fixed.mh is None else fixed.mh):.17g}",
            f"--mhp={(fixed.mA if fixed.mHp is None else fixed.mHp):.17g}",
            f"--sin-ba={fixed.sin_ba:.6g}",
            f"--tan-beta={fixed.tan_beta:.6g}",
            f"--lam6={fixed.lambda6:.6g}",
            f"--lam7={fixed.lambda7:.6g}",
            f"--seed-M2-center={seed_m2_center:.6g}",
            f"--seed-M2-halfwidth={seed_m2_halfwidth:.6g}",
            f"--out-points={output_csv}",
            f"--out-intervals={out_intervals}",
            f"--out-summary={out_summary}",
            f"--out-meta={out_meta}"
        ]

    def axis_metadata(self) -> Dict[str, Any]:
        """
        Axis metadata for M^2 via Phys_M2BandTracker.
        """
        return {
            "engine_name": self.engine_name,
            "executable_basename": self.executable_basename,
            "scan_axis": self.scan_axis.value,
            "axis_units": "GeV^2",
            "axis_label": "M2",
            "axis_description": (
                "Soft-breaking parameter M2 dynamically tracked via Phys_M2BandTracker. "
                "M2 = m12_sq/(sin(beta)*cos(beta)). "
                "Output is non-canonical; boundaries are validated only inside tested pilot domains."
            ),
            "csv_column_note": (
                "Generates points.csv, intervals.csv, summary.jsonl, and meta.json. "
                "The main output_csv acts as the points file."
            ),
        }

    def expected_csv_columns(self) -> List[str]:
        """Partial list of expected CSV columns from Phys_M2BandTracker."""
        return [
            "m_phi", "M2_input", "m12_sq_input", "m12_sq_out",
            "lam1_out", "lam2_out", "theory_ok", "triple_ok",
        ]
