"""
engines/lambda1.py — Engine adapter for PhysScanWithFixings.

CLI contract:
    ./PhysScanWithFixings \\
        mphi_min mphi_max N_mphi \\
        lam1_min lam1_max N_lam1 \\
        mA sin(b-a) tan(beta) lambda6 lambda7 \\
        output.csv

The second triplet (positions 4–6) is lambda_1 (dimensionless quartic
coupling).  This is NOT M2, NOT m12, NOT m12_sq.

Fixed params for this engine:  mA, sin_ba, tan_beta, lambda6, lambda7.
Lambda1 is the scan axis; it is stored in ScanGrid.axis_min/max/n_axis.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List

from dihiggs.app.orchestrator.engines.base import EngineAdapter, ScanAxis
from dihiggs.app.orchestrator.grid import ScanGrid
from dihiggs.app.orchestrator.models import FixedParams


class Lambda1Engine:
    """
    Engine adapter for ``PhysScanWithFixings``.

    Second scan axis = lambda_1 (dimensionless).

    Positional CLI layout (13 args):
        [executable]
        mphi_min  mphi_max  N_mphi       <- m_phi grid  (GeV)
        lam1_min  lam1_max  N_lam1       <- lambda_1 grid (dimensionless)
        mA  sin_ba  tan_beta  lambda6  lambda7
        output.csv
    """

    # Satisfy EngineAdapter protocol as a concrete class (not a runtime
    # protocol user, so we don't inherit — duck typing suffices).

    @property
    def engine_name(self) -> str:
        return "lambda1"

    @property
    def scan_axis(self) -> ScanAxis:
        return ScanAxis.LAMBDA1

    @property
    def executable_basename(self) -> str:
        return "PhysScanWithFixings"

    def build_command(
        self,
        executable: Path,
        grid: ScanGrid,
        fixed: FixedParams,
        output_csv: Path,
    ) -> List[str]:
        """
        Build the 13-token command for PhysScanWithFixings.

        Parameters
        ----------
        executable:
            Resolved path to the PhysScanWithFixings binary.
        grid:
            grid.axis_min/max/n_axis are the lambda_1 range.
        fixed:
            fixed.tan_beta is the tan(beta) value for this task.
            fixed.lambda1 MUST be None (lambda1 is the scan axis here).
        output_csv:
            Destination CSV path.

        Raises
        ------
        ValueError
            If ``fixed.lambda1`` is not None (would indicate a caller
            error where lambda1 has been placed in both fixed and grid).
        """
        if fixed.lambda1 is not None:
            raise ValueError(
                "Lambda1Engine: fixed.lambda1 must be None because lambda1 "
                "is the scan axis.  Pass lambda1 via ScanGrid.axis_min/max/n_axis."
            )

        return [
            str(executable),
            f"{grid.mphi_min:.6g}",
            f"{grid.mphi_max:.6g}",
            str(grid.n_mphi),
            f"{grid.axis_min:.6g}",    # lam1_min
            f"{grid.axis_max:.6g}",    # lam1_max
            str(grid.n_axis),           # N_lam1
            f"{fixed.mA:.6g}",
            f"{fixed.sin_ba:.6g}",
            f"{fixed.tan_beta:.6g}",
            f"{fixed.lambda6:.6g}",
            f"{fixed.lambda7:.6g}",
            str(output_csv),
        ]

    def axis_metadata(self) -> Dict[str, Any]:
        """
        Axis metadata for lambda_1.

        Returns
        -------
        dict
            Engine name, scan axis, units, label, and a note explaining
            that there is no m12/M2 ambiguity for this engine.
        """
        return {
            "engine_name": self.engine_name,
            "executable_basename": self.executable_basename,
            "scan_axis": self.scan_axis.value,
            "axis_units": "dimensionless",
            "axis_label": "lambda_1",
            "axis_description": (
                "Quartic coupling lambda_1 in the 2HDM scalar potential. "
                "Dimensionless; typical scan range [0, 12]."
            ),
            "csv_column_note": (
                "This engine does NOT produce M2 or m12_sq columns. "
                "Do not conflate lambda1 with M2/m12/m12_sq."
            ),
        }

    def expected_csv_columns(self) -> List[str]:
        """Partial list of expected CSV columns from PhysScanWithFixings."""
        # These are the first few columns; the full schema has ~29 columns.
        return [
            "mphi", "lam1", "mA", "sin_ba", "tan_beta",
            "lambda6", "lambda7", "ctau_m",
        ]
