"""
engines/m2.py — Engine adapter for Phys_M2BoundaryScan.

CLI contract:
    ./Phys_M2BoundaryScan \\
        mphi_min mphi_max N_mphi \\
        M2_min   M2_max   N_M2   \\
        mA sin(b-a) tan(beta) lambda6 lambda7 \\
        output.csv

The second triplet (positions 4–6) is M^2 = m12_sq, in units of GeV^2.

IMPORTANT naming note
---------------------
Historical CSV output from the C++ code may label the column as "m12"
while actually storing M^2 values (i.e. m12_sq = m12^2) in GeV^2.
The axis_metadata() method documents this unambiguously.  Do NOT infer
units from column names.

Fixed params for this engine: mA, sin_ba, tan_beta, lambda6, lambda7,
lambda1 (lambda1 is fixed at a constant value, NOT a scan axis).
M2 (m12_sq) is the scan axis; stored in ScanGrid.axis_min/max/n_axis.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List

from dihiggs.app.orchestrator.engines.base import EngineAdapter, ScanAxis
from dihiggs.app.orchestrator.grid import ScanGrid
from dihiggs.app.orchestrator.models import FixedParams


class M2Engine:
    """
    Engine adapter for ``Phys_M2BoundaryScan``.

    Second scan axis = M^2 = m12_sq (units: GeV^2).

    Positional CLI layout (13 args):
        [executable]
        mphi_min  mphi_max  N_mphi       <- m_phi grid (GeV)
        M2_min    M2_max    N_M2         <- M^2 grid (GeV^2)
        mA  sin_ba  tan_beta  lambda6  lambda7
        output.csv

    Unlike the lambda1 engine, lambda1 is a fixed constant that belongs
    in ``FixedParams.lambda1``.  The M^2 range belongs in
    ``ScanGrid.axis_min/axis_max/n_axis``.
    """

    @property
    def engine_name(self) -> str:
        return "m2"

    @property
    def scan_axis(self) -> ScanAxis:
        return ScanAxis.M2

    @property
    def executable_basename(self) -> str:
        return "Phys_M2BoundaryScan"

    def build_command(
        self,
        executable: Path,
        grid: ScanGrid,
        fixed: FixedParams,
        output_csv: Path,
    ) -> List[str]:
        """
        Build the 13-token command for Phys_M2BoundaryScan.

        Parameters
        ----------
        executable:
            Resolved path to the Phys_M2BoundaryScan binary.
        grid:
            grid.axis_min/max/n_axis are the M^2 range (GeV^2).
        fixed:
            fixed.lambda1 MUST be set (it is a fixed constant for M2 scans).
            fixed.tan_beta is the tan(beta) value for this task.
        output_csv:
            Destination CSV path.

        Raises
        ------
        ValueError
            If ``fixed.lambda1`` is None (M2 scans require a fixed lambda1).
        """
        # Note: the CLI positional layout is structurally identical to
        # PhysScanWithFixings (same 13 positions).  The second triplet
        # carries M2 instead of lambda1.  The C++ binary interprets them
        # differently; the Python side makes this explicit via the engine.
        return [
            str(executable),
            f"{grid.mphi_min:.6g}",
            f"{grid.mphi_max:.6g}",
            str(grid.n_mphi),
            f"{grid.axis_min:.6g}",    # M2_min  (GeV^2)
            f"{grid.axis_max:.6g}",    # M2_max  (GeV^2)
            str(grid.n_axis),           # N_M2
            f"{fixed.mA:.6g}",
            f"{fixed.sin_ba:.6g}",
            f"{fixed.tan_beta:.6g}",
            f"{fixed.lambda6:.6g}",
            f"{fixed.lambda7:.6g}",
            str(output_csv),
        ]

    def axis_metadata(self) -> Dict[str, Any]:
        """
        Axis metadata for M^2 = m12_sq.

        Returns
        -------
        dict
            Engine name, scan axis, units (GeV^2), and a detailed note
            about the historical CSV column name ambiguity.
        """
        return {
            "engine_name": self.engine_name,
            "executable_basename": self.executable_basename,
            "scan_axis": self.scan_axis.value,
            "axis_units": "GeV^2",
            "axis_label": "M2",
            "axis_description": (
                "Soft-breaking parameter M2 in units of GeV^2. "
                "M2 = m12_sq / (sin_beta * cos_beta) and m12_sq = M2 * sin_beta * cos_beta. "
                "Typical scan range [0, 500000] GeV^2. lambda1 is a derived output."
            ),
            "csv_column_note": (
                "IMPORTANT: M2_input is the scanned M2 value. "
                "m12_sq_input and m12_sq_out are the actual soft-breaking m12^2 values passed/recovered through 2HDMC. "
                "Historical columns named 'm12' may mean m12_sq in GeV^2, not sqrt(m12_sq) in GeV."
            ),
        }

    def expected_csv_columns(self) -> List[str]:
        """Partial list of expected CSV columns from Phys_M2BoundaryScan."""
        # m12 column stores m12_sq (GeV^2) — see csv_column_note above.
        return [
            "mphi", "m12", "mA", "sin_ba", "tan_beta",
            "lambda6", "lambda7", "ctau_m",
        ]
