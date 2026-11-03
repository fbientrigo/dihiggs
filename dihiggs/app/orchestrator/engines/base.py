"""
engines/base.py — EngineAdapter protocol and ScanAxis enum.

Any C++ backend is represented by a class that satisfies the EngineAdapter
Protocol.  The concrete implementations are:

  M2Engine         →  DihiggsPointV2Evaluator  (second axis = M^2 in GeV^2)
  M2TrackerEngine  →  Phys_M2BandTracker       (experimental)

Axis semantics contract
-----------------------
The distinction between lambda1, M2, m12, and m12_sq is critical:

  ScanAxis.M2       — axis is M^2 = m12^2/(sin(beta)cos(beta)) (GeV^2).
                      Historical CSVs sometimes call this column "m12"
                      while storing GeV^2 values (i.e. m12_sq).
                      Never infer units from column names alone; always
                      consult the engine's axis_metadata().

The orchestrator runner records axis_metadata() in every scan_meta.json so
that downstream analysis code always has an authoritative unit annotation.
"""

from __future__ import annotations

from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Protocol, runtime_checkable

from dihiggs.app.orchestrator.grid import ScanGrid
from dihiggs.app.orchestrator.models import FixedParams


class ScanAxis(Enum):
    """
    Identifies the physical quantity along the second scan axis.

    Values
    ------
    M2:
        Axis is M^2 = m12_sq/(sin(beta)cos(beta)), in GeV^2.
        Note: some historical CSV column names say "m12" but store GeV^2.
    """

    M2 = "M2"


@runtime_checkable
class EngineAdapter(Protocol):
    """
    Protocol satisfied by each physics engine adapter.

    Concrete implementations: M2Engine, M2TrackerEngine.
    """

    @property
    def engine_name(self) -> str:
        """
        Canonical, stable identifier for this engine.

        Used as part of the grid signature to prevent resume collisions.
        Examples: ``"m2"``, ``"m2_tracker"``.
        """
        ...

    @property
    def scan_axis(self) -> ScanAxis:
        """Physical quantity on the second scan axis."""
        ...

    @property
    def executable_basename(self) -> str:
        """
        Default basename of the compiled C++ binary.

        Examples: ``"DihiggsPointV2Evaluator"``, ``"Phys_M2BandTracker"``.
        """
        ...

    def build_command(
        self,
        executable: Path,
        grid: ScanGrid,
        fixed: FixedParams,
        output_csv: Path,
    ) -> List[str]:
        """
        Build the full command list for ``subprocess.run``.

        Parameters
        ----------
        executable:
            Resolved absolute path to the C++ binary.
        grid:
            Scan grid (mphi + second axis).
        fixed:
            Fixed physics parameters for this task.
        output_csv:
            Path where the C++ binary will write its results CSV.

        Returns
        -------
        List[str]
            Command as token list, ready for ``subprocess.run``.
        """
        ...

    def axis_metadata(self) -> Dict[str, Any]:
        """
        Return authoritative metadata describing the scan axis.

        The returned dict MUST include at least:
          - ``"engine_name"``   : str
          - ``"scan_axis"``     : str (e.g. ``"M2"``)
          - ``"axis_units"``    : str (e.g. ``"dimensionless"`` or ``"GeV^2"``)
          - ``"axis_label"``    : str (e.g. ``"lambda_1"`` or ``"M2 (m12_sq)"``).
          - ``"csv_column_note"``: str explaining any historical name ambiguity.

        This dict is written into every scan_meta.json so that downstream
        analysis always has an authoritative unit annotation.
        """
        ...

    def expected_csv_columns(self) -> List[str]:
        """
        Return a list of expected CSV column names (may be partial).

        Used for basic output validation; can return ``[]`` to skip validation.
        """
        ...
