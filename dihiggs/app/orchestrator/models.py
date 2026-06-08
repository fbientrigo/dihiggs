"""
models.py — Pure data containers for the scan orchestrator.

No physics logic; no I/O.  Only dataclasses and simple factory helpers.

Classes
-------
FixedParams
    Physics parameters that are held constant across the scan grid.
    Includes mA, sin_ba, tan_beta, lambda6, lambda7.

    NOTE: lambda1 is *not* here because for PhysScanWithFixings it is a
    scanned axis.  For Phys_M2BoundaryScan, lambda1 is a fixed scalar and
    therefore *is* included here (field ``lambda1``).  Use the appropriate
    engine's ``FixedParams`` factory helpers to build this correctly.

TaskSpec
    Everything needed to execute one C++ scan task.

TaskResult
    Record produced after a task finishes (success, failure, dry-run, skip).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional


# ---------------------------------------------------------------------------
# FixedParams
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class FixedParams:
    """
    Physics parameters that do not vary over the scan grid.

    Parameters
    ----------
    mA:
        CP-odd Higgs mass (also used as mHp), in GeV.
    sin_ba:
        sin(beta - alpha), alignment parameter.
    tan_beta:
        tan(beta), ratio of Higgs doublet VEVs.  For PhysScanWithFixings
        this is typically a *scan* dimension iterated externally; for
        Phys_M2BoundaryScan it is fixed per task.
    lambda6:
        Quartic coupling lambda_6 (dimensionless).
    lambda7:
        Quartic coupling lambda_7 (dimensionless).  Almost always 0.
    lambda1:
        Quartic coupling lambda_1.  Relevant for Phys_M2BoundaryScan only;
        for PhysScanWithFixings this is the scan axis (leave as None).

    Semantics note
    --------------
    Do NOT conflate lambda1, M2, m12, and m12_sq:
    - lambda1  is dimensionless.
    - M2 / m12_sq has units of GeV^2.
    - Historical CSVs may name a column ``m12`` while storing GeV^2 values.
    The engine adapter is the authoritative source of axis semantics.
    """

    mA: float
    sin_ba: float
    tan_beta: float
    lambda6: float
    lambda7: float
    lambda1: Optional[float] = None  # None for lambda1-engine tasks

    def as_dict(self) -> Dict[str, Any]:
        """Serialise to a plain dictionary (suitable for JSON)."""
        d: Dict[str, Any] = {
            "mA": self.mA,
            "sin_ba": self.sin_ba,
            "tan_beta": self.tan_beta,
            "lambda6": self.lambda6,
            "lambda7": self.lambda7,
        }
        if self.lambda1 is not None:
            d["lambda1"] = self.lambda1
        return d


# ---------------------------------------------------------------------------
# TaskSpec
# ---------------------------------------------------------------------------

@dataclass
class TaskSpec:
    """
    All information needed to execute (or dry-run) one scan task.

    Attributes
    ----------
    task_id:
        Unique identifier within a run (e.g. ``"tb_10000"``).
    engine_name:
        Name of the engine adapter used (e.g. ``"lambda1"`` or ``"m2"``).
    executable:
        Absolute path to the C++ binary.
    cmd:
        Full command list as it will be passed to ``subprocess.run``.
        Set by the runner after the engine builds it.
    fixed:
        Fixed physics parameters for this task.
    output_csv:
        Path where the C++ binary should write its CSV output.
    grid_sig:
        Pre-computed grid signature string.
    axis_metadata:
        Engine-specific metadata describing the scan axis semantics
        (axis name, units, engine name).
    """

    task_id: str
    engine_name: str
    executable: Path
    cmd: List[str]
    fixed: FixedParams
    output_csv: Path
    grid_sig: str
    axis_metadata: Dict[str, Any] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# TaskResult
# ---------------------------------------------------------------------------

@dataclass
class TaskResult:
    """
    Record produced after executing (or skipping/dry-running) one task.

    Attributes
    ----------
    task_id:
        Matches the originating ``TaskSpec.task_id``.
    status:
        One of: ``"done"``, ``"fail"``, ``"crash"``, ``"skip"``, ``"dry_run"``.
    grid_sig:
        Grid signature for this task.
    output_csv:
        Path of the expected CSV output.
    elapsed_s:
        Wall-clock seconds the subprocess ran (0 for skip/dry_run).
    total_attempts:
        Parsed from C++ stdout; ``None`` if not found or not run.
    triple_ok_points:
        Parsed from C++ stdout; ``None`` if not found or not run.
    returncode:
        Subprocess return code; ``None`` for skip/dry_run/crash.
    error_message:
        Human-readable failure reason (exception repr, or ``None``).
    stdout_path:
        Path to saved stdout file (only on failure).
    stderr_path:
        Path to saved stderr file (only on failure).
    extra:
        Catch-all for engine-specific or event-specific extra fields.
    """

    task_id: str
    status: str  # "done" | "fail" | "crash" | "skip" | "dry_run"
    grid_sig: str
    output_csv: Path
    elapsed_s: float = 0.0
    total_attempts: Optional[int] = None
    triple_ok_points: Optional[str] = None
    returncode: Optional[int] = None
    error_message: Optional[str] = None
    stdout_path: Optional[Path] = None
    stderr_path: Optional[Path] = None
    extra: Dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> Dict[str, Any]:
        """Serialise to a plain dictionary (suitable for JSONL / JSON)."""
        d: Dict[str, Any] = {
            "task_id": self.task_id,
            "status": self.status,
            "grid_signature": self.grid_sig,
            "output_csv": str(self.output_csv),
            "elapsed_s": self.elapsed_s,
            "total_attempts": self.total_attempts,
            "triple_ok_points": self.triple_ok_points,
            "returncode": self.returncode,
            "error_message": self.error_message,
        }
        if self.stdout_path is not None:
            d["stdout_path"] = str(self.stdout_path)
        if self.stderr_path is not None:
            d["stderr_path"] = str(self.stderr_path)
        if self.extra:
            d.update(self.extra)
        return d
