"""
dihiggs.app.orchestrator
========================
Modular scan harness supporting PhysScanWithFixings (lambda1 axis)
and Phys_M2BoundaryScan (M2/m12_sq axis).

Public re-exports
-----------------
- ScanRunner          : main execution class
- ScanGrid            : scan-grid descriptor
- FixedParams         : fixed physics parameters
- TaskResult          : per-task result record
- grid_signature      : stable grid hash
- Lambda1Engine       : adapter for PhysScanWithFixings
- M2Engine            : adapter for Phys_M2BoundaryScan
- ScanAxis            : enum distinguishing lambda1 vs M2 axis semantics
"""

from dihiggs.app.orchestrator.engines.base import EngineAdapter, ScanAxis
from dihiggs.app.orchestrator.engines.lambda1 import Lambda1Engine
from dihiggs.app.orchestrator.engines.m2 import M2Engine
from dihiggs.app.orchestrator.grid import ScanGrid, grid_signature
from dihiggs.app.orchestrator.models import FixedParams, TaskResult, TaskSpec
from dihiggs.app.orchestrator.runner import ScanRunner

__all__ = [
    "ScanRunner",
    "ScanGrid",
    "FixedParams",
    "TaskSpec",
    "TaskResult",
    "grid_signature",
    "Lambda1Engine",
    "M2Engine",
    "ScanAxis",
    "EngineAdapter",
]
