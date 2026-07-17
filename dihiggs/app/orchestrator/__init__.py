"""
dihiggs.app.orchestrator
========================
Modular scan harness with canonical Lambda1EvaluatorV2 and
DihiggsPointV2Evaluator paths plus explicit legacy/experimental adapters.

Public re-exports
-----------------
- ScanRunner          : main execution class
- ScanGrid            : scan-grid descriptor
- FixedParams         : fixed physics parameters
- TaskResult          : per-task result record
- grid_signature      : stable grid hash
- Lambda1Engine       : compatibility adapter for PhysScanWithFixings
- M2Engine            : canonical DihiggsPointV2Evaluator adapter
- run_lambda1_v2      : canonical input/manifest runner
- ScanAxis            : enum distinguishing lambda1 vs M2 axis semantics
"""

from dihiggs.app.orchestrator.engines.base import EngineAdapter, ScanAxis
from dihiggs.app.orchestrator.engines.lambda1 import Lambda1Engine
from dihiggs.app.orchestrator.engines.m2 import M2Engine
from dihiggs.app.orchestrator.lambda1_v2 import run_lambda1_v2
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
    "run_lambda1_v2",
    "ScanAxis",
    "EngineAdapter",
]
