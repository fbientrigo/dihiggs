"""
engines/gen_fixings.py — Engine adapter for GenScanWithFixings ("Stage 2").

CLI contract:
    ./GenScanWithFixings \
        --bronze-csv=<shard.csv> \
        --output-csv=<validated.csv> \
        [--calibration-n=50] \
        [--calibration-frac=0.10] \
        [--rng-seed=0]

Unlike the other engines, this one does not generate its own (m_phi, axis)
grid: it consumes a pre-generated chris/CalcLambda1ScanFixings bronze shard
CSV (``fixed.bronze_shard_csv``) and calibrates each row against 2HDMC. The
``ScanGrid`` passed to ``build_command`` is a placeholder required by
``ScanRunner``'s interface; its ``mphi_*``/``axis_*`` fields are ignored.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List

from dihiggs.app.orchestrator.engines.base import EngineAdapter, ScanAxis
from dihiggs.app.orchestrator.grid import ScanGrid
from dihiggs.app.orchestrator.models import FixedParams


class GenFixingsEngine:
    """
    Engine adapter for ``GenScanWithFixings`` (Stage 2 calibration/validation).

    Second scan axis = N/A (``ScanAxis.GEN_FIXINGS``); this engine consumes
    pre-generated bronze rows rather than generating its own grid.
    """

    @property
    def engine_name(self) -> str:
        return "gen_fixings"

    @property
    def scan_axis(self) -> ScanAxis:
        return ScanAxis.GEN_FIXINGS

    @property
    def executable_basename(self) -> str:
        return "GenScanWithFixings"

    def build_command(
        self,
        executable: Path,
        grid: ScanGrid,
        fixed: FixedParams,
        output_csv: Path,
    ) -> List[str]:
        """
        Build the command for GenScanWithFixings using named flags.

        ``grid`` is unused (placeholder required by the EngineAdapter
        protocol/ScanRunner); the actual scan domain is the bronze shard
        referenced by ``fixed.bronze_shard_csv``.
        """
        if not fixed.bronze_shard_csv:
            raise ValueError(
                "GenFixingsEngine requires fixed.bronze_shard_csv to be set."
            )

        cmd = [
            str(executable),
            f"--bronze-csv={fixed.bronze_shard_csv}",
            f"--output-csv={output_csv}",
        ]
        if fixed.calibration_n is not None:
            cmd.append(f"--calibration-n={fixed.calibration_n}")
        if fixed.calibration_frac is not None:
            cmd.append(f"--calibration-frac={fixed.calibration_frac:.6g}")
        if fixed.rng_seed is not None:
            cmd.append(f"--rng-seed={fixed.rng_seed}")
        return cmd

    def axis_metadata(self) -> Dict[str, Any]:
        """
        Axis metadata for GenScanWithFixings.

        This engine has no generated scan axis: it validates rows already
        produced by chris/CalcLambda1ScanFixings (the bronze shard).
        """
        return {
            "engine_name": self.engine_name,
            "executable_basename": self.executable_basename,
            "scan_axis": self.scan_axis.value,
            "axis_units": "n/a",
            "axis_label": "gen_fixings (bronze-shard calibration)",
            "axis_description": (
                "No generated grid: GenScanWithFixings reads pre-generated "
                "bronze rows from fixed.bronze_shard_csv (chris/"
                "CalcLambda1ScanFixings output), calibrates each row's "
                "reconstructed generic-basis lambdas against 2HDMC via a "
                "+/-10% N=50 random search, and emits validated points "
                "('silver' layer)."
            ),
            "csv_column_note": (
                "Output is the 29-column legacy PhysScanWithFixings schema "
                "plus calibration diagnostics, stability_ok, chris "
                "cross-check columns, and replay-safety metadata."
            ),
        }

    def expected_csv_columns(self) -> List[str]:
        """Expected CSV columns from GenScanWithFixings."""
        return [
            "m_phi", "mA", "alpha", "beta", "lambda6", "lambda7", "m12",
            "sin_ba", "tan_beta",
            "positivity_ok", "unitarity_ok", "perturbativity_ok",
            "width_bb", "width_tautau", "width_WW", "width_ZZ",
            "width_gaga", "width_Zga", "width_gg", "width_hh",
            "total_width", "br_gaga",
            "lam1", "computed_lam1", "lam2", "computed_lam2",
            "lam3", "lam4", "lam5",
            "mA_target", "mH_target", "mh_calibrated", "mHp_calibrated",
            "sba_calibrated", "calibration_score", "calibration_n_used",
            "variation_idx",
            "stability_ok",
            "chris_width_bb", "chris_width_tautau", "chris_width_gg",
            "chris_width_gaga", "chris_width_Zga", "chris_ctau_mm",
            "delta_width_bb", "ratio_width_bb",
            "delta_width_tautau", "ratio_width_tautau",
            "delta_width_gg", "ratio_width_gg",
            "delta_width_gaga", "ratio_width_gaga",
            "delta_width_Zga", "ratio_width_Zga",
            "delta_ctau_mm", "ratio_ctau_mm",
            "m12_2_used", "m12_2_gen_after_set", "delta_m12_2_gen_minus_used",
            "replay_semantics_version", "yukawa_type", "higgs_state",
            "model_api_path", "calc_engine", "git_commit", "git_dirty",
        ]
