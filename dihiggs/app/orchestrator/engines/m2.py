"""Canonical ``dihiggs.point.v2`` M² producer adapter."""

from pathlib import Path
from typing import Any, Dict, List

from dihiggs.app.orchestrator.engines.base import ScanAxis
from dihiggs.app.orchestrator.grid import ScanGrid
from dihiggs.app.orchestrator.models import FixedParams


class M2Engine:
    def __init__(self, campaign_id: str = "campaign", run_id: str = "run") -> None:
        self.campaign_id = campaign_id
        self.run_id = run_id

    def set_provenance(self, campaign_id: str, run_id: str) -> None:
        self.campaign_id, self.run_id = campaign_id, run_id

    @property
    def engine_name(self) -> str:
        return "m2"

    @property
    def scan_axis(self) -> ScanAxis:
        return ScanAxis.M2

    @property
    def executable_basename(self) -> str:
        return "DihiggsPointV2Evaluator"

    def build_command(
        self, executable: Path, grid: ScanGrid, fixed: FixedParams, output_csv: Path
    ) -> List[str]:
        value = lambda number: format(number, ".17g")
        mh = 125.13 if fixed.mh is None else fixed.mh
        mHp = fixed.mA if fixed.mHp is None else fixed.mHp
        yukawa_type = 1 if fixed.yukawa_type is None else fixed.yukawa_type
        if yukawa_type not in (1, 2, 3, 4):
            raise ValueError("yukawa_type must be one of 1, 2, 3, 4")
        return [
            str(executable),
            "--campaign-id", self.campaign_id,
            "--run-id", self.run_id,
            "--mh", value(mh),
            "--mH-min", value(grid.mphi_min),
            "--mH-max", value(grid.mphi_max),
            "--n-mH", str(grid.n_mphi),
            "--mA", value(fixed.mA),
            "--mHp", value(mHp),
            "--yukawa-type", str(yukawa_type),
            "--sin-ba", value(fixed.sin_ba),
            "--tan-beta", value(fixed.tan_beta),
            "--M2-min", value(grid.axis_min),
            "--M2-max", value(grid.axis_max),
            "--n-M2", str(grid.n_axis),
            "--lambda6", value(fixed.lambda6),
            "--lambda7", value(fixed.lambda7),
            "--output", str(output_csv),
        ]

    def axis_metadata(self) -> Dict[str, Any]:
        return {
            "engine_name": self.engine_name,
            "executable_basename": self.executable_basename,
            "schema_version": "dihiggs.point.v2",
            "scan_axis": self.scan_axis.value,
            "axis_units": "GeV^2",
            "axis_label": "M2",
            "axis_description": "M2 = m12_sq / (sin(beta) * cos(beta)).",
            "mass_convention": {
                "mh_GeV": 125.13,
                "source": "PDG 2026 Higgs listing",
                "source_url": "https://pdg.lbl.gov/encoder_listings/s126.pdf",
            },
            "acceptance": {
                "construction_ok": "exact THDM::set_param_phys return value",
                "numerical_ok": "finite reconstructed lambda1-lambda7, tan_beta, m12_sq and M2",
                "triple_ok_legacy": "positivity_reported_ok && unitarity_ok && perturbativity_ok",
                "theory_ok_v1": "triple_ok_legacy",
                "experimental_ok": "unevaluated (nan)",
            },
            "csv_column_note": "M2 and m12_sq are distinct; m12_sq=M2*sin(beta)*cos(beta).",
        }

    def expected_csv_columns(self) -> List[str]:
        return [
            "schema_version", "campaign_id", "run_id", "point_id", "mh_input_GeV",
            "mH_input_GeV", "M2_input_GeV2", "m12_sq_input_GeV2",
            "construction_ok", "numerical_ok", "triple_ok_legacy", "theory_ok_v1",
            "experimental_evaluated", "experimental_ok", "total_width_GeV", "ctau_mm",
        ]
