import csv
from decimal import Decimal, localcontext
import json
from pathlib import Path

CAMPAIGN = Path(__file__).resolve().parents[1] / "campaigns/t2_t4_20261008T043858Z"


def test_active_mass_table_rebuilds_quartics_and_recomputes_xstar():
    with (CAMPAIGN / "pilot_scan_tables_mh125_20.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 90
    with localcontext() as ctx:
        ctx.prec=100
        for row in rows:
            assert row["mh_GeV"] == "125.20"
            expected_z1 = Decimal(row["mh_GeV"]) ** 2 / Decimal(row["v_GeV"]) ** 2
            assert abs(Decimal(row["Z1"]) - expected_z1) < Decimal("1e-75")
        xstars = [row for row in rows if row["trajectory"] == "Xstar"]
        assert len(xstars) == 30
        assert all(abs(Decimal(row["Z7"])) < Decimal("1e-75") for row in xstars)
    detailed = [json.loads(line) for line in
        (CAMPAIGN / "active_mh125_20_constraints.jsonl").read_text().splitlines()]
    assert len(detailed) == len(rows)
    assert all(row["inputs"]["mh_GeV"] == "125.20" for row in detailed)
