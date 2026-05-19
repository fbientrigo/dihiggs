import csv
from pathlib import Path

from scripts.summarize_christopher_fixed_lam1_campaign import summarize_csv


def _write_csv(path: Path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def test_summary_handles_tiny_fake_scan_csv(tmp_path):
    csv_path = tmp_path / "scan_tb_10000.csv"
    fields = [
        "m_phi",
        "lam1",
        "lambda6",
        "tan_beta",
        "positivity_ok",
        "unitarity_ok",
        "perturbativity_ok",
        "total_width",
        "br_gaga",
    ]
    rows = [
        {
            "m_phi": "125.0",
            "lam1": "0.4",
            "lambda6": "0.01",
            "tan_beta": "10000",
            "positivity_ok": "1",
            "unitarity_ok": "1",
            "perturbativity_ok": "1",
            "total_width": "1.0",
            "br_gaga": "0.001",
        },
        {
            "m_phi": "130.0",
            "lam1": "0.4",
            "lambda6": "0.01",
            "tan_beta": "10000",
            "positivity_ok": "1",
            "unitarity_ok": "0",
            "perturbativity_ok": "1",
            "total_width": "2.0",
            "br_gaga": "0.002",
        },
    ]
    _write_csv(csv_path, fields, rows)
    out = summarize_csv(csv_path, expected_n_mphi=400)
    assert out["status"] == "LOW_YIELD"
    assert out["n_rows"] == 2
    assert out["n_triple_ok"] == 1


def test_summary_detects_missing_required_columns(tmp_path):
    csv_path = tmp_path / "scan_tb_10000.csv"
    fields = ["m_phi", "lam1", "lambda6", "tan_beta"]
    rows = [{"m_phi": "125", "lam1": "0.4", "lambda6": "0.01", "tan_beta": "10000"}]
    _write_csv(csv_path, fields, rows)
    out = summarize_csv(csv_path, expected_n_mphi=400)
    assert out["status"] == "MISSING_COLUMNS"


def test_summary_detects_constants(tmp_path):
    csv_path = tmp_path / "scan_tb_10000.csv"
    fields = [
        "m_phi",
        "lam1",
        "lambda6",
        "tan_beta",
        "positivity_ok",
        "unitarity_ok",
        "perturbativity_ok",
        "total_width",
        "br_gaga",
    ]
    rows = []
    for i in range(5):
        rows.append(
            {
                "m_phi": str(125.0 + i),
                "lam1": "1.0",
                "lambda6": "0.1",
                "tan_beta": "10000",
                "positivity_ok": "1",
                "unitarity_ok": "1",
                "perturbativity_ok": "1",
                "total_width": "1.0",
                "br_gaga": "0.001",
            }
        )
    _write_csv(csv_path, fields, rows)
    out = summarize_csv(csv_path, expected_n_mphi=5)
    assert out["status"] == "OK"
    assert out["lambda1"] == 1.0
    assert out["lambda6"] == 0.1
    assert out["tan_beta"] == 10000.0
