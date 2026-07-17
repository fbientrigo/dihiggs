import subprocess
import os
import json
import csv
from pathlib import Path
import pytest

def test_tracker_execution():
    if not Path("dihiggs/app/Phys_M2BandTracker").is_file():
        pytest.skip("experimental tracker is optional and not built")
    out_dir = "scripts/out"
    os.makedirs(out_dir, exist_ok=True)
    
    cmd = [
        "dihiggs/app/Phys_M2BandTracker",
        "--mphi-min=140.0",
        "--mphi-max=145.0",
        "--mphi-step=1.0",
        "--tan-beta=50.0",
        "--lam6=0.001",
        "--seed-M2-center=19500",
        "--out-points=scripts/out/points.csv",
        "--out-intervals=scripts/out/intervals.csv",
        "--out-summary=scripts/out/summary.jsonl",
        "--out-meta=scripts/out/meta.json",
        "--max-edge-iter=5"
    ]
    
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    
    # Verify outputs
    assert os.path.exists("scripts/out/points.csv")
    assert os.path.exists("scripts/out/intervals.csv")
    assert os.path.exists("scripts/out/summary.jsonl")
    assert os.path.exists("scripts/out/meta.json")
    
    with open("scripts/out/meta.json") as f:
        meta = json.load(f)
        assert "total_evaluations" in meta
        assert meta["total_evaluations"] > 0
        
    with open("scripts/out/intervals.csv") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        assert len(rows) == 6 # 140 to 145 inclusive
        assert "M2_width" in rows[0]
        
    print("SUCCESS: Tracker outputs validate against schema!")

if __name__ == "__main__":
    test_tracker_execution()
