import subprocess
import csv
import os
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


def test_disconnected_intervals_and_predicted_branch_selection(tmp_path):
    source = tmp_path / "tracker_interval_check.cpp"
    source.write_text(r'''
#include "M2IntervalDetector.hpp"
#include <cassert>
#include <vector>

int main() {
    std::vector<PointResult> points(8);
    for (int i = 0; i < 8; ++i) {
        points[i].m_phi = 200.0;
        points[i].M2_input = 100.0 * i;
    }
    points[1].triple_ok = points[2].triple_ok = true;
    points[5].triple_ok = points[6].triple_ok = true;
    const auto intervals = detect_intervals(points);
    assert(intervals.size() == 2);
    assert(intervals[0].M2_low == 100.0 && intervals[0].M2_high == 200.0);
    assert(intervals[1].M2_outer_low == 400.0 && intervals[1].M2_outer_high == 700.0);
    assert(select_interval_nearest(intervals, 575.0).M2_center == 550.0);
    assert(select_interval_nearest(intervals, 125.0).M2_center == 150.0);
}
''')
    executable = tmp_path / "tracker_interval_check"
    subprocess.run([
        "g++", "-std=c++17", f"-I{ROOT / 'dihiggs/include'}", f"-I{ROOT / '2hdmc/src'}",
        str(source), str(ROOT / "dihiggs/src/M2IntervalDetector.cpp"), "-o", str(executable),
    ], check=True)
    subprocess.run([str(executable)], check=True)


def test_refined_tracker_bounds_match_dense_canonical_grid_in_pilot_domain(tmp_path):
    tracker = ROOT / "dihiggs/app/Phys_M2BandTracker"
    producer = ROOT / "dihiggs/app/DihiggsPointV2Evaluator"
    if not tracker.is_file() or not producer.is_file():
        pytest.skip("build tracker and canonical producer first")
    points, intervals = tmp_path / "points.csv", tmp_path / "intervals.csv"
    env = {**os.environ, "OMP_NUM_THREADS": "1"}
    common = [
        "--mh=125.13", "--mphi-min=130", "--mphi-max=130", "--mphi-step=1",
        "--ma=300", "--mhp=300", "--sin-ba=0.999", "--tan-beta=50",
        "--lam6=0.1", "--lam7=0", "--m2-min=15000", "--m2-max=18000",
        "--seed-M2-center=16500", "--seed-M2-halfwidth=1500", "--seed-n-M2=61",
        "--edge-tol=0.1", "--max-edge-iter=20", f"--out-points={points}",
        f"--out-intervals={intervals}", f"--out-summary={tmp_path / 'summary.jsonl'}",
        f"--out-meta={tmp_path / 'meta.json'}",
    ]
    subprocess.run([str(tracker), *common], check=True, env=env, capture_output=True, text=True)
    with intervals.open(newline="") as handle:
        tracked = next(csv.DictReader(handle))

    dense = tmp_path / "dense.csv"
    subprocess.run([
        str(producer), "--campaign-id", "tracker-test", "--run-id", "dense", "--mh", "125.13",
        "--mH-min", "130", "--mH-max", "130", "--n-mH", "1", "--mA", "300",
        "--mHp", "300", "--yukawa-type", "1", "--sin-ba", "0.999", "--tan-beta", "50",
        "--M2-min", "15000", "--M2-max", "18000", "--n-M2", "301", "--lambda6", "0.1",
        "--lambda7", "0", "--output", str(dense),
    ], check=True, env=env, capture_output=True, text=True)
    with dense.open(newline="") as handle:
        accepted = [float(row["M2_input_GeV2"]) for row in csv.DictReader(handle)
                    if float(row["theory_ok_v1"]) == 1.0]
    assert float(tracked["M2_low"]) == pytest.approx(min(accepted), abs=10.1)
    assert float(tracked["M2_high"]) == pytest.approx(max(accepted), abs=10.1)
    assert float(tracked["M2_low"]) <= float(tracked["M2_high"])


def test_tracker_local_dense_fallback_recovers_pilot_band(tmp_path):
    tracker = ROOT / "dihiggs/app/Phys_M2BandTracker"
    if not tracker.is_file():
        pytest.skip("build tracker first")
    completed = subprocess.run([
        str(tracker), "--mh=125.13", "--mphi-min=130", "--mphi-max=130", "--mphi-step=1",
        "--ma=300", "--mhp=300", "--sin-ba=0.999", "--tan-beta=50", "--lam6=0.1",
        "--lam7=0", "--m2-min=15000", "--m2-max=18000", "--seed-M2-center=15050",
        "--seed-M2-halfwidth=50", "--seed-n-M2=5", "--fallback-dense-count=301",
        "--fallback-dense-pad=2000", "--allow-mass-step-halving=false",
        f"--out-points={tmp_path / 'points.csv'}", f"--out-intervals={tmp_path / 'intervals.csv'}",
        f"--out-summary={tmp_path / 'summary.jsonl'}", f"--out-meta={tmp_path / 'meta.json'}",
    ], check=True, env={**os.environ, "OMP_NUM_THREADS": "1"}, capture_output=True, text=True)
    assert "Local Dense Scan succeeded" in completed.stdout
    assert len((tmp_path / "intervals.csv").read_text().splitlines()) == 2
