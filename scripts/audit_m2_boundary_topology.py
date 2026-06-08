import csv
import sys
import argparse
from pathlib import Path
from collections import defaultdict

def main():
    parser = argparse.ArgumentParser(description="Audit M2 boundary topology.")
    parser.add_argument("--input", required=True, help="Input CSV file path.")
    parser.add_argument("--output-dir", required=True, help="Output directory path.")
    args = parser.parse_args()

    input_path = Path(args.input)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Dictionary to hold data: m_phi -> list of (M2_input, triple_ok)
    data = defaultdict(list)

    with open(input_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            m_phi = float(row['m_phi'])
            M2_input = float(row['M2_input'])
            pos_ok = float(row['positivity_ok']) == 1.0
            uni_ok = float(row['unitarity_ok']) == 1.0
            pert_ok = float(row['perturbativity_ok']) == 1.0
            triple_ok = pos_ok and uni_ok and pert_ok
            data[m_phi].append((M2_input, triple_ok))

    summary = []
    
    for m_phi in sorted(data.keys()):
        points = data[m_phi]
        points.sort(key=lambda x: x[0])  # Sort by M2_input

        intervals = []
        in_interval = False
        start_M2 = None
        last_M2 = None

        for m2, ok in points:
            if ok:
                if not in_interval:
                    in_interval = True
                    start_M2 = m2
                last_M2 = m2
            else:
                if in_interval:
                    intervals.append((start_M2, last_M2))
                    in_interval = False

        if in_interval:
            intervals.append((start_M2, last_M2))

        n_intervals = len(intervals)
        if n_intervals == 0:
            status = "no valid interval"
        elif n_intervals == 1:
            status = "one connected interval"
        else:
            status = "multiple disconnected intervals"

        summary.append({
            'm_phi': m_phi,
            'n_intervals': n_intervals,
            'status': status,
            'intervals': repr(intervals)
        })

    csv_out = output_dir / "m2_interval_summary.csv"
    md_out = output_dir / "m2_interval_summary.md"

    with open(csv_out, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['m_phi', 'n_intervals', 'status', 'intervals'])
        writer.writeheader()
        writer.writerows(summary)

    with open(md_out, 'w') as f:
        f.write("# M2 Interval Summary\n\n")
        f.write("| m_phi | Valid Intervals | Status | Interval Bounds |\n")
        f.write("|---|---|---|---|\n")
        for row in summary:
            f.write(f"| {row['m_phi']} | {row['n_intervals']} | {row['status']} | {row['intervals']} |\n")

    print(f"Topology summary written to {output_dir}")

if __name__ == '__main__':
    main()
