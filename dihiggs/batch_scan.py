import os
import subprocess
import argparse
import csv
from datetime import datetime

def run_scan(exe_path, params, output_dir):
    """Run one scan and log time and output."""
    fname = (
        f"out_mphi{params['mphi_min']}-{params['mphi_max']}_"
        f"Nmphi{params['N_mphi']}_m12{params['m12_min']}-{params['m12_max']}_"
        f"Nm12{params['N_m12']}_mA{params['mA']}_sba{params['sin_ba']}_"
        f"tanb{params['tan_beta']}_l6{params['lambda6']}_l7{params['lambda7']}.csv"
    )
    out_path = os.path.join(output_dir, fname)
    cmd = [
        exe_path,
        str(params['mphi_min']), str(params['mphi_max']), str(params['N_mphi']),
        str(params['m12_min']), str(params['m12_max']), str(params['N_m12']),
        str(params['mA']), str(params['sin_ba']), str(params['tan_beta']),
        str(params['lambda6']), str(params['lambda7']),
        out_path
    ]
    start = datetime.now()
    subprocess.run(cmd, check=True)
    end = datetime.now()
    elapsed_s = (end - start).total_seconds()
    return fname, elapsed_s

def main():
    parser = argparse.ArgumentParser(
        description="Run PhysScanWithFixings_speed_test for preset configs"
    )
    parser.add_argument("--exe", default="./app/PhysScanWithFixings_speed_test",
                        help="Path to the executable")
    parser.add_argument("--config", required=True,
                        help="CSV file with predefined parameter sets")
    parser.add_argument("--output-dir", default="results",
                        help="Directory to store CSV outputs")
    parser.add_argument("--summary", default="summary.csv",
                        help="Filename for summary CSV in output-dir")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    summary_path = os.path.join(args.output_dir, args.summary)

    with open(args.config, newline='') as cfg, \
         open(summary_path, "w", newline='') as summary_file:
        reader = csv.DictReader(cfg)
        writer = csv.writer(summary_file)
        writer.writerow([
            "output_file", "mphi_min", "mphi_max", "N_mphi",
            "m12_min", "m12_max", "N_m12", "mA", "sin_ba",
            "tan_beta", "lambda6", "lambda7", "time_s"
        ])

        for row in reader:
            # parse types
            params = {
                'mphi_min':   float(row['mphi_min']),
                'mphi_max':   float(row['mphi_max']),
                'N_mphi':     int(row['N_mphi']),
                'm12_min':    float(row['m12_min']),
                'm12_max':    float(row['m12_max']),
                'N_m12':      int(row['N_m12']),
                'mA':         float(row['mA']),
                'sin_ba':     float(row['sin_ba']),
                'tan_beta':   float(row['tan_beta']),
                'lambda6':    float(row['lambda6']),
                'lambda7':    float(row['lambda7'])
            }
            fname, t = run_scan(args.exe, params, args.output_dir)
            writer.writerow([
                fname,
                params['mphi_min'], params['mphi_max'], params['N_mphi'],
                params['m12_min'], params['m12_max'], params['N_m12'],
                params['mA'], params['sin_ba'], params['tan_beta'],
                params['lambda6'], params['lambda7'], f"{t:.2f}"
            ])
            summary_file.flush()

if __name__ == "__main__":
    main()

# Ejemplo de ejecución
#./batch_scan.py \
#  --exe ./app/PhysScanWithFixings_speed_test \
#  --config configs.csv \
#  --output-dir results
#
#