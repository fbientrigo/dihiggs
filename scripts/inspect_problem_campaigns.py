from pathlib import Path
import csv

SOURCE_ROOT = Path("/mnt/c/Users/Asus/cern_db/dihiggs_lake")
CAMPAIGNS = [
    "campaign=lc_test",
    "campaign=legacy_recovered_v1",
]

KEYS = ["positivity_ok", "unitarity_ok", "perturbativity_ok", "lam1"]

def find_scan_csvs(campaign_dir: Path) -> list[Path]:
    return sorted(campaign_dir.glob("**/tb_*/scan_tb_*.csv"))

for campaign_name in CAMPAIGNS:
    campaign_dir = SOURCE_ROOT / campaign_name
    print("=" * 120)
    print(campaign_name)

    if not campaign_dir.exists():
        print("  campaign directory not found")
        continue

    csvs = find_scan_csvs(campaign_dir)
    print(f"  scan_csv_count = {len(csvs)}")

    if not csvs:
        print("  no scan_tb_*.csv found")
        continue

    sample_csv = csvs[0]
    print(f"  sample_csv = {sample_csv}")

    with sample_csv.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames or []

        print(f"  n_columns = {len(header)}")
        print("  has_keys =")
        for key in KEYS:
            print(f"    - {key}: {'YES' if key in header else 'NO'}")

        print("  first_5_rows_selected_fields =")
        for i, row in enumerate(reader):
            if i >= 5:
                break
            view = {k: row.get(k, "") for k in KEYS}
            print(f"    row_{i}: {view}")