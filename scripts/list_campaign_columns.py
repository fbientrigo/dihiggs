from pathlib import Path
import csv

SOURCE_ROOT = Path("/mnt/c/Users/Asus/cern_db/dihiggs_lake")
LIMIT = 10

def list_campaign_dirs(root: Path) -> list[Path]:
    campaigns = []
    for entry in root.iterdir():
        if not entry.is_dir():
            continue
        if entry.name == "recomputed":
            continue
        if entry.name.startswith("campaign="):
            campaigns.append(entry)
    campaigns.sort(key=lambda p: p.name)
    return campaigns

def find_scan_csvs(campaign_dir: Path) -> list[Path]:
    return sorted(campaign_dir.glob("**/tb_*/scan_tb_*.csv"))

def read_header(csv_path: Path) -> list[str]:
    with csv_path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.reader(f)
        header = next(reader, [])
    return [h.strip() for h in header]

campaigns = list_campaign_dirs(SOURCE_ROOT)[:LIMIT]

print(f"Analizando primeras {len(campaigns)} campañas en: {SOURCE_ROOT}\n")

for idx, campaign_dir in enumerate(campaigns, start=1):
    csvs = find_scan_csvs(campaign_dir)
    print("=" * 100)
    print(f"[{idx}] {campaign_dir.name}")
    print(f"scan_csv_count = {len(csvs)}")

    if not csvs:
        print("No se encontraron scan_tb_*.csv")
        print()
        continue

    sample_csv = csvs[0]
    header = read_header(sample_csv)

    print(f"sample_csv = {sample_csv}")
    print(f"n_columns = {len(header)}")
    print("columns =")
    for col in header:
        print(f"  - {col}")

    needed = ["positivity_ok", "unitarity_ok", "perturbativity_ok", "triple_ok"]
    print("key_presence =")
    for key in needed:
        print(f"  - {key}: {'YES' if key in header else 'NO'}")
    print()