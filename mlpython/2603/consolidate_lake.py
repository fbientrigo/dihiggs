"""
consolidate_lake.py
====================
One-shot CSV → Parquet consolidator for the di-Higgs data lake.

Fixes the SIGBUS crash caused by mmap'ing 7 000+ CSV files over WSL2's 9P
filesystem.  Reads in small batches (no mmap), writes to native ext4 as a
single Parquet file.

Usage
-----
    # Build/refresh the consolidated Parquet (run from terminal or notebook)
    python consolidate_lake.py

    # From a notebook — just get the Parquet path
    from consolidate_lake import get_parquet_path
    parquet = get_parquet_path()           # returns Path, builds if needed

Smart-cache:  A lightweight manifest (JSON) records the set of CSVs and
campaign dirs.  On the next run the script compares the current snapshot
against the manifest and **only rebuilds** when something changed.
"""

from __future__ import annotations

import argparse
import gc
import glob
import hashlib
import json
import os
import sys
import time
from pathlib import Path
from typing import Any

import polars as pl

# ──────────────────────────────────────────────────────────────────────
# CONFIG
# ──────────────────────────────────────────────────────────────────────
DATA_LAKE_DIR = Path("/mnt/c/Users/Asus/cern_db/dihiggs_lake")

# Output lives on native ext4 → no 9P, no SIGBUS
OUTPUT_DIR = Path.home() / "cern_db" / "dihiggs_consolidated"
PARQUET_FILE = OUTPUT_DIR / "dihiggs_lake.parquet"
MANIFEST_FILE = OUTPUT_DIR / "manifest.json"

# How many CSVs to concat per write-batch (keeps peak RAM ≈ 400-600 MB)
BATCH_SIZE = 60

# Progress reporting
PROGRESS_EVERY = 100

# Schema overrides — float32 is enough for EDA and saves ~50 % RAM
SCHEMA_OVERRIDES: dict[str, pl.DataType] = {
    "m_phi": pl.Float32,
    "lam1": pl.Float32,
    "tan_beta": pl.Float32,
    "lambda6": pl.Float32,
    "positivity_ok": pl.Int8,
    "unitarity_ok": pl.Int8,
    "perturbativity_ok": pl.Int8,
    "br_gaga": pl.Float32,
    "total_width": pl.Float32,
}


# ──────────────────────────────────────────────────────────────────────
# HELPERS
# ──────────────────────────────────────────────────────────────────────
def _fmt(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.1f}s"
    if seconds < 3600:
        return f"{seconds / 60:.1f}min"
    return f"{seconds / 3600:.1f}h"


def _list_csvs(lake_dir: Path, frac: float = 1.0) -> list[str]:
    """Return sorted list of every CSV inside *lake_dir* (recursive), optionally sampling a fraction per campaign."""
    if frac >= 1.0:
        return sorted(glob.glob(str(lake_dir / "**" / "*.csv"), recursive=True))
    
    all_csvs = []
    campaigns = _list_campaigns(lake_dir)
    if not campaigns:
        # Fallback if no subdirectories
        csvs = sorted(glob.glob(str(lake_dir / "*.csv")))
        n = max(1, int(len(csvs) * frac)) if csvs else 0
        return csvs[:n]
        
    for camp in campaigns:
        camp_dir = lake_dir / camp
        camp_csvs = sorted(glob.glob(str(camp_dir / "**" / "*.csv"), recursive=True))
        n_take = max(1, int(len(camp_csvs) * frac)) if camp_csvs else 0
        all_csvs.extend(camp_csvs[:n_take])
    return all_csvs


def _list_campaigns(lake_dir: Path) -> list[str]:
    """Return first-level subdirectory names (campaign folders)."""
    return sorted(
        d.name
        for d in lake_dir.iterdir()
        if d.is_dir()
    )


def _fingerprint(csv_files: list[str], campaigns: list[str]) -> str:
    """Cheap hash combining file count, campaign names, and first+last paths."""
    blob = json.dumps({
        "n_files": len(csv_files),
        "campaigns": campaigns,
        "first_5": csv_files[:5],
        "last_5": csv_files[-5:],
    }, sort_keys=True)
    return hashlib.sha256(blob.encode()).hexdigest()


def _load_manifest() -> dict[str, Any] | None:
    if MANIFEST_FILE.exists():
        try:
            return json.loads(MANIFEST_FILE.read_text())
        except Exception:
            return None
    return None


def _save_manifest(
    fingerprint: str,
    n_files: int,
    campaigns: list[str],
    n_rows: int,
    elapsed: float,
) -> None:
    MANIFEST_FILE.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST_FILE.write_text(
        json.dumps(
            {
                "fingerprint": fingerprint,
                "n_files": n_files,
                "campaigns": campaigns,
                "n_rows": n_rows,
                "parquet": str(PARQUET_FILE),
                "elapsed_s": round(elapsed, 2),
                "created_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            },
            indent=2,
        )
    )


# ──────────────────────────────────────────────────────────────────────
# CORE — batch reader (no mmap, RAM-safe)
# ──────────────────────────────────────────────────────────────────────
def _read_csv_safe(path: str) -> pl.DataFrame | None:
    """Read a single CSV with controlled memory.  Returns None on failure."""
    try:
        return pl.read_csv(
            path,
            schema_overrides=SCHEMA_OVERRIDES,
            infer_schema_length=200,
            ignore_errors=True,
            try_parse_dates=False,
            low_memory=False,     # ← CRITICAL: avoid mmap over 9P!
            rechunk=False,
        )
    except Exception as exc:
        print(f"  ⚠ skip: {os.path.basename(path)} ({exc})")
        return None


def build_parquet(
    lake_dir: Path = DATA_LAKE_DIR,
    output: Path = PARQUET_FILE,
    batch_size: int = BATCH_SIZE,
    frac: float = 1.0,
) -> tuple[Path, int]:
    """
    Read all CSVs in *lake_dir* in small batches and write one Parquet file.

    Strategy (RAM-safe for 7.5 GB):
      Phase 1 — Write per-batch chunk parquets (peak RAM ≈ 1 batch)
      Phase 2 — Lazy-scan all chunks on ext4 → sink into final parquet

    Returns (parquet_path, total_rows).
    """
    csv_files = _list_csvs(lake_dir, frac=frac)
    n_files = len(csv_files)

    if n_files == 0:
        raise FileNotFoundError(f"No CSV files found under {lake_dir}")

    print(f"[consolidate] {n_files:,} CSVs found in {lake_dir}")
    print(f"[consolidate] Output → {output}")
    print(f"[consolidate] Batch size = {batch_size}\n")

    output.parent.mkdir(parents=True, exist_ok=True)

    # ── Phase 1: CSV batches → individual chunk parquets ──────────
    chunks_dir = output.parent / "_chunks"
    chunks_dir.mkdir(exist_ok=True)

    total_rows = 0
    chunk_paths: list[Path] = []
    t0 = time.perf_counter()

    for batch_start in range(0, n_files, batch_size):
        batch_end = min(batch_start + batch_size, n_files)
        batch_paths = csv_files[batch_start:batch_end]

        frames: list[pl.DataFrame] = []
        for fpath in batch_paths:
            df = _read_csv_safe(fpath)
            if df is not None and df.height > 0:
                frames.append(df)

        if not frames:
            continue

        chunk = pl.concat(frames, how="diagonal_relaxed", rechunk=True)
        total_rows += chunk.height

        chunk_path = chunks_dir / f"chunk_{batch_start:06d}.parquet"
        chunk.write_parquet(chunk_path, compression="zstd")
        chunk_paths.append(chunk_path)

        del frames, chunk
        gc.collect()

        # Progress
        done = batch_end
        if done % PROGRESS_EVERY < batch_size or done == n_files:
            elapsed = time.perf_counter() - t0
            rate = done / elapsed if elapsed > 0 else 0
            eta = (n_files - done) / rate if rate > 0 else 0
            print(
                f"  [{done:>6,}/{n_files:,}] "
                f"rows={total_rows:>12,}  "
                f"elapsed={_fmt(elapsed)}  "
                f"ETA≈{_fmt(eta)}"
            )

    # ── Phase 2: merge chunks into single Parquet (DuckDB, ext4-native)
    print(f"\\n[consolidate] Phase 2 — merging {len(chunk_paths)} chunks (streaming with DuckDB) …")

    if not chunk_paths:
        raise RuntimeError("No data was successfully read from the CSVs")

    import duckdb
    
    # We use DuckDB for the final merge to prevent OOM (Out Of Memory) issues 
    # when dealing with dozens of millions of rows.
    con = duckdb.connect(config={'memory_limit': '4GB'})
    
    # Using glob pattern for duckdb so it reads all chunks at once efficiently
    chunks_glob = str(chunks_dir / "*.parquet")
    
    # We write directly to the final .parquet file using ZSTD compression
    query = f"""
    COPY (
        SELECT * FROM read_parquet('{chunks_glob}')
    ) TO '{output}' (FORMAT PARQUET, COMPRESSION 'ZSTD')
    """
    con.execute(query)
    
    # Get final rows count quickly from DuckDB
    final_rows = con.execute(f"SELECT count(*) FROM read_parquet('{output}')").fetchone()[0]
    con.close()

    # Clean up chunk files
    for p in chunk_paths:
        p.unlink(missing_ok=True)
    if chunks_dir.exists() and not any(chunks_dir.iterdir()):
        chunks_dir.rmdir()

    elapsed_total = time.perf_counter() - t0
    print(f"\n[✔] Consolidated Parquet written: {output}")
    print(f"    Rows: {final_rows:,}  |  Time: {_fmt(elapsed_total)}")
    size_mb = output.stat().st_size / (1024 * 1024)
    print(f"    Size: {size_mb:,.1f} MB (zstd compressed)")

    return output, final_rows


# ──────────────────────────────────────────────────────────────────────
# PUBLIC API
# ──────────────────────────────────────────────────────────────────────
def needs_rebuild(lake_dir: Path = DATA_LAKE_DIR, frac: float = 1.0) -> bool:
    """
    Return True if the Parquet must be (re)built.

    Checks:
      1. Parquet file exists?
      2. Fingerprint of CSV list + campaigns matches manifest?
    """
    if not PARQUET_FILE.exists():
        return True

    manifest = _load_manifest()
    if manifest is None:
        return True

    csv_files = _list_csvs(lake_dir, frac=frac)
    campaigns = _list_campaigns(lake_dir)
    fp = _fingerprint(csv_files, campaigns)

    if fp != manifest.get("fingerprint"):
        print(f"[consolidate] Cambio detectado (fingerprint mismatch) → rebuild")
        return True

    return False


def ensure_parquet(
    lake_dir: Path = DATA_LAKE_DIR,
    force: bool = False,
    frac: float = 1.0,
) -> Path:
    """
    Ensure the consolidated Parquet exists and is up-to-date.
    Rebuilds only when necessary (or if *force=True*).

    Returns the Path to the Parquet file.
    """
    if force or needs_rebuild(lake_dir, frac=frac):
        parquet_path, n_rows = build_parquet(lake_dir, frac=frac)
        csv_files = _list_csvs(lake_dir, frac=frac)
        campaigns = _list_campaigns(lake_dir)
        fp = _fingerprint(csv_files, campaigns)
        _save_manifest(fp, len(csv_files), campaigns, n_rows, 0.0)
    else:
        parquet_path = PARQUET_FILE
        manifest = _load_manifest()
        n_rows = manifest.get("n_rows", 0) if manifest else 0
        print(f"[consolidate] Parquet up-to-date ✔  ({n_rows:,} rows)")
        print(f"              {parquet_path}")

    return parquet_path


def get_parquet_path(
    lake_dir: Path = DATA_LAKE_DIR,
    force: bool = False,
    frac: float = 1.0,
) -> Path:
    """
    Convenience wrapper: ensure Parquet exists and return its Path.

    Use this from notebooks:
        >>> from consolidate_lake import get_parquet_path
        >>> pq = get_parquet_path()
        >>> df = pl.scan_parquet(pq)
    """
    return ensure_parquet(lake_dir, force=force, frac=frac)


# ──────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────
def main() -> None:
    parser = argparse.ArgumentParser(description="One-shot CSV -> Parquet consolidator for di-Higgs lake.")
    parser.add_argument("--force", action="store_true", help="Force rebuild of Parquet file")
    parser.add_argument("--frac", type=float, default=1.0, help="Fraction of files to take per campaign (e.g. 0.1 for 10%)")
    args = parser.parse_args()

    parquet = ensure_parquet(force=args.force, frac=args.frac)
    print(f"\n{'='*60}")
    print(f"Parquet listo para usar:")
    print(f"  {parquet}")
    print(f"{'='*60}")
    print(f"\nEn notebook:")
    print(f"  from consolidate_lake import get_parquet_path")
    print(f"  df = pl.scan_parquet(get_parquet_path())")


if __name__ == "__main__":
    main()
