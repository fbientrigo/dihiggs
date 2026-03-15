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


def _list_csvs(lake_dir: Path) -> list[str]:
    """Return sorted list of every CSV inside *lake_dir* (recursive)."""
    return sorted(glob.glob(str(lake_dir / "**" / "*.csv"), recursive=True))


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
) -> tuple[Path, int]:
    """
    Read all CSVs in *lake_dir* in small batches and write one Parquet file.

    Strategy (RAM-safe for 7.5 GB):
      Phase 1 — Write per-batch chunk parquets (peak RAM ≈ 1 batch)
      Phase 2 — Lazy-scan all chunks on ext4 → sink into final parquet

    Returns (parquet_path, total_rows).
    """
    csv_files = _list_csvs(lake_dir)
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

        chunk = pl.concat(frames, how="vertical_relaxed", rechunk=True)
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

    # ── Phase 2: merge chunks into single Parquet (lazy, ext4-native) ─
    print(f"\n[consolidate] Phase 2 — merging {len(chunk_paths)} chunks …")

    if not chunk_paths:
        raise RuntimeError("No data was successfully read from the CSVs")

    merged = pl.scan_parquet(
        [str(p) for p in chunk_paths],
    ).collect(streaming=True)

    merged.write_parquet(output, compression="zstd")
    final_rows = merged.height
    del merged
    gc.collect()

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
def needs_rebuild(lake_dir: Path = DATA_LAKE_DIR) -> bool:
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

    csv_files = _list_csvs(lake_dir)
    campaigns = _list_campaigns(lake_dir)
    fp = _fingerprint(csv_files, campaigns)

    if fp != manifest.get("fingerprint"):
        print(f"[consolidate] Cambio detectado (fingerprint mismatch) → rebuild")
        return True

    return False


def ensure_parquet(
    lake_dir: Path = DATA_LAKE_DIR,
    force: bool = False,
) -> Path:
    """
    Ensure the consolidated Parquet exists and is up-to-date.
    Rebuilds only when necessary (or if *force=True*).

    Returns the Path to the Parquet file.
    """
    if force or needs_rebuild(lake_dir):
        parquet_path, n_rows = build_parquet(lake_dir)
        csv_files = _list_csvs(lake_dir)
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
) -> Path:
    """
    Convenience wrapper: ensure Parquet exists and return its Path.

    Use this from notebooks:
        >>> from consolidate_lake import get_parquet_path
        >>> pq = get_parquet_path()
        >>> df = pl.scan_parquet(pq)
    """
    return ensure_parquet(lake_dir, force=force)


# ──────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────
def main() -> None:
    force = "--force" in sys.argv
    parquet = ensure_parquet(force=force)
    print(f"\n{'='*60}")
    print(f"Parquet listo para usar:")
    print(f"  {parquet}")
    print(f"{'='*60}")
    print(f"\nEn notebook:")
    print(f"  from consolidate_lake import get_parquet_path")
    print(f"  df = pl.scan_parquet(get_parquet_path())")


if __name__ == "__main__":
    main()
