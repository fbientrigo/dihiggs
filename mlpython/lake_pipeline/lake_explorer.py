import os
import glob
import time
import gc
import json
import warnings
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

# =========================================================
# CONFIG
# =========================================================
DEFAULT_DATA_LAKE_DIR = "/mnt/c/Users/Asus/cern_db/dihiggs_lake"
PROJECT_CONFIG_FILE = "project_config.json"


def _load_data_lake_dir(default: str = DEFAULT_DATA_LAKE_DIR) -> str:
    this_file = Path(__file__).resolve()
    cfg_path = next(
        (parent / PROJECT_CONFIG_FILE for parent in [this_file.parent, *this_file.parents] if (parent / PROJECT_CONFIG_FILE).exists()),
        None,
    )
    if cfg_path is None:
        return default
    try:
        payload = json.loads(cfg_path.read_text(encoding="utf-8"))
        configured = payload.get("data_lake_dir") if isinstance(payload, dict) else None
        if isinstance(configured, str) and configured.strip():
            return configured
    except Exception as exc:
        print(f"[WARN] Config inválida en {cfg_path}: {exc}. Se usa ruta por defecto.")
    return default


DATA_LAKE_DIR = _load_data_lake_dir()
IMG_DIR = str(Path(__file__).resolve().parent.parent / "paper_img" / "reunion_marzo_streaming")

TARGET_COLS = [
    "m_phi", "lam1", "tan_beta", "lambda6",
    "positivity_ok", "unitarity_ok", "perturbativity_ok",
    "br_gaga", "total_width"
]

# Tipos compactos para reducir RAM.
# float32 basta para plots/estadísticas rápidas de reunión.
DTYPE_MAP = {
    "m_phi": "float32",
    "lam1": "float32",
    "tan_beta": "float32",
    "lambda6": "float32",
    "positivity_ok": "Int8",
    "unitarity_ok": "Int8",
    "perturbativity_ok": "Int8",
    "br_gaga": "float32",
    "total_width": "float32",
}

# Muestra máxima guardada por archivo para plots
MAX_VALID_SAMPLE_PER_FILE = 100

# Tope global de puntos usados para plots
MAX_PLOT_POINTS = 350_000

# Para el scatter de BR vs width
MAX_SCATTER_POINTS = 120_000

# Progreso
PROGRESS_EVERY = 250

RNG = np.random.default_rng(42)


def format_time(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.2f} seg"
    elif seconds < 3600:
        return f"{seconds/60:.2f} min"
    else:
        return f"{seconds/3600:.2f} hrs"


def read_csv_fast(path: str) -> pd.DataFrame | None:
    """
    Lee un CSV cargando solo las columnas útiles.
    Devuelve None si falla.
    """
    try:
        df = pd.read_csv(
            path,
            usecols=lambda c: c in TARGET_COLS,
            dtype=DTYPE_MAP,
            engine="c",
        )
        if df.empty:
            return None

        # Asegurarse de que estén las columnas mínimas necesarias
        needed = {
            "m_phi", "lam1", "tan_beta", "lambda6",
            "positivity_ok", "unitarity_ok", "perturbativity_ok",
            "br_gaga", "total_width"
        }
        if not needed.issubset(df.columns):
            return None

        return df

    except Exception:
        return None


def compact_sample(valid_df: pd.DataFrame, nmax: int) -> pd.DataFrame:
    """
    Guarda solo una muestra pequeña por archivo para no explotar RAM.
    """
    n = len(valid_df)
    if n == 0:
        return valid_df.iloc[:0].copy()
    if n <= nmax:
        return valid_df.copy()
    seed = int(RNG.integers(0, 2**31 - 1))
    return valid_df.sample(n=nmax, random_state=seed).copy()


def save_summary_txt(path: str, text: str) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)


def list_csv_files(data_lake_dir: str) -> list[str]:
    return sorted(glob.glob(os.path.join(data_lake_dir, "**/*.csv"), recursive=True))


def load_lake_df(
    data_lake_dir: str = DATA_LAKE_DIR,
    columns: list[str] | None = None,
) -> pd.DataFrame:
    """
    Load the data lake as a Pandas DataFrame via the consolidated Parquet.

    Uses ``consolidate_lake.ensure_parquet()`` under the hood so the read
    comes from native ext4 — no mmap over 9P.
    """
    from pathlib import Path
    try:
        import polars as _pl
        from consolidate_lake import ensure_parquet
        parquet_path = ensure_parquet(lake_dir=Path(data_lake_dir))
        cols = columns or TARGET_COLS
        lf = _pl.scan_parquet(parquet_path)
        available = set(lf.collect_schema().names())
        project = [c for c in cols if c in available]
        return lf.select(project).collect().to_pandas()
    except ImportError:
        print("[WARN] consolidate_lake not found, falling back to CSV scan")
        all_csv = list_csv_files(data_lake_dir)
        parts = [read_csv_fast(f) for f in all_csv]
        return pd.concat([p for p in parts if p is not None], ignore_index=True)


def extract_streaming_stats(
    all_csv_files: list[str],
    progress_every: int = PROGRESS_EVERY,
    max_valid_sample_per_file: int = MAX_VALID_SAMPLE_PER_FILE,
) -> dict[str, Any]:
    print("[*] Iniciando extracción STREAMING (sin cargar todo a RAM)...")
    t0 = time.perf_counter()

    total_rows = 0
    valid_rows = 0
    fail_pos = 0
    fail_uni = 0
    fail_pert = 0
    skipped_files = 0

    lambda6_counts = Counter()
    tanb_counts = Counter()

    sampled_valid_parts = []
    total_files = len(all_csv_files)

    for i, fpath in enumerate(all_csv_files, start=1):
        df = read_csv_fast(fpath)

        if df is None:
            skipped_files += 1
            continue

        total_rows += len(df)

        fail_pos += int((df["positivity_ok"] == 0).sum())
        fail_uni += int((df["unitarity_ok"] == 0).sum())
        fail_pert += int((df["perturbativity_ok"] == 0).sum())

        valid_mask = (
            (df["positivity_ok"] == 1)
            & (df["unitarity_ok"] == 1)
            & (df["perturbativity_ok"] == 1)
        )

        valid_df = df.loc[
            valid_mask,
            ["m_phi", "lam1", "tan_beta", "lambda6", "br_gaga", "total_width"],
        ].copy()

        n_valid = len(valid_df)
        valid_rows += n_valid

        if n_valid > 0:
            lambda6_counts.update(valid_df["lambda6"].round(6).value_counts().to_dict())
            tanb_counts.update(valid_df["tan_beta"].round(6).value_counts().to_dict())
            sampled_valid_parts.append(compact_sample(valid_df, max_valid_sample_per_file))

        del df, valid_df, valid_mask
        gc.collect()

        if (i % progress_every == 0) or (i == total_files):
            elapsed = time.perf_counter() - t0
            rate = i / elapsed if elapsed > 0 else 0.0
            eta = (total_files - i) / rate if rate > 0 else float("inf")
            print(
                f"[{i:>5,}/{total_files:,}] "
                f"elapsed={format_time(elapsed)} | "
                f"ETA={format_time(eta)} | "
                f"rows={total_rows:,} | valid={valid_rows:,} | skipped={skipped_files:,}"
            )

    elapsed_total = time.perf_counter() - t0
    print(f"\n[+] Extracción completada en {format_time(elapsed_total)}.")

    if sampled_valid_parts:
        sampled_valid = pd.concat(sampled_valid_parts, ignore_index=True)
    else:
        sampled_valid = pd.DataFrame(
            columns=["m_phi", "lam1", "tan_beta", "lambda6", "br_gaga", "total_width"]
        )

    del sampled_valid_parts
    gc.collect()

    return {
        "total_rows": total_rows,
        "valid_rows": valid_rows,
        "fail_pos": fail_pos,
        "fail_uni": fail_uni,
        "fail_pert": fail_pert,
        "skipped_files": skipped_files,
        "lambda6_counts": lambda6_counts,
        "tanb_counts": tanb_counts,
        "sampled_valid": sampled_valid,
        "elapsed_total": elapsed_total,
    }


def _extract_stats_from_df(
    full_df: pd.DataFrame,
    max_sample: int = MAX_VALID_SAMPLE_PER_FILE,
) -> dict[str, Any]:
    """Compute stats dict from a full DataFrame (loaded from Parquet)."""
    import time as _time

    t0 = _time.perf_counter()
    total_rows = len(full_df)

    fail_pos = int((full_df["positivity_ok"] == 0).sum()) if "positivity_ok" in full_df.columns else 0
    fail_uni = int((full_df["unitarity_ok"] == 0).sum()) if "unitarity_ok" in full_df.columns else 0
    fail_pert = int((full_df["perturbativity_ok"] == 0).sum()) if "perturbativity_ok" in full_df.columns else 0

    valid_mask = (
        (full_df.get("positivity_ok", 1) == 1)
        & (full_df.get("unitarity_ok", 1) == 1)
        & (full_df.get("perturbativity_ok", 1) == 1)
    )

    plot_cols = ["m_phi", "lam1", "tan_beta", "lambda6", "br_gaga", "total_width"]
    available_plot_cols = [c for c in plot_cols if c in full_df.columns]
    valid_df = full_df.loc[valid_mask, available_plot_cols].copy()
    valid_rows = len(valid_df)

    lambda6_counts = Counter()
    tanb_counts = Counter()
    if valid_rows > 0 and "lambda6" in valid_df.columns:
        lambda6_counts.update(valid_df["lambda6"].round(6).value_counts().to_dict())
    if valid_rows > 0 and "tan_beta" in valid_df.columns:
        tanb_counts.update(valid_df["tan_beta"].round(6).value_counts().to_dict())

    # Sample for plotting
    sample_size = min(valid_rows, max_sample * 500)
    if valid_rows > sample_size:
        sampled_valid = valid_df.sample(n=sample_size, random_state=42)
    else:
        sampled_valid = valid_df

    elapsed_total = _time.perf_counter() - t0

    return {
        "total_rows": total_rows,
        "valid_rows": valid_rows,
        "fail_pos": fail_pos,
        "fail_uni": fail_uni,
        "fail_pert": fail_pert,
        "skipped_files": 0,
        "lambda6_counts": lambda6_counts,
        "tanb_counts": tanb_counts,
        "sampled_valid": sampled_valid,
        "elapsed_total": elapsed_total,
    }




def build_summary_text(stats: dict[str, Any], total_files: int) -> str:
    eff = (100.0 * stats["valid_rows"] / stats["total_rows"]) if stats["total_rows"] > 0 else 0.0

    summary = []
    summary.append("=" * 60)
    summary.append("ESTADÍSTICAS GLOBALES DEL DATA LAKE")
    summary.append("=" * 60)
    summary.append(f"Archivos CSV encontrados      : {total_files:,}")
    summary.append(f"Archivos saltados             : {stats['skipped_files']:,}")
    summary.append(f"Puntos totales computados     : {stats['total_rows']:,}")
    summary.append(f"Puntos físicamente válidos    : {stats['valid_rows']:,}")
    summary.append(f"Eficiencia de muestreo        : {eff:.4f}%")
    summary.append("")
    summary.append("Desglose de fallos:")
    summary.append(f" - Fallaron positividad       : {stats['fail_pos']:,}")
    summary.append(f" - Fallaron unitariedad       : {stats['fail_uni']:,}")
    summary.append(f" - Fallaron perturbatividad   : {stats['fail_pert']:,}")
    summary.append("")
    summary.append("Top 5 lambda6 con más puntos válidos:")
    for k, v in stats["lambda6_counts"].most_common(5):
        summary.append(f" - lambda6 = {k}: {v:,}")
    summary.append("")
    summary.append("Top 5 tan_beta con más puntos válidos:")
    for k, v in stats["tanb_counts"].most_common(5):
        summary.append(f" - tan_beta = {k}: {v:,}")
    summary.append("")
    summary.append(f"Muestra retenida para plots   : {len(stats['sampled_valid']):,}")
    summary.append(f"Tiempo total                  : {format_time(stats['elapsed_total'])}")
    summary.append("=" * 60)
    return "\n".join(summary)


def plot_valid_samples(
    sampled_valid: pd.DataFrame,
    lambda6_counts: Counter,
    img_dir: str,
    max_scatter_points: int = MAX_SCATTER_POINTS,
) -> None:
    if len(sampled_valid) == 0:
        print("[-] No hay puntos válidos para graficar.")
        return

    print("[*] Generando gráficos...")

    plt.figure(figsize=(10, 6))
    hb = plt.hexbin(
        sampled_valid["lam1"],
        sampled_valid["m_phi"],
        gridsize=40,
        bins="log",
        mincnt=1,
        cmap="YlGnBu",
    )
    plt.colorbar(hb, label="log10(N puntos de la muestra)")
    plt.title(r"Distribución de puntos válidos (muestra): $\lambda_1$ vs $m_\phi$")
    plt.xlabel(r"$\lambda_1$")
    plt.ylabel(r"$m_\phi$ [GeV]")
    plt.tight_layout()
    out1 = os.path.join(img_dir, "1_espacio_fase_hexbin.png")
    plt.savefig(out1, dpi=250)
    plt.close()
    print(f" [+] Guardado: {out1}")

    top_lambda6 = [k for k, _ in lambda6_counts.most_common(10)]

    vdf = sampled_valid.copy()
    vdf["lambda6_round"] = vdf["lambda6"].round(6)
    vdf = vdf[vdf["lambda6_round"].isin(top_lambda6)].copy()

    groups = []
    labels = []

    for lv in sorted(vdf["lambda6_round"].dropna().unique()):
        arr = vdf.loc[vdf["lambda6_round"] == lv, "lam1"].dropna().to_numpy()
        if len(arr) > 4000:
            arr = RNG.choice(arr, size=4000, replace=False)
        if len(arr) >= 5:
            groups.append(arr)
            labels.append(f"{lv:g}")

    if groups:
        plt.figure(figsize=(max(10, 0.9 * len(labels) + 2), 6))
        plt.violinplot(groups, showmedians=True, showextrema=False)
        plt.xticks(np.arange(1, len(labels) + 1), labels, rotation=45)
        plt.title(r"Rango permitido de $\lambda_1$ condicionado por los $\lambda_6$ más poblados")
        plt.xlabel(r"$\lambda_6$")
        plt.ylabel(r"$\lambda_1$ (válidos, muestra)")
        plt.tight_layout()
        out2 = os.path.join(img_dir, "2_lambda6_vs_lam1_violin_top10.png")
        plt.savefig(out2, dpi=250)
        plt.close()
        print(f" [+] Guardado: {out2}")
    else:
        print(" [!] Plot 2 omitido: no hubo suficientes datos agrupables por lambda6.")

    df_br = sampled_valid[sampled_valid["br_gaga"] > 1e-15].copy()

    if len(df_br) > max_scatter_points:
        seed = int(RNG.integers(0, 2**31 - 1))
        df_br = df_br.sample(n=max_scatter_points, random_state=seed).copy()

    if len(df_br) > 0:
        plt.figure(figsize=(10, 6))
        plt.scatter(df_br["total_width"], df_br["br_gaga"], s=8, alpha=0.25)
        plt.xscale("log")
        plt.yscale("log")
        plt.axvline(x=1e-10, linestyle="--", color="red", label=r"Región long-lived ($\Gamma < 10^{-10}$ GeV)")
        plt.title(r"Correlación: anchura total vs BR($h \to \gamma\gamma$)")
        plt.xlabel(r"$\Gamma_{total}$ [GeV]")
        plt.ylabel(r"BR($h \to \gamma\gamma$)")
        plt.legend()
        plt.tight_layout()
        out3 = os.path.join(img_dir, "3_long_lived_particles.png")
        plt.savefig(out3, dpi=250)
        plt.close()
        print(f" [+] Guardado: {out3}")
    else:
        print(" [!] Plot 3 omitido: no hubo puntos con br_gaga > 1e-15.")


def run_pipeline(
    data_lake_dir: str = DATA_LAKE_DIR,
    img_dir: str = IMG_DIR,
    progress_every: int = PROGRESS_EVERY,
    max_valid_sample_per_file: int = MAX_VALID_SAMPLE_PER_FILE,
    max_plot_points: int = MAX_PLOT_POINTS,
    max_scatter_points: int = MAX_SCATTER_POINTS,
    save_summary: bool = True,
    make_plots: bool = True,
    use_parquet: bool = True,
) -> dict[str, Any]:
    os.makedirs(img_dir, exist_ok=True)

    if use_parquet:
        # Fast path: load from consolidated Parquet (ext4, no SIGBUS)
        print(f"[*] Loading data lake from consolidated Parquet …")
        full_df = load_lake_df(data_lake_dir)
        total_files = len(list_csv_files(data_lake_dir))
        print(f"[+] Loaded {len(full_df):,} rows ({total_files:,} source CSVs)\n")
        stats = _extract_stats_from_df(full_df, max_valid_sample_per_file)
    else:
        # Legacy path: iterate CSV by CSV
        print(f"[*] Escaneando el Data Lake en: {data_lake_dir}")
        all_csv_files = list_csv_files(data_lake_dir)
        total_files = len(all_csv_files)
        print(f"[+] Se encontraron {total_files:,} archivos CSV.\n")
        if total_files == 0:
            raise FileNotFoundError(f"No hay archivos CSV en: {data_lake_dir}")
        stats = extract_streaming_stats(
            all_csv_files,
            progress_every=progress_every,
            max_valid_sample_per_file=max_valid_sample_per_file,
        )

    sampled_valid = stats["sampled_valid"]
    if len(sampled_valid) > max_plot_points:
        seed = int(RNG.integers(0, 2**31 - 1))
        sampled_valid = sampled_valid.sample(n=max_plot_points, random_state=seed).copy()
        stats["sampled_valid"] = sampled_valid

    summary_text = build_summary_text(stats, total_files=total_files)
    print("\n" + summary_text + "\n")

    if save_summary:
        save_summary_txt(os.path.join(img_dir, "summary.txt"), summary_text)

    if make_plots:
        plot_valid_samples(
            sampled_valid=sampled_valid,
            lambda6_counts=stats["lambda6_counts"],
            img_dir=img_dir,
            max_scatter_points=max_scatter_points,
        )

    print("\n[🚀] Listo. Estadísticas y gráficos exportados.")
    print(f"[*] Revisa: {img_dir}")

    return {
        "total_files": total_files,
        "summary_text": summary_text,
        **stats,
    }


def main() -> None:
    run_pipeline()


if __name__ == "__main__":
    main()