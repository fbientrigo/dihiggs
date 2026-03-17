import os
import glob
import time
import gc
import warnings
from collections import Counter
from typing import Any

import numpy as np
import polars as pl
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

# =========================================================
# CONFIG
# =========================================================
DATA_LAKE_DIR = "/mnt/c/Users/Asus/cern_db/dihiggs_lake"
IMG_DIR = "/home/fabi/dihiggs/mlpython/paper_img/reunion_marzo_polars"

TARGET_COLS = [
    "m_phi",
    "lam1",
    "tan_beta",
    "lambda6",
    "positivity_ok",
    "unitarity_ok",
    "perturbativity_ok",
    "br_gaga",
    "total_width",
]

SCHEMA_OVERRIDES = {
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

# batches por archivo
BATCH_SIZE = 50_000

# conservar solo una muestra pequeña para plots
MAX_SAMPLE_PER_BATCH = 1500
MAX_PLOT_POINTS = 180_000
MAX_SCATTER_POINTS = 120_000

PROGRESS_EVERY = 200

RNG = np.random.default_rng(42)


def format_time(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.2f} seg"
    elif seconds < 3600:
        return f"{seconds/60:.2f} min"
    else:
        return f"{seconds/3600:.2f} hrs"


def safe_read_batches(path: str):
    """
    Devuelve un BatchedCsvReader de Polars o None si falla.
    """
    try:
        reader = pl.read_csv_batched(
            source=path,
            columns=TARGET_COLS,
            schema_overrides=SCHEMA_OVERRIDES,
            infer_schema_length=100,
            batch_size=BATCH_SIZE,
            low_memory=True,
            rechunk=False,
            ignore_errors=True,
            try_parse_dates=False,
        )
        return reader
    except Exception:
        return None


def sample_df(df: pl.DataFrame, nmax: int) -> pl.DataFrame:
    """
    Devuelve una muestra aleatoria pequeña para plots.
    """
    if df.height <= nmax:
        return df
    return df.sample(n=nmax, shuffle=True, seed=int(RNG.integers(0, 2**31 - 1)))


def update_counter_from_grouped(counter: Counter, grouped_df: pl.DataFrame, key_col: str, count_col: str = "len") -> None:
    """
    Actualiza un Counter Python desde un group_by.len() de Polars.
    """
    if grouped_df.height == 0:
        return
    for key, cnt in grouped_df.select([key_col, count_col]).iter_rows():
        counter[key] += int(cnt)


def list_csv_files(data_lake_dir: str) -> list[str]:
    return sorted(glob.glob(os.path.join(data_lake_dir, "**/*.csv"), recursive=True))


def scan_lake(
    data_lake_dir: str = DATA_LAKE_DIR,
    columns: list[str] | None = None,
    schema_overrides: dict | None = None,
) -> pl.LazyFrame:
    """
    Lazy-scan the data lake and return a single LazyFrame.

    Delegates to ``consolidate_lake.ensure_parquet()`` so the data is read
    from a single Parquet on native ext4 — no mmap over 9P, no SIGBUS.
    The parquet is rebuilt automatically only when campaigns/files change.
    """
    from pathlib import Path
    try:
        from consolidate_lake import ensure_parquet
        parquet_path = ensure_parquet(lake_dir=Path(data_lake_dir))
        lf = pl.scan_parquet(parquet_path)
    except ImportError:
        # Fallback: consolidate_lake not available → legacy CSV scan
        print("[WARN] consolidate_lake not found, falling back to CSV scan")
        lf = _scan_lake_csv_legacy(data_lake_dir, schema_overrides)

    if columns is not None:
        available = set(lf.collect_schema().names())
        project = [c for c in columns if c in available]
        if project:
            lf = lf.select(project)

    return lf


def _scan_lake_csv_legacy(
    data_lake_dir: str = DATA_LAKE_DIR,
    schema_overrides: dict | None = None,
) -> pl.LazyFrame:
    """Legacy CSV scanner — kept as fallback.  Prone to SIGBUS on WSL2."""
    csv_files = list_csv_files(data_lake_dir)
    if not csv_files:
        raise FileNotFoundError(f"No CSV files found in: {data_lake_dir}")

    if schema_overrides is None:
        schema_overrides = SCHEMA_OVERRIDES

    lazy_frames: list[pl.LazyFrame] = []
    for path in csv_files:
        try:
            lf = pl.scan_csv(
                path,
                schema_overrides=schema_overrides,
                infer_schema_length=200,
                try_parse_dates=False,
                ignore_errors=True,
                low_memory=True,
            )
            lazy_frames.append(lf)
        except Exception:
            continue

    if not lazy_frames:
        raise RuntimeError(
            f"Could not scan any CSV files from: {data_lake_dir}"
        )

    return pl.concat(lazy_frames, how="vertical_relaxed")


def extract_streaming_stats(
    all_csv_files: list[str],
    progress_every: int = PROGRESS_EVERY,
    batches_per_read: int = 8,
    max_sample_per_batch: int = MAX_SAMPLE_PER_BATCH,
) -> dict[str, Any]:
    print("[*] Iniciando extracción con Polars (batched, memoria controlada)...")
    t0 = time.perf_counter()

    total_rows = 0
    valid_rows = 0
    fail_pos = 0
    fail_uni = 0
    fail_pert = 0
    skipped_files = 0

    lambda6_counts = Counter()
    tanb_counts = Counter()

    plot_samples = []
    total_files = len(all_csv_files)

    for i, fpath in enumerate(all_csv_files, start=1):
        reader = safe_read_batches(fpath)

        if reader is None:
            skipped_files += 1
            continue

        try:
            batches = reader.next_batches(batches_per_read)
        except Exception:
            skipped_files += 1
            continue

        while batches:
            for batch in batches:
                if batch is None or batch.height == 0:
                    continue

                total_rows += batch.height

                stats = batch.select(
                    [
                        (pl.col("positivity_ok") == 0).sum().alias("fail_pos"),
                        (pl.col("unitarity_ok") == 0).sum().alias("fail_uni"),
                        (pl.col("perturbativity_ok") == 0).sum().alias("fail_pert"),
                    ]
                )

                fail_pos += int(stats[0, "fail_pos"])
                fail_uni += int(stats[0, "fail_uni"])
                fail_pert += int(stats[0, "fail_pert"])

                valid = batch.filter(
                    (pl.col("positivity_ok") == 1)
                    & (pl.col("unitarity_ok") == 1)
                    & (pl.col("perturbativity_ok") == 1)
                ).select(["m_phi", "lam1", "tan_beta", "lambda6", "br_gaga", "total_width"])

                n_valid = valid.height
                valid_rows += n_valid

                if n_valid > 0:
                    lambda6_grouped = (
                        valid.with_columns(pl.col("lambda6").round(6)).group_by("lambda6").len()
                    )
                    tanb_grouped = (
                        valid.with_columns(pl.col("tan_beta").round(6)).group_by("tan_beta").len()
                    )

                    update_counter_from_grouped(lambda6_counts, lambda6_grouped, "lambda6")
                    update_counter_from_grouped(tanb_counts, tanb_grouped, "tan_beta")

                    plot_samples.append(sample_df(valid, max_sample_per_batch))

                    del lambda6_grouped, tanb_grouped

                del stats, valid, batch
                gc.collect()

            try:
                batches = reader.next_batches(batches_per_read)
            except Exception:
                break

        del reader
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

    if plot_samples:
        plot_df = pl.concat(plot_samples, rechunk=False)
    else:
        plot_df = pl.DataFrame(
            schema={
                "m_phi": pl.Float32,
                "lam1": pl.Float32,
                "tan_beta": pl.Float32,
                "lambda6": pl.Float32,
                "br_gaga": pl.Float32,
                "total_width": pl.Float32,
            }
        )

    del plot_samples
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
        "plot_df": plot_df,
        "elapsed_total": elapsed_total,
    }


def build_summary_text(stats: dict[str, Any], total_files: int) -> str:
    eff = (100.0 * stats["valid_rows"] / stats["total_rows"]) if stats["total_rows"] > 0 else 0.0

    summary_lines = [
        "=" * 60,
        "ESTADÍSTICAS GLOBALES DEL DATA LAKE",
        "=" * 60,
        f"Archivos CSV encontrados      : {total_files:,}",
        f"Archivos saltados             : {stats['skipped_files']:,}",
        f"Puntos totales computados     : {stats['total_rows']:,}",
        f"Puntos físicamente válidos    : {stats['valid_rows']:,}",
        f"Eficiencia de muestreo        : {eff:.4f}%",
        "",
        "Desglose de fallos:",
        f" - Fallaron positividad       : {stats['fail_pos']:,}",
        f" - Fallaron unitariedad       : {stats['fail_uni']:,}",
        f" - Fallaron perturbatividad   : {stats['fail_pert']:,}",
        "",
        "Top 5 lambda6 con más puntos válidos:",
    ]

    for k, v in stats["lambda6_counts"].most_common(5):
        summary_lines.append(f" - lambda6 = {k}: {v:,}")

    summary_lines.append("")
    summary_lines.append("Top 5 tan_beta con más puntos válidos:")
    for k, v in stats["tanb_counts"].most_common(5):
        summary_lines.append(f" - tan_beta = {k}: {v:,}")

    summary_lines.extend(
        [
            "",
            f"Muestra retenida para plots   : {stats['plot_df'].height:,}",
            f"Tiempo total                  : {format_time(stats['elapsed_total'])}",
            "=" * 60,
        ]
    )

    return "\n".join(summary_lines)


def plot_valid_samples(
    plot_df: pl.DataFrame,
    lambda6_counts: Counter,
    img_dir: str,
    max_scatter_points: int = MAX_SCATTER_POINTS,
) -> None:
    if plot_df.height == 0:
        print("[-] No hubo puntos válidos para graficar.")
        return

    pdf = plot_df.to_pandas()
    gc.collect()

    print("[*] Generando gráficos...")

    plt.figure(figsize=(10, 6))
    hb = plt.hexbin(
        pdf["lam1"],
        pdf["m_phi"],
        gridsize=40,
        bins="log",
        mincnt=1,
        cmap="YlGnBu",
    )
    plt.colorbar(hb, label="log10(N puntos en muestra)")
    plt.title(r"Distribución de puntos válidos (muestra): $\lambda_1$ vs $m_\phi$")
    plt.xlabel(r"$\lambda_1$")
    plt.ylabel(r"$m_\phi$ [GeV]")
    plt.tight_layout()
    out1 = os.path.join(img_dir, "1_espacio_fase_hexbin.png")
    plt.savefig(out1, dpi=250)
    plt.close()
    print(f" [+] Guardado: {out1}")

    top_lambda6 = [k for k, _ in lambda6_counts.most_common(10)]
    vdf = pdf[np.isin(np.round(pdf["lambda6"], 6), top_lambda6)].copy()

    groups = []
    labels = []

    for lv in sorted(np.unique(np.round(vdf["lambda6"].dropna(), 6))):
        arr = vdf.loc[np.round(vdf["lambda6"], 6) == lv, "lam1"].dropna().to_numpy()
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
        print(" [!] Plot 2 omitido: no hubo suficientes grupos por lambda6.")

    df_br = pdf[pdf["br_gaga"] > 1e-15].copy()

    if len(df_br) > max_scatter_points:
        idx = RNG.choice(len(df_br), size=max_scatter_points, replace=False)
        df_br = df_br.iloc[idx].copy()

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
    batches_per_read: int = 8,
    max_sample_per_batch: int = MAX_SAMPLE_PER_BATCH,
    max_plot_points: int = MAX_PLOT_POINTS,
    max_scatter_points: int = MAX_SCATTER_POINTS,
    save_summary: bool = True,
    make_plots: bool = True,
) -> dict[str, Any]:
    print(f"[*] Escaneando el Data Lake en: {data_lake_dir}")
    all_csv_files = list_csv_files(data_lake_dir)
    total_files = len(all_csv_files)
    print(f"[+] Se encontraron {total_files:,} archivos CSV.\n")

    if total_files == 0:
        raise FileNotFoundError(f"No se encontraron CSVs en: {data_lake_dir}")

    os.makedirs(img_dir, exist_ok=True)

    stats = extract_streaming_stats(
        all_csv_files,
        progress_every=progress_every,
        batches_per_read=batches_per_read,
        max_sample_per_batch=max_sample_per_batch,
    )

    plot_df = stats["plot_df"]
    if plot_df.height > max_plot_points:
        plot_df = plot_df.sample(
            n=max_plot_points,
            shuffle=True,
            seed=int(RNG.integers(0, 2**31 - 1)),
        )
        stats["plot_df"] = plot_df

    summary_text = build_summary_text(stats, total_files=total_files)
    print("\n" + summary_text + "\n")

    if save_summary:
        with open(os.path.join(img_dir, "summary.txt"), "w", encoding="utf-8") as f:
            f.write(summary_text)

    if make_plots:
        plot_valid_samples(
            plot_df=stats["plot_df"],
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