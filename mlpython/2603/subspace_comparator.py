import os
import sys
import numpy as np
import polars as pl
import matplotlib.pyplot as plt

import polar_lake_explorer as ple
import consolidate_lake as cl

# =========================
# CONFIGURATION
# =========================
DATA_LAKE_PATH = "/mnt/c/Users/Asus/cern_db/dihiggs_lake"
TEMP_PARQUET_FILE = "temp_subspace.parquet"
IMG_DIR = "subspace_comparisons_logs"
LOG_FILE_PATH = os.path.join(IMG_DIR, "execution_log.txt")

def append_log(msg):
    print(msg)
    with open(LOG_FILE_PATH, "a") as f:
        f.write(msg + "\n")

def _resolve_col(schema_names, aliases):
    lower_map = {c.lower(): c for c in schema_names}
    for a in aliases:
        if a.lower() in lower_map:
            return lower_map[a.lower()]
    for c in schema_names:
        cl_ = c.lower()
        if any(tok in cl_ for tok in aliases):
            return c
    return None

def _numeric_columns(schema):
    numeric_dtypes = {
        pl.Int8, pl.Int16, pl.Int32, pl.Int64,
        pl.UInt8, pl.UInt16, pl.UInt32, pl.UInt64,
        pl.Float32, pl.Float64, pl.Decimal
    }
    cols = []
    for c, dt in schema.items():
        if dt in numeric_dtypes:
            cols.append(c)
    return cols

def _hist_counts_lazy(lf: pl.LazyFrame, col: str, bins: int, label: str = "", is_logx: bool = False):
    schema = lf.collect_schema()
    is_float = schema[col] in (pl.Float32, pl.Float64)
    valid_expr = pl.col(col).is_not_null()
    if is_float:
        valid_expr = valid_expr & ~pl.col(col).is_nan()
        
    if is_logx:
        valid_expr = valid_expr & (pl.col(col) > 0)

    mm = lf.select(
        pl.col(col).filter(valid_expr).min().alias("min_v"),
        pl.col(col).filter(valid_expr).max().alias("max_v"),
        pl.len().alias("total"),
        valid_expr.sum().alias("valid")
    ).collect()

    vmin = mm["min_v"].item()
    vmax = mm["max_v"].item()
    total_c = mm["total"].item()
    valid_c = mm["valid"].item()

    if total_c > valid_c:
        append_log(f"  [!] {label}: {total_c - valid_c} invalid/non-positive/NaN/Null values skipped in column '{col}'.")

    if vmin is None or vmax is None or (is_logx and vmin <= 0):
        append_log(f"  [-] {label}: No valid data found for column '{col}'" + (" (requires strictly positive for logx)" if is_logx else "") + ".")
        return np.array([]), np.array([])

    lf_valid = lf.filter(valid_expr)

    if vmax == vmin:
        edges = np.array([vmin, vmax + (vmax*1e-6 if vmax else 1e-12)])
        counts = np.array([valid_c], dtype=float)
        return edges, counts

    if is_logx:
        log_vmin = np.log10(vmin)
        log_vmax = np.log10(vmax)
        bin_expr = (
            (((pl.col(col).log10() - log_vmin) / (log_vmax - log_vmin)) * bins)
            .floor()
            .cast(pl.Int32)
            .clip(0, bins - 1)
            .alias("_bin")
        )
        edges = np.logspace(log_vmin, log_vmax, bins + 1)
    else:
        bin_expr = (
            (((pl.col(col) - vmin) / (vmax - vmin)) * bins)
            .floor()
            .cast(pl.Int32)
            .clip(0, bins - 1)
            .alias("_bin")
        )
        edges = np.linspace(vmin, vmax, bins + 1)

    h = (
        lf_valid.select(bin_expr)
        .group_by("_bin")
        .len()
        .sort("_bin")
        .collect()
    )

    counts = np.zeros(bins, dtype=float)
    if h.height > 0:
        idx = h["_bin"].to_numpy()
        cnt = h["len"].to_numpy()
        counts[idx] = cnt

    return edges, counts

def plot_single_comparison(df_total_lazy: pl.LazyFrame, df_subset_lazy: pl.LazyFrame, col: str, bins=60, title="", logy=False, logx=False):
    edges_t, cnt_t = _hist_counts_lazy(df_total_lazy, col, bins=bins, label="All Data", is_logx=logx)
    edges_s, cnt_s = _hist_counts_lazy(df_subset_lazy, col, bins=bins, label="Subset", is_logx=logx)

    fig, ax = plt.subplots(figsize=(8, 6))

    if len(cnt_t) == 0 or len(cnt_s) == 0:
        msg = f"  [-] Skipping plot for {col} due to lack of valid data (Total data valid: {'Yes' if len(cnt_t)>0 else 'No'}, Subset valid: {'Yes' if len(cnt_s)>0 else 'No'})."
        append_log(msg)
        ax.set_title(f"{col} (sin datos válidos)")
        ax.axis("off")
        return fig

    edges = edges_t
    if len(edges_t) != len(edges_s) or not np.allclose(edges_t, edges_s):
        vmin, vmax = edges[0], edges[-1]
        b = len(edges) - 1
        
        schema = df_subset_lazy.collect_schema()
        valid_expr = pl.col(col).is_not_null()
        if schema[col] in (pl.Float32, pl.Float64):
            valid_expr = valid_expr & ~pl.col(col).is_nan()

        if logx:
            valid_expr = valid_expr & (pl.col(col) > 0)
            log_vmin = np.log10(vmin)
            log_vmax = np.log10(vmax)
            bin_expr = (
                (((pl.col(col).log10() - log_vmin) / (log_vmax - log_vmin)) * b)
                .floor()
                .cast(pl.Int32)
                .clip(0, b - 1)
                .alias("_bin")
            )
        else:
            bin_expr = (
                (((pl.col(col) - vmin) / (vmax - vmin)) * b)
                .floor()
                .cast(pl.Int32)
                .clip(0, b - 1)
                .alias("_bin")
            )

        hs = (
            df_subset_lazy.filter(valid_expr)
            .select(bin_expr)
            .group_by("_bin")
            .len()
            .sort("_bin")
            .collect()
        )
        cnt_s = np.zeros(b, dtype=float)
        if hs.height > 0:
            cnt_s[hs["_bin"].to_numpy()] = hs["len"].to_numpy()

    if cnt_t.sum() > 0:
        cnt_t = cnt_t / cnt_t.sum()
    if cnt_s.sum() > 0:
        cnt_s = cnt_s / cnt_s.sum()

    centers = 0.5 * (edges[:-1] + edges[1:])

    if logx:
        width = (edges[1:] - edges[:-1])
        ax.bar(edges[:-1], cnt_t, width=width, alpha=0.4, label="All Data", align="edge", color="#1f77b4")
        if hasattr(ax, 'stairs'):
            ax.stairs(cnt_s, edges, color="#ff7f0e", label="Subset", linewidth=2.5)
        else:
            ax.hist(edges[:-1], bins=edges, weights=cnt_s, histtype='step', color="#ff7f0e", label="Subset", linewidth=2.5)
        ax.set_xscale("log")
    else:
        width = (edges[1] - edges[0]) * 0.45
        ax.bar(centers - width / 2, cnt_t, width=width, alpha=0.55, label="All Data", color="#1f77b4")
        ax.bar(centers + width / 2, cnt_s, width=width, alpha=0.55, label="Subset", color="#ff7f0e")

    if logy:
        ax.set_yscale("log")

    ax.set_title(title)
    ax.set_xlabel(col)
    ax.set_ylabel("Density" + (" (Log Scale)" if logy else ""))
    ax.grid(alpha=0.2)
    ax.legend()
    fig.tight_layout()
    return fig

import argparse
from pathlib import Path

def main():
    os.makedirs(IMG_DIR, exist_ok=True)
    if os.path.exists(LOG_FILE_PATH):
        os.remove(LOG_FILE_PATH)
        
    parser = argparse.ArgumentParser(description="Subspace Comparator")
    parser.add_argument("--br", type=float, default=None, help="Threshold for Branching Ratio (BR > th_BR)")
    parser.add_argument("--dw", type=float, default=None, help="Threshold for Total Decay Width (DW < th_DW)")
    parser.add_argument("--force", action="store_true", help="Omit the checking of the parquet and assume it exists/updated")
    parser.add_argument("--vars", type=str, default=None, help="Variables to compare separated by commas, or 'all'")
    parser.add_argument("--logy", action="store_true", help="Make y-axis logarithmic for all plots")
    parser.add_argument("--logx-width", action="store_true", help="Make x-axis logarithmic for 'width_*' variables")
    args = parser.parse_args()

    append_log("==============================================")
    append_log("   Subspace Comparator Interactive Script   ")
    append_log("==============================================")

    if args.force:
        parquet_path = cl.PARQUET_FILE
        append_log(f"[*] Force flag active: Skipping parquet checking, using {parquet_path}")
        if not parquet_path.exists():
            append_log(f"[!] Parquet does not exist at {parquet_path}. Run without --force to build it.")
            sys.exit(1)
    else:
        append_log(f"[*] Ensuring consolidated Parquet for the Data Lake exists...")
        parquet_path = cl.get_parquet_path(lake_dir=Path(DATA_LAKE_PATH))
        
    append_log(f"[*] Loading data lake from Parquet: {parquet_path}")
    df_lake = pl.scan_parquet(parquet_path)
    schema = df_lake.collect_schema()
    schema_names = schema.names()
    append_log(f"[+] Reached lake schema. Found {len(schema_names)} columns.")

    base_filter_expr = None
    expected_filters = ["perturbativity_ok", "unitarity_ok", "positivity_ok"]
    found_filters = []
    
    for req_col in expected_filters:
        col_matched = _resolve_col(schema_names, [req_col])
        if col_matched:
            found_filters.append(col_matched)
            # Support both integer 1 and boolean True
            expr = (pl.col(col_matched) == 1) | (pl.col(col_matched) == True)
            if base_filter_expr is None:
                base_filter_expr = expr
            else:
                base_filter_expr = base_filter_expr & expr

    if base_filter_expr is not None:
        append_log(f"\n[*] Applying mandatory theoretical filters: {', '.join(found_filters)}...")
        count_before = df_lake.select(pl.len()).collect().item()
        df_lake = df_lake.filter(base_filter_expr)
        count_after = df_lake.select(pl.len()).collect().item()
        append_log(f"  -> Rows before filtering: {count_before}")
        append_log(f"  -> Rows after filtering : {count_after}")
        append_log(f"  -> Filtered out         : {count_before - count_after} rows")
    else:
        append_log("\n[!] Mandatory theoretical filters (perturbativity, unitarity, positivity) not found in schema. Proceeding without them.")

    br_col = _resolve_col(schema_names, ["br", "branching_ratio", "br_gaga"])
    tdw_col = _resolve_col(schema_names, ["total_decay_width", "decay_width", "tdw", "total_width"])

    if br_col is None or tdw_col is None:
        append_log("[!] Could not automatically resolve Branching Ratio or Total Decay Width columns.")
        sys.exit(1)
    
    append_log(f"\n[*] Detected Branching Ratio column: {br_col}")
    append_log(f"[*] Detected Total Decay Width column: {tdw_col}")

    if args.br is not None:
        th_br = args.br
        append_log(f"[*] Using provided threshold for Branching Ratio: {th_br}")
    else:
        print("\n[?] Enter threshold for Branching Ratio (BR > th_BR). Default [1e-1]:")
        raw_br = input().strip()
        th_br = float(raw_br) if raw_br else 1e-1

    if args.dw is not None:
        th_dw = args.dw
        append_log(f"[*] Using provided threshold for Total Decay Width: {th_dw}")
    else:
        print("[?] Enter threshold for Total Decay Width (DW < th_DW). Default [1e-11]:")
        raw_dw = input().strip()
        th_dw = float(raw_dw) if raw_dw else 1e-11

    append_log(f"\n[*] Filtering subspace where {br_col} > {th_br} AND {tdw_col} < {th_dw}...")
    
    df_subset = df_lake.filter(
        (pl.col(br_col) > th_br) & (pl.col(tdw_col) < th_dw)
    )

    append_log(f"[*] Extracting full statistics and writing temporal parquet to {TEMP_PARQUET_FILE}...")
    df_subset.sink_parquet(TEMP_PARQUET_FILE)
    append_log(f"[+] Temporal parquet created: {TEMP_PARQUET_FILE}")

    df_subset_lazy = pl.scan_parquet(TEMP_PARQUET_FILE)
    
    num_cols = _numeric_columns(schema)
    
    if args.vars is not None:
        raw_vars = args.vars
        append_log(f"\n[*] Using provided variables to compare: {raw_vars}")
    else:
        print(f"\n[+] Available numeric columns to compare: {', '.join(num_cols)}")
        print("\n[?] Enter variables to compare separated by commas (e.g. m_phi, lam1, tan_beta). Or 'all' for all numeric variables:")
        raw_vars = input().strip()
        if not raw_vars:
            append_log("[-] No variables provided, exiting without plotting.")
            sys.exit(0)

    if raw_vars.lower() == 'all':
        vars_to_compare = num_cols
    else:
        vars_to_compare = [v.strip() for v in raw_vars.split(",") if v.strip() in num_cols]
    
    if not vars_to_compare:
        append_log("[-] Valid variables not found, exiting.")
        sys.exit(0)
    
    append_log(f"\n[*] Generating plots into '{IMG_DIR}/' ...")

    all_lf = df_lake
    sub_lf = df_subset_lazy

    for v in vars_to_compare:
        append_log(f"  -> Plotting comparison for {v}")
        title = f"Comparison of {v}: All vs Subset (BR > {th_br}, DW < {th_dw})"
        
        apply_logx = False
        if args.logx_width and ("width" in v.lower()):
            apply_logx = True

        fig = plot_single_comparison(all_lf, sub_lf, col=v, bins=60, title=title, logy=args.logy, logx=apply_logx)
        
        out_path = os.path.join(IMG_DIR, f"compare_{v}.png")
        fig.savefig(out_path, dpi=200)
        plt.close(fig)

    append_log(f"\n[🚀] Analysis complete. All comparisons and logs saved in directory: {IMG_DIR}.")
    
if __name__ == "__main__":
    main()
