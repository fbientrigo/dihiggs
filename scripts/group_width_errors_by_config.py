#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd

CHANNEL_COLS = {
    "bb": "rel_err_width_bb",
    "cc": "rel_err_width_cc",
    "gg": "rel_err_width_gg",
    "gaga": "rel_err_width_gaga",
    "Zga": "rel_err_width_Zga",
    "total": "rel_err_width_total",
}


def parse_group_cols(spec: str) -> list[str]:
    cols=[c.strip() for c in spec.split(',') if c.strip()]
    if not cols:
        raise ValueError('group-cols must not be empty')
    return cols


def build_grouped_report(df: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    for c in group_cols:
        if c not in df.columns:
            raise ValueError(f"missing group column: {c}")

    if "triple_ok" in df.columns:
        filt = df["triple_ok"] == True
        work = df.loc[filt].copy()
    else:
        work = df.copy()

    present = {ch: col for ch, col in CHANNEL_COLS.items() if col in work.columns}

    gb = work.groupby(group_cols, dropna=False)
    rows=[]
    for keys, g in gb:
        if not isinstance(keys, tuple):
            keys=(keys,)
        row={group_cols[i]:keys[i] for i in range(len(group_cols))}
        row['n']=int(len(g))
        for ch,col in present.items():
            s=pd.to_numeric(g[col], errors='coerce').replace([np.inf,-np.inf],np.nan).dropna()
            row[f'mean_rel_{ch}']=float(s.mean()) if len(s) else np.nan
            row[f'median_rel_{ch}']=float(s.median()) if len(s) else np.nan
            row[f'max_rel_{ch}']=float(s.max()) if len(s) else np.nan
        rows.append(row)
    out=pd.DataFrame(rows)
    sort_cols=[c for c in group_cols]
    return out.sort_values(sort_cols).reset_index(drop=True)


def main()->int:
    p=argparse.ArgumentParser(description='Group merged comparison by physical configuration and summarize channel relative errors.')
    p.add_argument('csv', help='Merged comparison CSV path')
    p.add_argument('--group-cols', default='lambda1,lambda6,tan_beta', help='Comma-separated group columns')
    p.add_argument('--out-csv', required=True, help='Output CSV path')
    args=p.parse_args()

    csv_path=Path(args.csv)
    if not csv_path.exists():
        raise SystemExit(f'CSV not found: {csv_path}')
    group_cols=parse_group_cols(args.group_cols)
    df=pd.read_csv(csv_path)
    out=build_grouped_report(df, group_cols)

    out_csv=Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_csv,index=False)

    print('Grouped width-relative-error report')
    print(out.to_string(index=False))
    print(f"\nWrote grouped CSV: {out_csv}")
    return 0


if __name__=='__main__':
    raise SystemExit(main())
