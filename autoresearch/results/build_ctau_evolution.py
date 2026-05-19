#!/usr/bin/env python3
import csv
import json
from pathlib import Path
from datetime import datetime

ROOT = Path('/home/fabi/dihiggs')
RUNS_ROOT = ROOT / 'scripts' / 'out' / 'gemini_harness' / 'boulder_explore_ladder'
OUT_ROOT = ROOT / 'autoresearch' / 'results'
LOGS = OUT_ROOT / 'logs'
ALL_POINTS_CSV = LOGS / 'ctau_discoveries.csv'
FRONTIER_CSV = LOGS / 'ctau_frontier_events.csv'
PLOT_PNG = OUT_ROOT / 'ctau_evolution.png'


def parse_utc(s: str):
    if not s:
        return None
    s = s.strip()
    for fmt in ('%Y-%m-%dT%H:%M:%SZ', '%Y%m%dT%H%M%SZ'):
        try:
            return datetime.strptime(s, fmt)
        except ValueError:
            pass
    return None


def load_baseline_event():
    boulder_path = ROOT / 'boulder.json'
    if not boulder_path.exists():
        return None
    boulder = json.loads(boulder_path.read_text())
    best = boulder.get('best_known_point', {})
    if not best:
        return None

    baseline_quota = 100000
    recomputed = ROOT / 'scripts' / 'out' / 'gemini_harness' / 'orchestrate_10x10k' / 'final_summary_v3_recomputed.json'
    if recomputed.exists():
        try:
            js = json.loads(recomputed.read_text())
            baseline_quota = int(js.get('attempted_points') or js.get('planned_points') or baseline_quota)
        except Exception:
            pass

    return {
        'run_tag': 'baseline_orchestrate_10x10k',
        'utc': '',
        'stage': 'validated_baseline',
        'best_ctau_m': float(best.get('best_ctau_m')),
        'lambda1': best.get('lambda1'),
        'lambda6': best.get('lambda6'),
        'tan_beta': best.get('tan_beta'),
        'mA': best.get('mA'),
        'mphi': best.get('mphi'),
        'source_csv': best.get('source_csv', ''),
        'row_number': best.get('row_number', ''),
        'row_hash': best.get('row_hash', ''),
        'attempted_points_stage': baseline_quota,
        'attempted_points_run_total': baseline_quota,
    }


def run_stage_maps(stage_summary_path: Path):
    attempted_by_stage = {}
    cumulative_by_stage = {}
    total = 0
    with stage_summary_path.open(newline='') as f:
        rdr = csv.DictReader(f)
        for r in rdr:
            st = r.get('stage', '')
            att = int(float(r.get('attempted_points') or 0))
            attempted_by_stage[st] = att
            total += att
            cumulative_by_stage[st] = total
    return attempted_by_stage, cumulative_by_stage, total


def collect_discovery_events():
    events = []

    baseline = load_baseline_event()
    if baseline is not None:
        events.append(baseline)

    run_dirs = sorted([p for p in RUNS_ROOT.iterdir() if p.is_dir()])
    for run_dir in run_dirs:
        stage_summary = run_dir / 'stage_summary.csv'
        ledger = run_dir / 'best_points_ledger.jsonl'
        if not stage_summary.exists() or not ledger.exists():
            continue

        attempted_by_stage, cumulative_by_stage, run_total = run_stage_maps(stage_summary)

        with ledger.open() as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                row = json.loads(line)
                ctau = row.get('best_ctau_m')
                if ctau is None:
                    continue
                stage = row.get('stage', '')
                hp = row.get('best_hyperparameters') or {}
                events.append({
                    'run_tag': row.get('run_tag', run_dir.name),
                    'utc': row.get('utc', ''),
                    'stage': stage,
                    'best_ctau_m': float(ctau),
                    'lambda1': hp.get('lambda1'),
                    'lambda6': hp.get('lambda6'),
                    'tan_beta': hp.get('tan_beta'),
                    'mA': hp.get('mA'),
                    'mphi': hp.get('mphi'),
                    'source_csv': row.get('source_csv', ''),
                    'row_number': row.get('row_number', ''),
                    'row_hash': row.get('row_hash', ''),
                    'attempted_points_stage': attempted_by_stage.get(stage, ''),
                    'attempted_points_run_total': run_total,
                    'attempted_points_until_stage_in_run': cumulative_by_stage.get(stage, ''),
                    'accepted_attempted_ratio': row.get('accepted_attempted_ratio', ''),
                })

    # chronological order
    def sort_key(e):
        dt = parse_utc(e.get('utc', ''))
        if dt is None:
            # keep baseline first
            return datetime.min
        return dt

    events.sort(key=sort_key)

    # global cumulative quota
    global_quota = 0
    for e in events:
        if e.get('stage') == 'validated_baseline':
            global_quota = int(e.get('attempted_points_stage') or 0)
            e['attempted_points_global_cumulative'] = global_quota
            continue
        att = int(e.get('attempted_points_stage') or 0)
        global_quota += att
        e['attempted_points_global_cumulative'] = global_quota

    return events


def frontier_only(events):
    out = []
    best = float('-inf')
    for e in events:
        c = e['best_ctau_m']
        if c > best:
            best = c
            ee = dict(e)
            ee['frontier_ctau_m'] = best
            out.append(ee)
    return out


def write_csv(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text('')
        return

    # stable union of all keys in row order of appearance
    fields = []
    seen = set()
    for r in rows:
        for k in r.keys():
            if k not in seen:
                seen.add(k)
                fields.append(k)

    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        w.writeheader()
        for r in rows:
            w.writerow(r)


def make_plot(frontier_rows):
    try:
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f'[warn] matplotlib unavailable: {e}')
        return

    x = [int(r['attempted_points_global_cumulative']) for r in frontier_rows]
    y = [float(r['best_ctau_m']) for r in frontier_rows]
    # bubble size reflects stage quota used for discovery
    s = [max(40, int(r.get('attempted_points_stage') or 0) / 200) for r in frontier_rows]

    plt.figure(figsize=(12, 7))
    plt.plot(x, y, linewidth=2, alpha=0.8, label='frontier ctau')
    plt.scatter(x, y, s=s, alpha=0.55)

    for r, xi, yi in zip(frontier_rows, x, y):
        stage = r.get('stage', '')
        tag = r.get('run_tag', '')
        lbl = f"{tag}:{stage}" if stage != 'validated_baseline' else 'validated_baseline'
        plt.annotate(lbl, (xi, yi), textcoords='offset points', xytext=(5, 5), fontsize=7)

    plt.title('DiHiggs cτ frontier evolution vs cumulative attempted quota')
    plt.xlabel('Cumulative attempted points quota')
    plt.ylabel('Best cτ (m)')
    plt.grid(True, alpha=0.25)
    plt.tight_layout()
    PLOT_PNG.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(PLOT_PNG, dpi=150)
    plt.close()


def main():
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)

    events = collect_discovery_events()
    frontier = frontier_only(events)

    write_csv(ALL_POINTS_CSV, events)
    write_csv(FRONTIER_CSV, frontier)
    make_plot(frontier)

    print(f'events={len(events)}')
    print(f'frontier_events={len(frontier)}')
    print(f'all_points_csv={ALL_POINTS_CSV}')
    print(f'frontier_csv={FRONTIER_CSV}')
    print(f'plot_png={PLOT_PNG}')


if __name__ == '__main__':
    main()
