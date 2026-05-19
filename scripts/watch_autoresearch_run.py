from __future__ import annotations

import argparse
import json
import os
import signal
import sys
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from autoresearch_heartbeat import read_heartbeat


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace('+00:00', 'Z')


def parse_iso(ts: str | None) -> float | None:
    if not ts:
        return None
    try:
        return datetime.fromisoformat(ts.replace('Z', '+00:00')).timestamp()
    except Exception:
        return None


def pid_alive(pid: int | None) -> bool:
    if not pid or pid <= 0:
        return False
    try:
        os.kill(pid, 0)
        return True
    except OSError:
        return False


def load_env_file(path: Path) -> None:
    if not path.exists():
        return
    for line in path.read_text(encoding='utf-8').splitlines():
        line = line.strip()
        if not line or line.startswith('#') or '=' not in line:
            continue
        k, v = line.split('=', 1)
        os.environ.setdefault(k.strip(), v.strip())


def send_telegram(msg: str) -> tuple[bool, str]:
    tok = os.environ.get('TELEGRAM_BOT_TOKEN')
    cid = os.environ.get('TELEGRAM_CHAT_ID')
    if not tok or not cid:
        return False, 'missing TELEGRAM_BOT_TOKEN or TELEGRAM_CHAT_ID'
    url = f'https://api.telegram.org/bot{tok}/sendMessage'
    data = urllib.parse.urlencode({'chat_id': cid, 'text': msg}).encode('utf-8')
    try:
        req = urllib.request.Request(url, data=data, method='POST')
        with urllib.request.urlopen(req, timeout=15) as r:
            body = r.read().decode('utf-8', errors='replace')
        return True, body[:400]
    except Exception as e:
        return False, str(e)


def atomic_write(path: Path, payload: dict[str, Any]) -> None:
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    os.replace(tmp, path)


def determine_status(hb: dict[str, Any] | None, stale_minutes: int) -> tuple[str, dict[str, Any]]:
    if hb is None:
        return 'missing_heartbeat', {'reason': 'heartbeat.json not found'}

    hstatus = str(hb.get('status') or 'unknown')
    pid = hb.get('pid')
    pid_present = pid is not None and str(pid).strip() != ''
    updated_ts = parse_iso(hb.get('updated_at_utc'))
    now_ts = time.time()
    stale = updated_ts is None or (now_ts - updated_ts) > stale_minutes * 60
    alive = pid_alive(int(pid)) if pid_present else None

    if hstatus in {'completed', 'failed'}:
        return hstatus, {'stale': stale, 'pid_alive': alive}
    if hb.get('hard_stop_reason'):
        return 'failed', {'reason': 'hard_stop_reason set', 'hard_stop_reason': hb.get('hard_stop_reason')}
    if hstatus == 'running':
        if pid_present and not alive:
            return 'pid_dead_without_final', {'pid': pid}
        if stale:
            return 'stalled', {'reason': 'heartbeat stale while non-final running'}
        return 'running', {'pid': pid}
    if stale and alive:
        return 'stalled', {'reason': 'non-final status stale with alive pid'}
    if stale and not alive:
        return 'pid_dead_without_final', {'pid': pid}
    return hstatus, {'pid_alive': alive, 'stale': stale}


def should_notify(prev: str | None, curr: str, hb: dict[str, Any] | None) -> bool:
    if prev is None:
        return False
    if prev == curr:
        return False
    if prev == 'running' and curr in {'completed', 'failed', 'stalled'}:
        return True
    if curr == 'pid_dead_without_final':
        return True
    if hb and hb.get('contract_satisfied') is True:
        return True
    if hb and hb.get('hard_stop_reason'):
        return True
    return False


def run_once(run_dir: Path, stale_minutes: int, notify_test: bool = False) -> int:
    hb_path = run_dir / 'heartbeat.json'
    events_path = run_dir / 'watchdog_events.jsonl'
    status_path = run_dir / 'watchdog_status.json'

    hb = read_heartbeat(hb_path)
    prev = None
    if status_path.exists():
        try:
            prev = json.loads(status_path.read_text(encoding='utf-8')).get('status')
        except Exception:
            prev = None

    curr, details = determine_status(hb, stale_minutes)
    status = {
        'utc': now_iso(),
        'run_dir': str(run_dir),
        'status': curr,
        'previous_status': prev,
        'details': details,
        'heartbeat_present': hb is not None,
        'heartbeat_run_id': (hb or {}).get('run_id'),
        'heartbeat_updated_at_utc': (hb or {}).get('updated_at_utc'),
    }
    atomic_write(status_path, status)
    with events_path.open('a', encoding='utf-8') as h:
        h.write(json.dumps(status, sort_keys=True) + '\n')

    notify_msg = None
    if notify_test:
        notify_msg = f'[watchdog notify-test] run_dir={run_dir}'
    elif should_notify(prev, curr, hb):
        notify_msg = f"[watchdog] transition {prev} -> {curr} run_id={(hb or {}).get('run_id')} hard_stop={(hb or {}).get('hard_stop_reason')} contract_satisfied={(hb or {}).get('contract_satisfied')}"

    notify_result = None
    if notify_msg:
        ok, info = send_telegram(notify_msg)
        notify_result = {'sent': ok, 'info': info}
        with events_path.open('a', encoding='utf-8') as h:
            h.write(json.dumps({'utc': now_iso(), 'event': 'telegram_notify', 'message': notify_msg, 'result': notify_result}, sort_keys=True) + '\n')

    attempted = (hb or {}).get('attempted_points')
    best = (hb or {}).get('best_ctau_m')
    print(f"status={curr} prev={prev} attempted={attempted} best_ctau_m={best} updated={(hb or {}).get('updated_at_utc')}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description='Watchdog for autoresearch heartbeat')
    ap.add_argument('--run-dir', required=True)
    ap.add_argument('--once', action='store_true')
    ap.add_argument('--interval-seconds', type=int, default=300)
    ap.add_argument('--stale-minutes', type=int, default=15)
    ap.add_argument('--telegram-env', default=None)
    ap.add_argument('--notify-test', action='store_true')
    args = ap.parse_args()

    run_dir = Path(args.run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    if args.telegram_env:
        load_env_file(Path(args.telegram_env).expanduser())

    if args.once:
        return run_once(run_dir, args.stale_minutes, notify_test=args.notify_test)

    while True:
        rc = run_once(run_dir, args.stale_minutes, notify_test=args.notify_test)
        if rc != 0:
            return rc
        time.sleep(max(1, args.interval_seconds))


if __name__ == '__main__':
    raise SystemExit(main())
