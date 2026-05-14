from __future__ import annotations

import argparse
import json
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
import sys
import uuid

try:
    from scripts.gemini_contract_templates import materialize_template
except ModuleNotFoundError:
    sys.path.insert(0, '/home/fabi/dihiggs/scripts')
    from gemini_contract_templates import materialize_template

ROOT = Path('/home/fabi/dihiggs')
HARNESS_DIR = ROOT / 'scripts' / 'out' / 'gemini_harness'

SCHEMA_ACTIONS = {"run_template", "run_contract", "summarize_only", "stop"}
FORBIDDEN_KEYWORDS = {
    "edit_cpp", "edit_scoring", "edit_triple_ok", "delete_artifacts", "widen_envelope", "broad_scan",
    "rm -rf", "sudo", "bash -c", "--force",
}


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace('+00:00', 'Z')


def load_json(path: Path) -> dict[str, Any]:
    obj = json.loads(path.read_text(encoding='utf-8'))
    if not isinstance(obj, dict):
        raise ValueError(f'{path} must contain a JSON object')
    return obj


def _read_tail(path: Path, max_lines: int = 80) -> str:
    if not path.exists():
        return ''
    lines = path.read_text(encoding='utf-8', errors='replace').splitlines()
    return '\n'.join(lines[-max_lines:])


def build_prompt(config: dict[str, Any], worker_state: dict[str, Any]) -> str:
    best_ctau = worker_state.get('best_ctau_m', config['baseline_best_ctau_m'])
    remaining = max(0, int(config['max_points_total']) - int(worker_state.get('executed_points_total', 0)))
    cycle_id = worker_state.get('cycle_id', 'unknown_cycle')
    last_ctx = worker_state.get('last_iteration_summary_for_prompt')
    schema = {
        "action": "run_template | summarize_only | stop",
        "template_name": "ctau_mphi_mA_refine|mphi_lower_edge_probe|mA_upper_edge_probe|control_box",
        "template_overrides": {
            "mphi_min": 180, "mphi_max": 220, "n_mphi": 40,
            "mA_values": [450, 475, 500],
            "tan_beta_center": 124416,
            "lambda6_center": 0.0019683,
            "tan_beta_local_multipliers": [0.98, 1.0, 1.02],
            "lambda6_local_multipliers": [0.97, 1.0, 1.03],
        },
        "rationale": "string",
        "safety_checks": ["..."],
        "stop_conditions": ["..."],
    }
    return (
        "You are Gemini Worker in bounded mode. Output JSON only. No markdown.\n"
        f"Cycle ID: {cycle_id}\n"
        f"Current best ctau_m: {best_ctau}\n"
        f"Target ctau_m: {config['target_ctau_m']}\n"
        f"Remaining executed-point budget: {remaining}\n"
        f"Envelope: {json.dumps(config['envelope'])}\n"
        f"Allowed actions: {config['allowed_actions']}\n"
        "Do NOT output full contracts; use run_template with template_overrides only.\n"
        "Keep lambda1 fixed at 1.0 and n_lam1=1.\n"
        f"Last uncontaminated iteration summary (if any): {json.dumps(last_ctx) if isinstance(last_ctx, dict) else 'none'}\n"
        f"Schema: {json.dumps(schema)}\n"
        f"Context GEMINI.md tail:\n{_read_tail(ROOT/'GEMINI.md', 120)}\n"
    )


def build_repair_prompt(config: dict[str, Any], worker_state: dict[str, Any], plancheck: dict[str, Any]) -> str:
    return (
        "Plancheck failed. Choose one valid template only and output corrected JSON only.\n"
        f"failure_reason: {', '.join(plancheck.get('plancheck_rejection_reasons', []))}\n"
        f"min_points_per_subrun: {config['min_points_per_subrun']}\n"
        f"max_points_per_subrun: {config['max_points_per_subrun']}\n"
        f"planned_points: {plancheck.get('planned_points')}\n"
        "Do not write full contract manually.\n"
        "Increase points using n_mphi or allowed local brackets.\n"
        "Keep lambda1 fixed at 1.0 (n_lam1=1).\n"
        "Forbidden: envelope widening, broad scan, C++ edits, scoring/triple_ok edits.\n"
        "Templates: ctau_mphi_mA_refine, mphi_lower_edge_probe, mA_upper_edge_probe, control_box.\n"
    )


def call_gemini(prompt: str) -> tuple[int, str, str]:
    proc = subprocess.run(['gemini', '-p', prompt], text=True, capture_output=True, cwd=str(ROOT))
    return proc.returncode, proc.stdout, proc.stderr


def parse_json_response(text: str) -> tuple[dict[str, Any] | None, str | None]:
    try:
        obj = json.loads(text.strip())
    except Exception as exc:
        return None, f'invalid_json: {exc}'
    if not isinstance(obj, dict):
        return None, 'invalid_json: top-level not object'
    return obj, None


def validate_gemini_json(obj: dict[str, Any], config: dict[str, Any], worker_state: dict[str, Any]) -> dict[str, Any]:
    errors: list[str] = []
    action = obj.get('action')
    if action not in SCHEMA_ACTIONS or action not in set(config['allowed_actions']):
        errors.append('action not allowlisted')
    if action == 'run_contract' and not bool(config.get('allow_direct_run_contract', False)):
        errors.append('run_contract disabled; use run_template')

    txt = json.dumps(obj).lower()
    for kw in FORBIDDEN_KEYWORDS.union(set(config['forbidden_actions'])):
        if kw.lower() in txt:
            errors.append(f'forbidden keyword/action detected: {kw}')
            break

    required = ['action', 'rationale', 'safety_checks', 'stop_conditions']
    for k in required:
        if k not in obj:
            errors.append(f'missing key: {k}')

    expected_points = None
    if action == 'run_template':
        if 'template_name' not in obj:
            errors.append('missing key: template_name')
        if 'template_overrides' in obj and not isinstance(obj.get('template_overrides'), dict):
            errors.append('template_overrides must be object')

    if action in {'run_template', 'run_contract'}:
        remaining = int(config['max_points_total']) - int(worker_state.get('executed_points_total', 0))
        if remaining <= 0:
            errors.append('no remaining budget')

    return {
        'valid': len(errors) == 0,
        'errors': errors,
        'action': action,
        'expected_points': expected_points,
        'json_schema_valid': not any(e.startswith('missing key') for e in errors),
        'envelope_valid': True,
        'policy_valid': not any('forbidden' in e or 'disabled' in e for e in errors),
    }


def append_jsonl(path: Path, rec: dict[str, Any]) -> None:
    with path.open('a', encoding='utf-8') as h:
        h.write(json.dumps(rec, sort_keys=True) + '\n')


def write_worker_progress(state: dict[str, Any], out_dir: Path) -> None:
    lines = ['# Gemini worker progress', '', f"Generated: {now_iso()}", f"Cycle ID: {state.get('cycle_id')}", f"Status: {state.get('status')}", f"Best ctau_m: {state.get('best_ctau_m')}"]
    for k in ('proposed_points_total', 'plancheck_points_total', 'executed_points_total', 'accepted_physics_points_total', 'isolated_summary_rows_total', 'contaminated_summary_count', 'rejected_points_total'):
        lines.append(f"{k}: {state.get(k, 0)}")
    lines.append(f"Iterations completed: {state.get('iterations_completed')}")
    lines.append('')
    for ev in state.get('history', [])[-20:]:
        lines.append(f"- iter {ev.get('iteration')}: {ev.get('result')} action={ev.get('action')} valid={ev.get('valid_json')}")
    (out_dir / 'worker_progress.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')


def _extract_plancheck_summary(plan_dir: Path) -> dict[str, Any]:
    reasons: list[str] = []
    planned_subcampaigns = 0
    planned_points = 0
    insufficient_points = False
    sub_path = plan_dir / 'subcampaigns.jsonl'
    if sub_path.exists():
        for line in sub_path.read_text(encoding='utf-8', errors='replace').splitlines():
            if not line.strip():
                continue
            try:
                rec = json.loads(line)
            except Exception:
                continue
            planned_subcampaigns += 1
            planned_points += int(rec.get('estimated_raw_points', 0) or 0)
            if rec.get('reason') == 'insufficient_points' or rec.get('status') == 'rejected':
                insufficient_points = True
    st_path = plan_dir / 'adaptive_state.json'
    if st_path.exists():
        try:
            st = load_json(st_path)
            if st.get('reason') == 'insufficient_points':
                insufficient_points = True
            for item in st.get('iterations', []):
                if isinstance(item, dict) and (item.get('reason') == 'insufficient_points' or item.get('status') == 'rejected'):
                    insufficient_points = True
            if planned_points == 0:
                planned_points = int(st.get('total_points_executed', 0) or 0)
        except Exception:
            reasons.append('failed_to_parse_adaptive_state')
    return {'planned_subcampaigns': planned_subcampaigns, 'planned_points': planned_points, 'insufficient_points': insufficient_points, 'reasons': reasons}


def _run_plancheck(contract_path: Path, plan_dir: Path) -> tuple[bool, dict[str, Any], str]:
    if plan_dir.exists():
        for p in plan_dir.glob("**/*"):
            if p.is_file():
                p.unlink()
    plan_dir.mkdir(parents=True, exist_ok=True)
    cmd = ['python', '-m', 'autoresearch.harness.bounded_adaptive_search', '--contract', str(contract_path), '--output-dir', str(plan_dir), '--plan-only']
    proc = subprocess.run(cmd, cwd=str(ROOT), text=True, capture_output=True)
    summary = _extract_plancheck_summary(plan_dir)
    return proc.returncode == 0, summary, proc.stderr[-2000:]


def _budget_snapshot(cfg: dict[str, Any], state: dict[str, Any]) -> dict[str, int]:
    before = int(state.get('executed_points_total', 0))
    max_total = int(cfg['max_points_total'])
    return {'executed_before': before, 'remaining_before': max(0, max_total - before), 'max_total': max_total}


def _count_executed_points(run_dir: Path) -> int:
    st = run_dir / 'adaptive_state.json'
    if not st.exists():
        return 0
    try:
        return int(load_json(st).get('total_points_executed', 0) or 0)
    except Exception:
        return 0


def _read_iteration_isolated_rows(run_dir: Path) -> int:
    st = run_dir / 'adaptive_state.json'
    if not st.exists():
        return 0
    try:
        payload = load_json(st)
        iters = payload.get('iterations', [])
        if not isinstance(iters, list) or not iters:
            return 0
        last = iters[-1]
        if not isinstance(last, dict):
            return 0
        summ = last.get('summary', {})
        if not isinstance(summ, dict):
            return 0
        iso = summ.get('isolated_summary', {})
        if not isinstance(iso, dict):
            return 0
        return int(iso.get('isolated_n_rows', 0) or 0)
    except Exception:
        return 0


def _materialize_from_action(obj: dict[str, Any], cfg: dict[str, Any], it: int, cycle_id: str, iteration_id: str) -> tuple[dict[str, Any], int]:
    env = cfg['envelope']
    if obj.get('action') == 'run_template':
        tm = materialize_template(
            template_name=str(obj.get('template_name')),
            overrides=obj.get('template_overrides') if isinstance(obj.get('template_overrides'), dict) else {},
            envelope=env,
            min_points=int(cfg['min_points_per_subrun']),
            max_points=int(cfg['max_points_per_subrun']),
            iteration=it,
            max_total_points=int(cfg['max_points_total']),
            tan_beta_center_default=float(cfg.get('tan_beta_center', 124416.0)),
            lambda6_center_default=float(cfg.get('lambda6_center', 0.0019683)),
            cycle_id=cycle_id,
            iteration_id=iteration_id,
        )
        return tm.contract, int(tm.real_points)
    raise ValueError('run_contract direct path disabled in v1.2 unless explicitly enabled')


def run_worker(config_path: Path, *, dry_run: bool, iterations_override: int | None) -> int:
    cfg = load_json(config_path)
    out_dir = HARNESS_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / 'candidates').mkdir(parents=True, exist_ok=True)

    state_path = out_dir / 'worker_state.json'
    events_path = out_dir / 'worker_events.jsonl'
    cycle_id = f"cycle_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}_{uuid.uuid4().hex[:6]}"
    state = {
        'status': 'running', 'started_utc': now_iso(), 'best_ctau_m': cfg['baseline_best_ctau_m'], 'points_used': 0,
        'proposed_points_total': 0, 'plancheck_points_total': 0, 'executed_points_total': 0,
        'accepted_physics_points_total': 0, 'rejected_points_total': 0,
        'isolated_summary_rows_total': 0, 'contaminated_summary_count': 0,
        'iterations_completed': 0, 'invalid_json_count': 0, 'no_improvement_count': 0, 'history': [],
        'cycle_id': cycle_id,
        'last_iteration_summary_for_prompt': None,
    }

    for k in ('proposed_points_total', 'plancheck_points_total', 'executed_points_total', 'accepted_physics_points_total', 'rejected_points_total', 'isolated_summary_rows_total', 'contaminated_summary_count'):
        state.setdefault(k, 0)
    state.setdefault('cycle_id', cycle_id)
    state.setdefault('last_iteration_summary_for_prompt', None)

    max_iter = int(iterations_override or cfg['max_iterations'])
    deadline = time.time() + int(cfg['max_walltime_minutes']) * 60

    for it in range(1, max_iter + 1):
        iteration_id = f"iter{it:02d}_{uuid.uuid4().hex[:6]}"
        if time.time() > deadline:
            state['status'] = 'stopped_walltime'; break

        prompt = build_prompt(cfg, state)
        (out_dir / 'latest_gemini_prompt.md').write_text(prompt, encoding='utf-8')
        _, out, _ = call_gemini(prompt)
        obj, parse_err = parse_json_response(out)
        if parse_err:
            state['invalid_json_count'] = int(state.get('invalid_json_count', 0)) + 1
            strict = prompt + '\nSTRICT: JSON object only.'
            _, out2, _ = call_gemini(strict)
            obj2, parse_err2 = parse_json_response(out2)
            if parse_err2 is not None:
                report = {'valid': False, 'errors': [parse_err], 'json_schema_valid': False, 'envelope_valid': False, 'policy_valid': False,
                          'plancheck_valid': False, 'planned_subcampaigns': 0, 'planned_points': 0,
                          'plancheck_rejection_reasons': ['invalid_json'], 'budget_before': _budget_snapshot(cfg, state),
                          'budget_after_if_executed': _budget_snapshot(cfg, state), 'execution_allowed': False}
                (out_dir / 'latest_validation_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
                state['history'].append({'iteration': it, 'result': 'invalid_json', 'action': None, 'valid_json': False})
                state['status'] = 'stopped_invalid_json_twice'
                break
            obj = obj2

        assert obj is not None
        (out_dir / 'latest_gemini_response.json').write_text(json.dumps(obj, indent=2), encoding='utf-8')
        validation = validate_gemini_json(obj, cfg, state)
        if not validation['valid']:
            report = {**validation, 'plancheck_valid': False, 'planned_subcampaigns': 0, 'planned_points': 0,
                      'plancheck_rejection_reasons': ['json_or_policy_validation_failed'], 'budget_before': _budget_snapshot(cfg, state),
                      'budget_after_if_executed': _budget_snapshot(cfg, state), 'execution_allowed': False}
            (out_dir / 'latest_validation_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
            state['history'].append({'iteration': it, 'result': 'validation_failed', 'action': obj.get('action'), 'valid_json': True})
            state['status'] = 'stopped_validation_failed'
            break

        action = obj.get('action')
        if action in {'stop', 'summarize_only'}:
            report = {**validation, 'plancheck_valid': False, 'planned_subcampaigns': 0, 'planned_points': 0,
                      'plancheck_rejection_reasons': [], 'budget_before': _budget_snapshot(cfg, state),
                      'budget_after_if_executed': _budget_snapshot(cfg, state), 'execution_allowed': False}
            (out_dir / 'latest_validation_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
            state['history'].append({'iteration': it, 'result': 'proposal_only' if action == 'summarize_only' else 'worker_stop', 'action': action, 'valid_json': True})
            state['iterations_completed'] += 1
            state['status'] = 'stopped_by_worker' if action == 'stop' else state.get('status', 'running')
            if action == 'stop':
                break
            continue

        try:
            contract, proposed_points = _materialize_from_action(obj, cfg, it, str(state.get('cycle_id')), iteration_id)
        except Exception as exc:
            report = {**validation, 'plancheck_valid': False, 'planned_subcampaigns': 0, 'planned_points': 0,
                      'plancheck_rejection_reasons': [f'template_materialization_failed:{exc}'], 'budget_before': _budget_snapshot(cfg, state),
                      'budget_after_if_executed': _budget_snapshot(cfg, state), 'execution_allowed': False}
            (out_dir / 'latest_validation_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
            state['history'].append({'iteration': it, 'result': 'plan_rejected', 'action': action, 'valid_json': True})
            state['status'] = 'stopped_plan_rejected'
            break

        state['proposed_points_total'] = int(state.get('proposed_points_total', 0)) + max(0, int(proposed_points))
        runtime = contract.setdefault('runtime', {})
        base_campaign = str(runtime.get('campaign', f'gemini_worker_iter{it:02d}'))
        runtime['campaign'] = f"{base_campaign}_{state.get('cycle_id')}_{iteration_id}"
        cand = out_dir / 'candidates' / f'iter{it:02d}_{iteration_id}_candidate_contract.json'
        cand.write_text(json.dumps(contract, indent=2), encoding='utf-8')

        plan_dir = out_dir / f'plancheck_{state.get("cycle_id")}_{iteration_id}'
        ok, ps, _ = _run_plancheck(cand, plan_dir)
        planned_subcampaigns = int(ps.get('planned_subcampaigns', 0))
        planned_points = int(ps.get('planned_points', 0))
        state['plancheck_points_total'] = int(state.get('plancheck_points_total', 0)) + max(0, planned_points)

        template_real_points = int(proposed_points)
        contract_estimated_points = int(contract.get('strategy', {}).get('expected_raw_points', contract.get('template_metadata', {}).get('template_real_points', template_real_points)))

        budget_before = _budget_snapshot(cfg, state)
        budget_after = {'executed_before': budget_before['executed_before'], 'remaining_before': budget_before['remaining_before'],
                        'executed_after': budget_before['executed_before'] + max(0, planned_points),
                        'remaining_after': max(0, budget_before['max_total'] - (budget_before['executed_before'] + max(0, planned_points))),
                        'max_total': budget_before['max_total']}

        reasons = list(ps.get('reasons', []))
        if not ok: reasons.append('plan_only_command_failed')
        if planned_subcampaigns <= 0: reasons.append('planned_subcampaigns<=0')
        if planned_points < int(cfg['min_points_per_subrun']): reasons.append('planned_points_below_min')
        if planned_points > int(cfg['max_points_per_subrun']): reasons.append('planned_points_above_max')
        if ps.get('insufficient_points'): reasons.append('insufficient_points_rejection')
        if not (template_real_points == contract_estimated_points == planned_points): reasons.append('template_plan_mismatch')
        if budget_before['executed_before'] + planned_points > int(cfg['max_points_total']): reasons.append('planned_points_exceed_remaining_budget')

        report = {**validation, 'plancheck_valid': len(reasons) == 0, 'planned_subcampaigns': planned_subcampaigns, 'planned_points': planned_points,
                  'template_real_points': template_real_points, 'contract_estimated_points': contract_estimated_points,
                  'plancheck_rejection_reasons': sorted(set(reasons)), 'budget_before': budget_before,
                  'budget_after_if_executed': budget_after, 'execution_allowed': len(reasons) == 0 and (not dry_run)}

        if not report['plancheck_valid']:
            state['rejected_points_total'] = int(state.get('rejected_points_total', 0)) + max(0, planned_points)
            append_jsonl(events_path, {'timestamp': now_iso(), 'iteration': it, 'event': 'plan_rejected', 'planned_points': planned_points, 'reasons': report['plancheck_rejection_reasons']})
            repair_prompt = build_repair_prompt(cfg, state, report)
            (out_dir / 'latest_gemini_repair_prompt.md').write_text(repair_prompt, encoding='utf-8')
            _, out_r, _ = call_gemini(repair_prompt)
            obj_r, err_r = parse_json_response(out_r)
            report['repair_attempted'] = True
            if err_r is not None or not isinstance(obj_r, dict):
                report['plancheck_rejection_reasons'] = sorted(set(report['plancheck_rejection_reasons'] + ['repair_invalid_json']))
            else:
                report['plancheck_rejection_reasons'] = sorted(set(report['plancheck_rejection_reasons'] + ['repair_json_invalid_or_non_run_template']))
            (out_dir / 'latest_validation_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
            state['history'].append({'iteration': it, 'result': 'plan_rejected', 'action': action, 'valid_json': True})
            state['status'] = 'stopped_plan_rejected'
            break

        run_result = 'dry_run_skipped_execution'
        if not dry_run:
            run_dir = out_dir / f'worker_run_{state.get("cycle_id")}_{iteration_id}'
            cpath = out_dir / f'worker_contract_iter{it:02d}_{iteration_id}.json'
            cpath.write_text(json.dumps(contract, indent=2), encoding='utf-8')
            cmd = ['python', '-m', 'autoresearch.harness.bounded_adaptive_search', '--contract', str(cpath), '--output-dir', str(run_dir), '--execute', '--max-iterations', '1']
            proc = subprocess.run(cmd, cwd=str(ROOT), text=True, capture_output=True)
            run_result = 'executed' if proc.returncode == 0 else 'execution_failed'
            append_jsonl(events_path, {'timestamp': now_iso(), 'iteration': it, 'iteration_id': iteration_id, 'cycle_id': state.get('cycle_id'), 'event': 'run_contract', 'return_code': proc.returncode, 'stdout_tail': proc.stdout[-2000:], 'stderr_tail': proc.stderr[-2000:]})
            executed_points = _count_executed_points(run_dir)
            isolated_rows = _read_iteration_isolated_rows(run_dir)
            contaminated = isolated_rows > executed_points if executed_points > 0 else False
            state['isolated_summary_rows_total'] = int(state.get('isolated_summary_rows_total', 0)) + max(0, isolated_rows)
            state['executed_points_total'] = int(state.get('executed_points_total', 0)) + max(0, executed_points)
            state['accepted_physics_points_total'] = int(state.get('accepted_physics_points_total', 0)) + max(0, executed_points)
            state['points_used'] = int(state.get('executed_points_total', 0))
            if contaminated:
                state['contaminated_summary_count'] = int(state.get('contaminated_summary_count', 0)) + 1
                state['last_iteration_summary_for_prompt'] = None
                report['contaminated_summary'] = True
                report['isolated_summary_rows'] = isolated_rows
                report['executed_points_iteration'] = executed_points
                report['contamination_reason'] = 'isolated_rows_exceed_executed_points'
                state['history'].append({'iteration': it, 'iteration_id': iteration_id, 'result': 'contaminated_summary', 'action': action, 'valid_json': True})
                state['status'] = 'stopped_contamination_detected'
                (out_dir / 'latest_validation_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
                break
            report['contaminated_summary'] = False
            report['isolated_summary_rows'] = isolated_rows
            report['executed_points_iteration'] = executed_points
            state['last_iteration_summary_for_prompt'] = {
                'iteration': it,
                'iteration_id': iteration_id,
                'executed_points': executed_points,
                'isolated_summary_rows': isolated_rows,
                'planned_points': planned_points,
                'template_name': obj.get('template_name'),
            }

        (out_dir / 'latest_validation_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
        state['iterations_completed'] += 1
        state['history'].append({'iteration': it, 'result': run_result, 'action': action, 'valid_json': True})
        if dry_run:
            append_jsonl(events_path, {'timestamp': now_iso(), 'iteration': it, 'event': 'dry_run_validated', 'planned_points': planned_points})

    if state.get('status') == 'running':
        state['status'] = 'completed_loop'
    state['updated_utc'] = now_iso()
    state_path.write_text(json.dumps(state, indent=2), encoding='utf-8')
    write_worker_progress(state, out_dir)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description='Gemini bounded worker loop wrapper')
    ap.add_argument('--config', required=True)
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--iterations', type=int)
    args = ap.parse_args()
    return run_worker(Path(args.config), dry_run=bool(args.dry_run), iterations_override=args.iterations)


if __name__ == '__main__':
    raise SystemExit(main())
