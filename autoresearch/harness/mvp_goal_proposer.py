from __future__ import annotations

from collections.abc import Mapping, Sequence
from itertools import product
from typing import cast

from autoresearch.benchmarks.score import compute_coverage
from autoresearch.harness.coverage_contract_bridge import coverage_history_from_state_records
from autoresearch.harness.mvp_graph import MvpGraphState


def _mapping(value: object) -> Mapping[str, object]:
    if isinstance(value, Mapping):
        return cast(Mapping[str, object], value)
    return {}


def _extract_axes_binned(record: Mapping[str, object]) -> dict[str, object]:
    direct_axes = _mapping(record.get("axes_binned"))
    if direct_axes:
        return dict(direct_axes)
    payload_axes = _mapping(_mapping(record.get("payload")).get("axes_binned"))
    if payload_axes:
        return dict(payload_axes)
    return {}


def _identity_ref(record: Mapping[str, object], *, state_type: str) -> dict[str, object] | None:
    campaign_id = record.get("campaign_id")
    arm_id = record.get("arm_id")
    run_dir_fingerprint = record.get("run_dir_fingerprint")
    attempt_id = record.get("attempt_id")
    if campaign_id in (None, "") or arm_id in (None, "") or run_dir_fingerprint in (None, "") or attempt_id in (None, ""):
        return None
    return {
        "state_type": state_type,
        "campaign_id": str(campaign_id),
        "arm_id": str(arm_id),
        "run_dir_fingerprint": str(run_dir_fingerprint),
        "attempt_id": str(attempt_id),
    }


def _graph_refs_for_attempts(graph: MvpGraphState | None, attempt_ids: Sequence[str]) -> list[dict[str, object]]:
    if graph is None or not attempt_ids:
        return []

    attempt_id_set = set(attempt_ids)
    refs: list[dict[str, object]] = []

    for key, record in sorted(graph.nodes.items()):
        attempt_id = record.get("attempt_id")
        if isinstance(attempt_id, str) and attempt_id in attempt_id_set:
            refs.append(
                {
                    "state_type": "graph_node",
                    "entity_type": key[0],
                    "campaign_id": str(record.get("campaign_id", "")),
                    "arm_id": str(record.get("arm_id", "")),
                    "run_dir_fingerprint": str(record.get("run_dir_fingerprint", "")),
                    "attempt_id": attempt_id,
                }
            )

    for key, record in sorted(graph.edges.items()):
        source_attempt_id = record.get("source_attempt_id")
        target_attempt_id = record.get("target_attempt_id")
        attempt_match = (
            isinstance(source_attempt_id, str)
            and source_attempt_id in attempt_id_set
            or isinstance(target_attempt_id, str)
            and target_attempt_id in attempt_id_set
            or isinstance(record.get("attempt_id"), str)
            and str(record.get("attempt_id")) in attempt_id_set
        )
        if attempt_match:
            refs.append(
                {
                    "state_type": "graph_edge",
                    "entity_type": key[0],
                    "campaign_id": str(record.get("campaign_id", "")),
                    "arm_id": str(record.get("arm_id", "")),
                    "source_attempt_id": str(record.get("source_attempt_id", "")),
                    "target_attempt_id": str(record.get("target_attempt_id", "")),
                }
            )

    return refs


def _axis_specs(config: Mapping[str, object]) -> list[Mapping[str, object]]:
    multi_axis = _mapping(_mapping(config.get("metrics")).get("multi_axis"))
    coverage_axes = multi_axis.get("coverage_axes")
    if isinstance(coverage_axes, Sequence) and not isinstance(coverage_axes, (str, bytes)):
        specs: list[Mapping[str, object]] = [cast(Mapping[str, object], spec) for spec in coverage_axes if isinstance(spec, Mapping)]
        if specs:
            return specs

    collapse_axes = multi_axis.get("collapse_axes")
    if not isinstance(collapse_axes, Sequence) or isinstance(collapse_axes, (str, bytes)):
        collapse_axes = ("tb", "lam1_bin")

    return [{"name": str(axis_name)} for axis_name in collapse_axes]


def _axis_values(axis_spec: Mapping[str, object], config: Mapping[str, object]) -> list[object]:
    axis_name = str(axis_spec.get("name", "")).strip()
    kind = str(axis_spec.get("kind", "")).strip()
    if kind == "categorical":
        domain = axis_spec.get("domain")
        if isinstance(domain, Sequence) and not isinstance(domain, (str, bytes)):
            return list(domain)

    if kind == "discrete":
        domain_size = axis_spec.get("domain_size", 0)
        if isinstance(domain_size, int) and domain_size > 0:
            return list(range(domain_size))

    search = _mapping(config.get("search"))
    if axis_name == "tb":
        tb_values = search.get("tb_values")
        if isinstance(tb_values, Sequence) and not isinstance(tb_values, (str, bytes)):
            return list(tb_values)
    if axis_name == "lam1_bin":
        n_bins = _mapping(search.get("lam1")).get("n_bins")
        if isinstance(n_bins, int) and n_bins > 0:
            return list(range(n_bins))
    if axis_name == "mphi_bin":
        n_bins = _mapping(search.get("mphi")).get("n_bins")
        if isinstance(n_bins, int) and n_bins > 0:
            return list(range(n_bins))

    return []


def _cell_id(axes_binned: Mapping[str, object], axis_order: Sequence[str]) -> str:
    return "|".join(f"{axis_name}={axes_binned[axis_name]}" for axis_name in axis_order)


def _covered_cell_ids(coverage_state: Sequence[Mapping[str, object]], axis_order: Sequence[str]) -> set[str]:
    covered: set[str] = set()
    for row in coverage_state:
        if str(row.get("artifact_status", "complete")) == "partial":
            continue
        axes_binned = _extract_axes_binned(row)
        if all(axis_name in axes_binned for axis_name in axis_order):
            covered.add(_cell_id(axes_binned, axis_order))
    return covered


def _scheduled_cell_ids(scheduled_goals: Sequence[Mapping[str, object]], axis_order: Sequence[str]) -> tuple[set[str], int]:
    scheduled: set[str] = set()
    suppressed_duplicates = 0
    for goal in scheduled_goals:
        axes_binned = _extract_axes_binned(goal)
        if not all(axis_name in axes_binned for axis_name in axis_order):
            continue
        cell_id = _cell_id(axes_binned, axis_order)
        if cell_id in scheduled:
            suppressed_duplicates += 1
            continue
        scheduled.add(cell_id)
    return scheduled, suppressed_duplicates


def _frontier_support(
    candidate_axes: Mapping[str, object],
    coverage_state: Sequence[Mapping[str, object]],
    discovery_records: Sequence[Mapping[str, object]],
    *,
    axis_specs: Sequence[Mapping[str, object]],
    axis_order: Sequence[str],
    graph: MvpGraphState | None,
) -> tuple[int, list[dict[str, object]], list[str]]:
    discovery_ref_by_attempt_id = {
        str(record.get("attempt_id")): ref
        for record in discovery_records
        for ref in [_identity_ref(record, state_type="discovery_records")]
        if ref is not None and isinstance(record.get("attempt_id"), str)
    }
    supporting_refs: list[dict[str, object]] = []
    supporting_attempt_ids: list[str] = []
    supporting_cell_ids: list[str] = []

    for record in coverage_state:
        if str(record.get("artifact_status", "complete")) == "partial":
            continue
        axes_binned = _extract_axes_binned(record)
        if not all(axis_name in axes_binned for axis_name in axis_order):
            continue

        distance = 0
        valid_neighbor = True
        for axis_spec in axis_specs:
            axis_name = str(axis_spec.get("name", ""))
            kind = str(axis_spec.get("kind", ""))
            candidate_value = candidate_axes.get(axis_name)
            existing_value = axes_binned.get(axis_name)
            if kind == "discrete":
                if not isinstance(candidate_value, int) or not isinstance(existing_value, int):
                    valid_neighbor = False
                    break
                distance += abs(candidate_value - existing_value)
            else:
                if candidate_value != existing_value:
                    valid_neighbor = False
                    break

        if not valid_neighbor or distance != 1:
            continue

        source_ref = _identity_ref(record, state_type="coverage_state")
        if source_ref is not None:
            supporting_refs.append(source_ref)
        attempt_id = record.get("attempt_id")
        if isinstance(attempt_id, str) and attempt_id:
            supporting_attempt_ids.append(attempt_id)
            discovery_ref = discovery_ref_by_attempt_id.get(attempt_id)
            if discovery_ref is not None:
                supporting_refs.append(discovery_ref)
        supporting_cell_ids.append(_cell_id(axes_binned, axis_order))

    supporting_refs.extend(_graph_refs_for_attempts(graph, supporting_attempt_ids))
    return len(supporting_cell_ids), supporting_refs, sorted(set(supporting_cell_ids))


def _latest_coverage_ref(coverage_state: Sequence[Mapping[str, object]]) -> dict[str, object] | None:
    for row in reversed(list(coverage_state)):
        ref = _identity_ref(row, state_type="coverage_state")
        if ref is not None:
            return ref
    return None


def propose_parameter_space_goals(
    *,
    coverage_state: Sequence[Mapping[str, object]],
    discovery_records: Sequence[Mapping[str, object]],
    graph: MvpGraphState | None,
    config: Mapping[str, object],
    scheduled_goals: Sequence[Mapping[str, object]] = (),
    max_new_runs: int | None = None,
) -> dict[str, object]:
    axis_specs = _axis_specs(config)
    axis_order = [str(spec.get("name", "")).strip() for spec in axis_specs if str(spec.get("name", "")).strip()]
    if not axis_order:
        raise ValueError("No proposer axes configured")

    axis_domains = {axis_name: _axis_values(spec, config) for axis_name, spec in zip(axis_order, axis_specs, strict=True)}
    if any(not domain for domain in axis_domains.values()):
        missing_axes = sorted(axis_name for axis_name, domain in axis_domains.items() if not domain)
        raise ValueError(f"Missing proposer domain values for axes: {', '.join(missing_axes)}")

    if max_new_runs is None:
        limits = _mapping(config.get("limits"))
        configured = limits.get("max_new_run_dirs_per_round")
        max_new_runs = int(configured) if isinstance(configured, int) else 3
    max_new_runs = max(0, int(max_new_runs))

    coverage_history = coverage_history_from_state_records(coverage_state)
    current_coverage = compute_coverage(coverage_history, config)
    visited_per_axis: dict[str, set[object]] = {axis_name: set() for axis_name in axis_order}
    for item in coverage_history:
        axes_binned = _extract_axes_binned(item)
        for axis_name in axis_order:
            if axis_name in axes_binned:
                visited_per_axis[axis_name].add(axes_binned[axis_name])

    covered_cells = _covered_cell_ids(coverage_state, axis_order)
    scheduled_cells, duplicate_scheduled_goals = _scheduled_cell_ids(scheduled_goals, axis_order)
    latest_coverage_ref = _latest_coverage_ref(coverage_state)

    suppressed_as_covered = 0
    suppressed_as_scheduled = 0
    proposals: list[dict[str, object]] = []
    candidate_slots = product(*(axis_domains[axis_name] for axis_name in axis_order))

    for axis_values in candidate_slots:
        axes_binned = dict(zip(axis_order, axis_values, strict=True))
        cell_id = _cell_id(axes_binned, axis_order)
        if cell_id in covered_cells:
            suppressed_as_covered += 1
            continue
        if cell_id in scheduled_cells:
            suppressed_as_scheduled += 1
            continue

        coverage_after = compute_coverage([*coverage_history, {"axes_binned": axes_binned}], config)
        coverage_gain = max(0.0, coverage_after - current_coverage)
        frontier_count, frontier_refs, frontier_cells = _frontier_support(
            axes_binned,
            coverage_state,
            discovery_records,
            axis_specs=axis_specs,
            axis_order=axis_order,
            graph=graph,
        )
        if coverage_gain <= 0.0 and frontier_count <= 0:
            continue

        novel_axes = sorted(axis_name for axis_name in axis_order if axes_binned[axis_name] not in visited_per_axis[axis_name])
        source_refs: list[dict[str, object]] = []
        if latest_coverage_ref is not None:
            source_refs.append(latest_coverage_ref)
        source_refs.extend(frontier_refs)

        confidence = "high" if coverage_gain > 0.0 and frontier_count > 0 else "medium"
        if coverage_gain <= 0.0 and frontier_count > 0:
            confidence = "frontier_only"

        signals: list[str] = []
        if coverage_gain > 0.0:
            signals.append("coverage_gap")
        if frontier_count > 0:
            signals.append("discovery_frontier")

        proposals.append(
            {
                "goal_id": f"goal::{cell_id}",
                "goal_type": "parameter_space_cell",
                "cell_id": cell_id,
                "axes_binned": axes_binned,
                "sort_key": (
                    -coverage_gain,
                    -frontier_count,
                    -len(novel_axes),
                    cell_id,
                ),
                "rationale": {
                    "signals": signals,
                    "confidence": confidence,
                    "coverage_gap": {
                        "current_coverage_fraction": current_coverage,
                        "projected_coverage_fraction": coverage_after,
                        "coverage_gain": coverage_gain,
                        "novel_axes": novel_axes,
                    },
                    "discovery_frontier": {
                        "frontier_neighbor_count": frontier_count,
                        "neighbor_cell_ids": frontier_cells,
                    },
                    "source_refs": source_refs,
                },
            }
        )

    proposals.sort(key=lambda item: cast(tuple[object, ...], item["sort_key"]))
    selected = proposals[:max_new_runs]
    for rank, proposal in enumerate(selected, start=1):
        rationale = cast(dict[str, object], proposal["rationale"])
        rationale["budget"] = {
            "max_new_runs": max_new_runs,
            "selected_rank": rank,
            "remaining_after_selection": max_new_runs - rank,
            "duplicate_suppression": {
                "scheduled_duplicates": duplicate_scheduled_goals,
                "covered_regions": suppressed_as_covered,
                "scheduled_regions": suppressed_as_scheduled,
            },
        }
        del proposal["sort_key"]

    saturation_reason = "coverage_saturated"
    if selected:
        saturation_reason = "budget_available"

    return {
        "status": "success" if selected else "saturated",
        "max_new_runs": max_new_runs,
        "coverage_fraction": current_coverage,
        "candidate_cells_considered": len(proposals),
        "suppressed_duplicates": {
            "scheduled_duplicates": duplicate_scheduled_goals,
            "covered_regions": suppressed_as_covered,
            "scheduled_regions": suppressed_as_scheduled,
        },
        "fallback_policy": None if selected else saturation_reason,
        "proposals": selected,
    }
