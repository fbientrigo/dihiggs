from __future__ import annotations

# pyright: reportUnknownVariableType=false, reportUnknownMemberType=false, reportUnnecessaryIsInstance=false

from collections.abc import Mapping, Sequence
from typing import Final

METRIC_NAMES: Final[tuple[str, str, str, str]] = ("yield", "coverage", "diversity", "composite")
DEFAULT_THRESHOLDS: Final[dict[str, float]] = {
    "coverage": 0.01,
    "diversity": 0.05,
    "yield": 0.005,
    "composite": 0.01,
}
EPSILON: Final[float] = 1e-12


def _as_int(value: object, default: int = 0) -> int:
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        return int(value)
    if isinstance(value, str):
        try:
            return int(value)
        except ValueError:
            return default
    return default


def _as_float(value: object, default: float = 0.0) -> float:
    if isinstance(value, bool):
        return float(value)
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        try:
            return float(value)
        except ValueError:
            return default
    return default


def _as_mapping(value: object) -> Mapping[str, object]:
    if isinstance(value, Mapping):
        return value
    return {}


class ConvergenceDetector:
    def __init__(self, config: Mapping[str, object]):
        self.config: dict[str, object] = dict(config)
        self._convergence_cfg: Mapping[str, object] = _as_mapping(self.config.get("convergence", self.config))
        self.window_rounds: int = max(_as_int(self._convergence_cfg.get("window_rounds", 5), 5), 1)
        self.confirmation_required: int = max(
            _as_int(self._convergence_cfg.get("confirmation_required", 2), 2),
            1,
        )
        threshold_cfg = _as_mapping(self._convergence_cfg.get("thresholds"))
        self.thresholds: dict[str, float] = {
            metric: _as_float(threshold_cfg.get(metric, DEFAULT_THRESHOLDS[metric]), DEFAULT_THRESHOLDS[metric])
            for metric in DEFAULT_THRESHOLDS
        }
        configured_warmup = self._convergence_cfg.get("warmup_rounds")
        self._explicit_warmup_rounds: int | None = (
            max(_as_int(configured_warmup, 0), 0) if configured_warmup is not None else None
        )
        arms_obj = self.config.get("arms")
        arms = arms_obj if isinstance(arms_obj, list) else []
        self._configured_arm_ids: list[str] = [
            str(arm.get("id"))
            for arm in arms
            if isinstance(arm, Mapping) and arm.get("id")
        ]

    def evaluate(
        self,
        snapshot_history: Sequence[Mapping[str, object]],
        arm_pulls: Mapping[str, object],
    ) -> dict[str, object]:
        snapshots = [snapshot for snapshot in snapshot_history if isinstance(snapshot, Mapping)]
        arm_pull_counts = {str(arm_id): _as_int(count) for arm_id, count in arm_pulls.items()}
        warmup_rounds = self._required_warmup_rounds(arm_pull_counts)
        missing_arms = self._missing_arm_ids(arm_pull_counts)

        if len(snapshots) < warmup_rounds:
            return self._result(
                state="warmup_incomplete",
                reason=(
                    f"Warmup incomplete: {len(snapshots)}/{warmup_rounds} rounds observed; "
                    f"convergence requires max(2*arms, 8) rounds before applying thresholds "
                    f"coverage<={self.thresholds['coverage']}, diversity<={self.thresholds['diversity']}, "
                    f"yield<={self.thresholds['yield']}, composite<={self.thresholds['composite']}."
                ),
                warmup_rounds=warmup_rounds,
                plateau_confirmations=0,
                empty_confirmations=0,
            )

        if missing_arms:
            return self._result(
                state="warmup_incomplete",
                reason=(
                    f"Warmup arm coverage incomplete: each arm must be pulled at least once before convergence; "
                    f"missing pulls for {', '.join(missing_arms)}."
                ),
                warmup_rounds=warmup_rounds,
                plateau_confirmations=0,
                empty_confirmations=0,
            )

        plateau_confirmations = self._count_consecutive_plateau_confirmations(snapshots)
        empty_confirmations = self._count_consecutive_empty_confirmations(snapshots)

        if empty_confirmations >= self.confirmation_required:
            discoveries = _as_int(snapshots[-1].get("discoveries", 0))
            return self._result(
                state="converged_empty_exhaustion",
                reason=(
                    f"Converged by empty-round exhaustion after {empty_confirmations} consecutive confirmations: "
                    f"events_emitted==0 for the latest rounds (latest discoveries={discoveries} kept for diagnostics)."
                ),
                warmup_rounds=warmup_rounds,
                plateau_confirmations=plateau_confirmations,
                empty_confirmations=empty_confirmations,
            )

        if plateau_confirmations >= self.confirmation_required:
            deltas = self._current_plateau_deltas(snapshots)
            return self._result(
                state="converged_plateau",
                reason=(
                    f"Converged by plateau after {plateau_confirmations} consecutive confirmations over a "
                    f"{self.window_rounds}-round trailing window; {self._format_threshold_reason(deltas)}."
                ),
                warmup_rounds=warmup_rounds,
                plateau_confirmations=plateau_confirmations,
                empty_confirmations=empty_confirmations,
            )

        if empty_confirmations > 0:
            discoveries = _as_int(snapshots[-1].get("discoveries", 0))
            return self._result(
                state="plateau_warning",
                reason=(
                    f"Empty-round warning {empty_confirmations}/{self.confirmation_required}: events_emitted==0 "
                    f"is the exhaustion basis; latest discoveries={discoveries} retained for diagnostics."
                ),
                warmup_rounds=warmup_rounds,
                plateau_confirmations=plateau_confirmations,
                empty_confirmations=empty_confirmations,
            )

        if plateau_confirmations > 0:
            deltas = self._current_plateau_deltas(snapshots)
            return self._result(
                state="plateau_warning",
                reason=(
                    f"Plateau warning {plateau_confirmations}/{self.confirmation_required} over a "
                    f"{self.window_rounds}-round trailing window; {self._format_threshold_reason(deltas)}."
                ),
                warmup_rounds=warmup_rounds,
                plateau_confirmations=plateau_confirmations,
                empty_confirmations=empty_confirmations,
            )

        return self._result(
            state="running",
            reason=(
                f"Running: warmup complete and every arm has been pulled; no confirmed plateau yet under "
                f"coverage<={self.thresholds['coverage']}, diversity<={self.thresholds['diversity']}, "
                f"yield<={self.thresholds['yield']}, composite<={self.thresholds['composite']}."
            ),
            warmup_rounds=warmup_rounds,
            plateau_confirmations=0,
            empty_confirmations=0,
        )

    def _required_warmup_rounds(self, arm_pulls: Mapping[str, int]) -> int:
        inferred_arm_count = len(self._configured_arm_ids) or len(arm_pulls)
        computed_warmup = max(2 * inferred_arm_count, 8)
        if self._explicit_warmup_rounds is None:
            return computed_warmup
        return max(self._explicit_warmup_rounds, computed_warmup if self._configured_arm_ids else self._explicit_warmup_rounds)

    def _missing_arm_ids(self, arm_pulls: Mapping[str, int]) -> list[str]:
        if self._configured_arm_ids:
            required_arm_ids = self._configured_arm_ids
        else:
            required_arm_ids = list(arm_pulls)
        return [arm_id for arm_id in required_arm_ids if arm_pulls.get(arm_id, 0) < 1]

    def _count_consecutive_empty_confirmations(self, snapshots: Sequence[Mapping[str, object]]) -> int:
        confirmations = 0
        for snapshot in reversed(snapshots):
            if _as_int(snapshot.get("events_emitted", 0)) == 0:
                confirmations += 1
                continue
            break
        return confirmations

    def _count_consecutive_plateau_confirmations(self, snapshots: Sequence[Mapping[str, object]]) -> int:
        confirmations = 0
        for end_index in range(len(snapshots) - 1, -1, -1):
            if not self._is_plateau_confirmation(snapshots, end_index):
                break
            confirmations += 1
        return confirmations

    def _is_plateau_confirmation(
        self,
        snapshots: Sequence[Mapping[str, object]],
        end_index: int,
    ) -> bool:
        return self._plateau_deltas_at(snapshots, end_index) is not None

    def _current_plateau_deltas(self, snapshots: Sequence[Mapping[str, object]]) -> dict[str, float] | None:
        return self._plateau_deltas_at(snapshots, len(snapshots) - 1)

    def _plateau_deltas_at(
        self,
        snapshots: Sequence[Mapping[str, object]],
        end_index: int,
    ) -> dict[str, float] | None:
        snapshot = snapshots[end_index]
        rolling_deltas = _as_mapping(snapshot.get("rolling_deltas"))
        deltas: dict[str, float] = {}
        for metric in METRIC_NAMES:
            delta_key = f"{metric}_delta"
            if delta_key in rolling_deltas:
                delta_value = rolling_deltas.get(delta_key)
                if delta_value is None:
                    return None
                deltas[metric] = _as_float(delta_value)
            else:
                prior_index = end_index - self.window_rounds
                if prior_index < 0:
                    return None
                deltas[metric] = self._metric_value(snapshots[end_index], metric) - self._metric_value(
                    snapshots[prior_index],
                    metric,
                )

        if all(deltas[metric] <= self.thresholds[metric] + EPSILON for metric in self.thresholds):
            return deltas
        return None

    @staticmethod
    def _metric_value(snapshot: Mapping[str, object], metric: str) -> float:
        rolling_metrics = _as_mapping(snapshot.get("rolling_metrics"))
        if metric in rolling_metrics:
            return _as_float(rolling_metrics.get(metric))
        return _as_float(snapshot.get(metric))

    def _format_threshold_reason(self, deltas: dict[str, float] | None) -> str:
        if deltas is None:
            return (
                f"thresholds are coverage<={self.thresholds['coverage']}, diversity<={self.thresholds['diversity']}, "
                f"yield<={self.thresholds['yield']}, composite<={self.thresholds['composite']}"
            )
        return (
            f"deltas coverage={deltas['coverage']:.4f}<={self.thresholds['coverage']}, "
            f"diversity={deltas['diversity']:.4f}<={self.thresholds['diversity']}, "
            f"yield={deltas['yield']:.4f}<={self.thresholds['yield']}, "
            f"composite={deltas['composite']:.4f}<={self.thresholds['composite']} "
            f"(thresholds: coverage<={self.thresholds['coverage']}, diversity<={self.thresholds['diversity']}, "
            f"yield<={self.thresholds['yield']}, composite<={self.thresholds['composite']})"
        )

    def _result(
        self,
        *,
        state: str,
        reason: str,
        warmup_rounds: int,
        plateau_confirmations: int,
        empty_confirmations: int,
    ) -> dict[str, object]:
        return {
            "state": state,
            "reason": reason,
            "warmup_rounds": warmup_rounds,
            "window_rounds": self.window_rounds,
            "confirmation_required": self.confirmation_required,
            "plateau_confirmations": plateau_confirmations,
            "empty_confirmations": empty_confirmations,
            "thresholds": dict(self.thresholds),
        }
