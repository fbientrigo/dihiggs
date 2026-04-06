"""Stateful alert engine with deduplication and resolution tracking.

Schema version: 1

Alert types:
- preflight_blocked: Preflight check failed (CRITICAL)
- dataset_missing: HB_DATASET or HS_DATASET env vars not set (WARNING)
- round_failed: Round completed with nonzero_exit subprocess status (WARNING)
- round_crashed: Round subprocess crashed or killed (CRITICAL)
- round_timed_out: Round subprocess exceeded timeout (WARNING)
- round_stalled: Multiple consecutive empty rounds without progress (WARNING)
- plateau_warning: Metrics show plateau pattern (INFO)
- converged: Campaign reached convergence (INFO)
- telemetry_corruption: Event log or telemetry corruption detected (CRITICAL)

Severity levels: critical, warning, info

Design:
- Stateful dedupe with configurable cooldown window (default: 900 seconds)
- Tracks first_seen, last_seen, resolved_at timestamps per alert
- Emits to alerts.jsonl and stderr
- SQLite-ready persistence shape for future integration
"""
from __future__ import annotations

import hashlib
import json
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Literal

AlertSeverity = Literal["critical", "warning", "info"]
AlertType = Literal[
    "preflight_blocked",
    "dataset_missing",
    "round_failed",
    "round_crashed",
    "round_timed_out",
    "round_stalled",
    "plateau_warning",
    "converged",
    "telemetry_corruption",
]

# Severity mapping for each alert type
ALERT_SEVERITY_MAP: dict[AlertType, AlertSeverity] = {
    "preflight_blocked": "critical",
    "dataset_missing": "warning",
    "round_failed": "warning",
    "round_crashed": "critical",
    "round_timed_out": "warning",
    "round_stalled": "warning",
    "plateau_warning": "info",
    "converged": "info",
    "telemetry_corruption": "critical",
}


@dataclass
class Alert:
    """Alert record with deduplication and resolution tracking."""

    alert_type: AlertType
    severity: AlertSeverity
    message: str
    context: dict[str, object]
    first_seen: str
    last_seen: str
    resolved_at: str | None = None
    fingerprint: str = field(default="", init=False)

    def __post_init__(self) -> None:
        """Generate stable fingerprint for deduplication."""
        # Fingerprint based on alert_type + sorted context keys/values
        context_str = json.dumps(self.context, sort_keys=True, separators=(",", ":"))
        fingerprint_input = f"{self.alert_type}|{context_str}"
        self.fingerprint = hashlib.sha256(fingerprint_input.encode("utf-8")).hexdigest()[:16]


@dataclass
class AlertEngine:
    """Stateful alert engine with deduplication and cooldown."""

    config: dict[str, object]
    active_alerts: dict[str, Alert] = field(default_factory=dict)
    resolved_alerts: dict[str, Alert] = field(default_factory=dict)
    dedupe_window_sec: int = field(init=False, default=900)

    def __post_init__(self) -> None:
        """Extract configuration."""
        dedupe_sec = self.config.get("dedupe_window_sec", 900)
        if isinstance(dedupe_sec, int):
            self.dedupe_window_sec = dedupe_sec
        elif isinstance(dedupe_sec, float):
            self.dedupe_window_sec = int(dedupe_sec)
        elif isinstance(dedupe_sec, str):
            self.dedupe_window_sec = int(dedupe_sec)
        else:
            self.dedupe_window_sec = 900
    def _utc_now(self) -> str:
        """Get current UTC timestamp in ISO 8601 format."""
        return datetime.now(timezone.utc).isoformat()

    def _is_within_cooldown(self, alert: Alert, now_str: str) -> bool:
        """Check if alert is within cooldown window."""
        try:
            last_seen_dt = datetime.fromisoformat(alert.last_seen.replace("Z", "+00:00"))
            now_dt = datetime.fromisoformat(now_str.replace("Z", "+00:00"))
            elapsed = (now_dt - last_seen_dt).total_seconds()
            return elapsed < self.dedupe_window_sec
        except (ValueError, AttributeError):
            return False

    def _create_alert(
        self,
        alert_type: AlertType,
        message: str,
        context: dict[str, object] | None = None,
    ) -> Alert:
        """Create a new alert with timestamps."""
        now = self._utc_now()
        severity = ALERT_SEVERITY_MAP[alert_type]
        return Alert(
            alert_type=alert_type,
            severity=severity,
            message=message,
            context=context or {},
            first_seen=now,
            last_seen=now,
        )

    def _check_preflight_blocked(self, state: dict[str, object]) -> Alert | None:
        """Check for preflight failures."""
        preflight_overall = state.get("preflight_overall")
        if preflight_overall == "fail":
            return self._create_alert(
                "preflight_blocked",
                "Preflight check failed - campaign execution blocked",
                {"preflight_overall": preflight_overall},
            )
        return None

    def _check_dataset_missing(self, state: dict[str, object]) -> Alert | None:
        """Check for missing dataset environment variables."""
        preflight_checks = state.get("preflight_checks", [])
        if isinstance(preflight_checks, list):
            for check in preflight_checks:
                if isinstance(check, dict):
                    if check.get("check") == "datasets" and check.get("status") == "warn":
                        return self._create_alert(
                            "dataset_missing",
                            f"Dataset environment variables not properly configured: {check.get('message', '')}",
                            {"check": check},
                        )
        return None

    def _check_round_failures(self, state: dict[str, object]) -> list[Alert]:
        """Check for round execution failures."""
        alerts = []
        subprocess_status = state.get("subprocess_status")
        round_index = state.get("round_index", -1)

        if subprocess_status == "nonzero_exit":
            alerts.append(
                self._create_alert(
                    "round_failed",
                    f"Round {round_index} completed with non-zero exit code",
                    {"round_index": round_index, "subprocess_status": subprocess_status},
                )
            )
        elif subprocess_status == "timeout":
            alerts.append(
                self._create_alert(
                    "round_timed_out",
                    f"Round {round_index} exceeded timeout limit",
                    {"round_index": round_index, "subprocess_status": subprocess_status},
                )
            )
        elif subprocess_status in ("crashed", "killed", "error"):
            alerts.append(
                self._create_alert(
                    "round_crashed",
                    f"Round {round_index} crashed or was killed",
                    {"round_index": round_index, "subprocess_status": subprocess_status},
                )
            )

        return alerts

    def _check_round_stalled(self, state: dict[str, object]) -> Alert | None:
        """Check for consecutive empty rounds indicating stall."""
        consecutive_empty = state.get("consecutive_empty_rounds", 0)
        max_empty_rounds = state.get("max_empty_rounds", 3)

        # Type guard: ensure both values are integers
        if isinstance(consecutive_empty, int) and isinstance(max_empty_rounds, int):
            if consecutive_empty >= max_empty_rounds - 1 and consecutive_empty > 0:
                return self._create_alert(
                    "round_stalled",
                    f"Campaign stalled with {consecutive_empty} consecutive empty rounds",
                    {"consecutive_empty_rounds": consecutive_empty, "max_empty_rounds": max_empty_rounds},
                )
        return None

    def _check_plateau_warning(self, state: dict[str, object]) -> Alert | None:
        """Check for plateau warning signal."""
        if state.get("plateau_warning") is True:
            metrics = state.get("rolling_metrics", {})
            return self._create_alert(
                "plateau_warning",
                "Campaign metrics showing plateau pattern",
                {"rolling_metrics": metrics, "plateau_warning": True},
            )
        return None

    def _check_converged(self, state: dict[str, object]) -> Alert | None:
        """Check for convergence signal."""
        convergence_state = state.get("convergence_state", "")
        if isinstance(convergence_state, str) and convergence_state.startswith("converged"):
            return self._create_alert(
                "converged",
                f"Campaign converged: {convergence_state}",
                {"convergence_state": convergence_state},
            )
        return None

    def _check_telemetry_corruption(self, state: dict[str, object]) -> Alert | None:
        """Check for telemetry corruption signals."""
        corruption_detected = state.get("telemetry_corruption", False)
        if corruption_detected:
            return self._create_alert(
                "telemetry_corruption",
                "Telemetry corruption detected in event log or database",
                {"telemetry_corruption": True},
            )
        return None

    def evaluate(self, state: dict[str, object]) -> dict[str, object]:
        """Evaluate current campaign state and produce alerts.

        Args:
            state: Campaign state snapshot with fields:
                - campaign_state: RUNNING, CONVERGED, FAILED, etc.
                - preflight_overall: pass/warn/fail
                - preflight_checks: list of check results
                - subprocess_status: success/timeout/nonzero_exit/crashed
                - round_index: current round number
                - consecutive_empty_rounds: count of consecutive empty rounds
                - max_empty_rounds: configured threshold
                - plateau_warning: boolean plateau signal
                - convergence_state: convergence detector state
                - telemetry_corruption: boolean corruption signal
                - active_alerts: list of currently active alerts (for state restoration)
                - rolling_metrics: current metric values

        Returns:
            dict with:
                - new_alerts: list of newly triggered alerts
                - active_alerts: list of all currently active alerts
                - resolved_alerts: list of alerts resolved in this evaluation
        """
        now = self._utc_now()
        new_alerts = []
        resolved_in_this_eval = []

        # Restore active alerts from state if provided
        if "active_alerts" in state and not self.active_alerts:
            active_list = state.get("active_alerts", [])
            if isinstance(active_list, list):
                for alert_dict in active_list:
                    if isinstance(alert_dict, dict):
                        alert = Alert(
                            alert_type=alert_dict["alert_type"],
                            severity=alert_dict["severity"],
                            message=alert_dict["message"],
                            context=alert_dict.get("context", {}),
                            first_seen=alert_dict["first_seen"],
                            last_seen=alert_dict["last_seen"],
                            resolved_at=alert_dict.get("resolved_at"),
                        )
                        # Use fingerprint from restored alert dict if available
                        fp = alert_dict.get("fingerprint", alert.fingerprint)
                        # Update last_seen to now for proper cooldown checking on restore
                        alert.last_seen = now
                        self.active_alerts[fp] = alert

        # Collect all potential alerts
        potential_alerts: list[Alert] = []

        # Preflight checks
        alert = self._check_preflight_blocked(state)
        if alert:
            potential_alerts.append(alert)

        alert = self._check_dataset_missing(state)
        if alert:
            potential_alerts.append(alert)

        # Round failures
        potential_alerts.extend(self._check_round_failures(state))

        # Stall detection
        alert = self._check_round_stalled(state)
        if alert:
            potential_alerts.append(alert)

        # Plateau and convergence
        alert = self._check_plateau_warning(state)
        if alert:
            potential_alerts.append(alert)

        alert = self._check_converged(state)
        if alert:
            potential_alerts.append(alert)

        # Telemetry corruption
        alert = self._check_telemetry_corruption(state)
        if alert:
            potential_alerts.append(alert)

        # Process potential alerts with deduplication
        for alert in potential_alerts:
            existing = self.active_alerts.get(alert.fingerprint)

            if existing:
                # Alert already active - check cooldown
                if self._is_within_cooldown(existing, now):
                    # Still in cooldown, update last_seen but don't emit new alert
                    existing.last_seen = now
                else:
                    # Cooldown expired, re-emit alert
                    existing.last_seen = now
                    new_alerts.append(alert)
            else:
                # New alert
                self.active_alerts[alert.fingerprint] = alert
                new_alerts.append(alert)

        # Check for resolved alerts (active alerts no longer in potential set)
        potential_fingerprints = {a.fingerprint for a in potential_alerts}
        for fingerprint, alert in list(self.active_alerts.items()):
            if fingerprint not in potential_fingerprints:
                # Alert condition cleared - mark as resolved
                alert.resolved_at = now
                self.resolved_alerts[fingerprint] = alert
                resolved_in_this_eval.append(alert)
                del self.active_alerts[fingerprint]

        # Convert to serializable format
        return {
            "new_alerts": [
                {
                    "alert_type": a.alert_type,
                    "severity": a.severity,
                    "message": a.message,
                    "context": a.context,
                    "first_seen": a.first_seen,
                    "last_seen": a.last_seen,
                    "fingerprint": a.fingerprint,
                }
                for a in new_alerts
            ],
            "active_alerts": [
                {
                    "alert_type": a.alert_type,
                    "severity": a.severity,
                    "message": a.message,
                    "context": a.context,
                    "first_seen": a.first_seen,
                    "last_seen": a.last_seen,
                    "fingerprint": a.fingerprint,
                }
                for a in self.active_alerts.values()
            ],
            "resolved_alerts": [
                {
                    "alert_type": a.alert_type,
                    "severity": a.severity,
                    "message": a.message,
                    "context": a.context,
                    "first_seen": a.first_seen,
                    "last_seen": a.last_seen,
                    "resolved_at": a.resolved_at,
                    "fingerprint": a.fingerprint,
                }
                for a in resolved_in_this_eval
            ],
        }

    def emit_alerts(self, alerts: list[dict[str, object]], outdir: str | Path) -> None:
        """Emit new alerts to alerts.jsonl and stderr.

        Args:
            alerts: List of alert dictionaries to emit
            outdir: Campaign output directory for alerts.jsonl
        """
        if not alerts:
            return

        outpath = Path(outdir)
        outpath.mkdir(parents=True, exist_ok=True)
        alerts_file = outpath / "alerts.jsonl"

        for alert in alerts:
            # Emit to stderr
            severity_raw = alert.get("severity", "info")
            severity_str = str(severity_raw).upper() if severity_raw else "INFO"
            alert_type = alert.get("alert_type", "unknown")
            message = alert.get("message", "")
            print(f"[{severity_str}] {alert_type}: {message}", file=sys.stderr)

            # Append to alerts.jsonl
            with open(alerts_file, "a", encoding="utf-8") as f:
                f.write(json.dumps(alert, separators=(",", ":")) + "\n")
