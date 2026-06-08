"""
parsing.py — Parse C++ stdout for scan statistics.

Handles two output tokens that the C++ executables may print:

  Total Attempts: <int>
  TRIPLE_OK_POINTS <int_or_str>

These are optional: if the binary does not emit them, the parser returns
None.  Both PhysScanWithFixings and Phys_M2BoundaryScan share the same
output format for these two tokens.
"""

from __future__ import annotations

from typing import Optional, Tuple


def parse_cpp_stats(stdout: str) -> Tuple[Optional[int], Optional[str]]:
    """
    Parse key statistics from C++ stdout.

    Parameters
    ----------
    stdout:
        Full stdout string captured from the C++ subprocess.

    Returns
    -------
    total_attempts:
        Integer value following ``"Total Attempts:"`` if found; else ``None``.
    triple_ok_points:
        Last whitespace-delimited token on the ``"TRIPLE_OK_POINTS"`` line
        if found; else ``None``.  Returned as a string to preserve any
        suffix the binary might append (e.g. ``"42"`` or ``"42/100"``).

    Notes
    -----
    * Lines are matched by substring containment, not strict regex, to
      tolerate minor formatting variation between binary versions.
    * Only the last occurrence of each token is kept (matching legacy
      behaviour in orchestrate_scans.py L649–678).
    """
    total_attempts: Optional[int] = None
    triple_ok_points: Optional[str] = None

    for line in stdout.splitlines():
        stripped = line.strip()

        if "Total Attempts:" in stripped:
            try:
                total_attempts = int(stripped.split(":")[-1].strip())
            except (ValueError, IndexError):
                pass

        if "TRIPLE_OK_POINTS" in stripped:
            parts = stripped.split()
            if parts:
                triple_ok_points = parts[-1]

    return total_attempts, triple_ok_points
