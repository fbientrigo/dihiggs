# Replay of the Q,S continuation pilot under active m_h = 125.20 GeV

The Q,S continuation pilot sequence (150 → 175 → 200 GeV) has been replayed
under the active project convention **`m_h = 125.20 GeV`**, with full report
recorded in [`continuation_report.json`](continuation_report.json).

Historical evidence produced under the superseded `125.13 GeV` convention is
preserved honestly in
[`continuation_report_mh12513_historical.json`](continuation_report_mh12513_historical.json).

## Replay Summary (m_h = 125.20 GeV)

| Step | Target m_H2 | m_A = m_Hp (GeV) | M2 (GeV^2) | point_id | theory_ok_v1 | rejection_stage | c*tau (mm) | g_hH2H2 (GeV) | BR(bb) |
|---|---:|---:|---:|---|:-:|:-:|---:|---:|---:|
| Gate 2 (Anchor) | 150.0 | 450.000 | 22499.9999995 | `point_98c841e915d3605a` | 1 | accepted | 4.32621 | 63.6626 | 0.75674 |
| Gate 3 | 175.0 | 458.939 | 30624.9999995 | `point_bac9d1637f23ed5e` | 1 | accepted | 3.62750 | 63.6626 | 0.71921 |
| Gate 4 | 200.0 | 469.042 | 39999.9999995 | `point_e5ac522ac26a920b` | 1 | accepted | 3.05775 | 63.6626 | 0.67546 |

**Verdict: `DIRECT_QS_CONTINUATION_VALID_TO_200_GEV`** (VALID_TO_200).
All three points evaluate as theory-valid on their first attempt without requiring local search or rescue.
