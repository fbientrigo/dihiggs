# Issue #62 Phase A inventory

Cutoff: `6bfad7662fd87750d838bf2fe0bd7ac00ee2326a`, the corrected Yukawa
installation-order fix. This inventory distinguishes historical evidence from
current physics output; it does not promote any pre-fix number.

| Artifact class | Classification | Exact replay inputs | Current consumer | Action |
|---|---|---|---|---|
| Five canonical point-v2 anchors | `ALREADY_POST_FIX` | Yes | bounded verification | Replay into this campaign’s new path |
| Lambda1 v2 Yukawa pilot | `ALREADY_POST_FIX` | Yes | verification report | Retain; no #62 replacement required |
| Lambda1 v1 golden fixture | `QUARANTINE` | Legacy case vectors only | characterization tests | Exclude from physics claims; preserve history |
| Lifetime-recovery audit CSV | `UNRESOLVABLE_PROVENANCE` | No: LFS/external row manifest absent | None | Do not infer coordinates |
| Autoresearch cτ logs | `UNRESOLVABLE_PROVENANCE` | No: LFS/external inputs absent | None; workflow frozen | Quarantine |
| Legacy scan matrices/revisions (215 members) | `UNRESOLVABLE_PROVENANCE` | No: payloads/coordinates/manifests absent | None | Retrieve primary evidence first |
| Deleted scratch/derived observables | `OBSOLETE` | No | None | Git history only |

No maintained downstream consumer currently reads an uncorrected pre-fix
observable. The bounded replay is therefore the smallest reproducible campaign;
the unavailable broad scan remains explicitly blocked rather than fabricated.
The machine-readable record is `PHASE_A_INVENTORY.json`.
