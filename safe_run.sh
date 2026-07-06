#!/usr/bin/env bash
# safe_run.sh — polite, resource-capped launcher ("execution server" front-end).
# =============================================================================
# Wrap ANY long runner so it can never take the host down and always leaves the
# machine interactively reachable. It does four things automatically:
#
#   1. CPU headroom   — caps worker threads at (nproc - RESERVE) via every common
#                       env knob (OpenMP / OpenBLAS / MKL / numexpr / numba), so
#                       at least RESERVE cores stay free for the OS and sshd.
#   2. Hard limits    — runs the command inside a transient systemd --user scope
#                       with CPUQuota + MemoryHigh/MemoryMax/MemorySwapMax, so a
#                       runaway is throttled/killed *inside its own cgroup*
#                       instead of freezing the box. Degrades gracefully if a
#                       controller isn't delegated.
#   3. Polite         — nice + idle ionice, so interactive work always preempts.
# Everything is tee'd to a timestamped log under DIHIGGS_RUN_LOGS.
#
# WHY THIS EXISTS: jun 2026 a `--threads 12` gen_fixings run on a 12-core host
# pinned every core, starved sshd (port 22 timed out), stopped logging mid-run,
# and the box had to be power-cycled. tmux did NOT help: tmux survives an SSH
# *disconnect*, but the tmux server lives on the frozen host, so a hard CPU/IO
# starvation freeze takes it down too. The fix is headroom + cgroup caps so the
# host never starves in the first place.
#
# USAGE
#   ./safe_run.sh [options] -- <command> [args...]
#
# OPTIONS
#   --reserve N        cores to keep free for the OS        (default 2)
#   --threads N        override worker thread count         (default nproc-reserve)
#   --nice N           niceness 0..19                       (default 19)
#   --mem-high-pct P   soft memory throttle, % of RAM       (default 70)
#   --mem-max-pct P    hard memory cap, % of RAM            (default 85)
#   --swap-max SIZE    per-job swap cap (e.g. 2G, 0)        (default 2G)
#   --name NAME        label for the scope unit + log file
#   --log-dir DIR      where to write logs                 (default $DIHIGGS_RUN_LOGS
#                                                            or ~/dihiggs/data/_runlogs)
#   --no-cgroup        skip systemd scope; nice/ionice + thread cap only
#   -h | --help        this help
#
# EXAMPLES
#   # A gen_fixings scan, hard-capped:
#   ./safe_run.sh --name lam1eq1_var1000 \
#       -- python -m dihiggs.app.orchestrator --engine gen_fixings \
#            --exec dihiggs/app/GenScanWithFixings --campaign lam1eq1_var1000 \
#            --bronze-csv DATA/.../bronze_all.csv --calibration-n 1000 \
#            --tanbeta 10 --outdir data/lam1eq1_var1000/lake
#
#   # Wrap any other runner politely:
#   ./safe_run.sh --name nightscan -- bash run_quarantine_night.sh
# =============================================================================
set -uo pipefail

# ---- defaults ---------------------------------------------------------------
RESERVE_CORES=2
THREADS=""
NICE_LEVEL=19
MEM_HIGH_PCT=70
MEM_MAX_PCT=85
SWAP_MAX="2G"
USE_CGROUP=1
NAME=""
LOG_DIR="${DIHIGGS_RUN_LOGS:-$HOME/dihiggs/data/_runlogs}"

usage() { sed -n '2,55p' "$0" | sed 's/^# \{0,1\}//'; }

# ---- parse options up to `--` ----------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --reserve)      RESERVE_CORES="$2"; shift 2;;
    --threads)      THREADS="$2"; shift 2;;
    --nice)         NICE_LEVEL="$2"; shift 2;;
    --mem-high-pct) MEM_HIGH_PCT="$2"; shift 2;;
    --mem-max-pct)  MEM_MAX_PCT="$2"; shift 2;;
    --swap-max)     SWAP_MAX="$2"; shift 2;;
    --name)         NAME="$2"; shift 2;;
    --log-dir)      LOG_DIR="$2"; shift 2;;
    --no-cgroup)    USE_CGROUP=0; shift;;
    -h|--help)      usage; exit 0;;
    --)             shift; break;;
    *) echo "safe_run: unknown option '$1'" >&2; usage; exit 2;;
  esac
done

if [[ $# -eq 0 ]]; then
  echo "safe_run: no command given (did you forget '--'?)" >&2
  usage; exit 2
fi
CMD=( "$@" )

# ---- compute resources ------------------------------------------------------
NPROC="$(nproc)"
if [[ -z "$THREADS" ]]; then
  THREADS=$(( NPROC - RESERVE_CORES ))
  (( THREADS < 1 )) && THREADS=1
fi
CPU_QUOTA=$(( THREADS * 100 ))   # percent; 10 cores -> 1000%
MEM_TOTAL_KB="$(awk '/MemTotal/{print $2}' /proc/meminfo)"
MEM_TOTAL_MB=$(( MEM_TOTAL_KB / 1024 ))
MEM_HIGH_MB=$(( MEM_TOTAL_MB * MEM_HIGH_PCT / 100 ))
MEM_MAX_MB=$(( MEM_TOTAL_MB * MEM_MAX_PCT / 100 ))

# ---- thread caps (C/C++ OpenMP and Python BLAS/numexpr/numba) ---------------
export OMP_NUM_THREADS="$THREADS"
export OMP_THREAD_LIMIT="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"
export NUMEXPR_NUM_THREADS="$THREADS"
export VECLIB_MAXIMUM_THREADS="$THREADS"
export NUMBA_NUM_THREADS="$THREADS"

# ---- logging ----------------------------------------------------------------
mkdir -p "$LOG_DIR"
ts="$(date +%Y%m%d_%H%M%S)"
slug="${NAME:-run}"; slug="${slug//[^A-Za-z0-9_.-]/_}"
LOG_FILE="$LOG_DIR/${ts}_${slug}.log"
SCOPE_NAME="dihiggs-${slug}-${ts}-$$"

log() { echo "[safe_run $(date '+%F %T')] $*" | tee -a "$LOG_FILE"; }

# ---- decide on cgroup hardening + which controllers are usable --------------
declare -a SCOPE_PROPS=()
HARDENED=0
build_scope_props() {
  local cfile="/sys/fs/cgroup/user.slice/user-$(id -u).slice/user@$(id -u).service/cgroup.controllers"
  local ctrls=""
  [[ -r "$cfile" ]] && ctrls=" $(cat "$cfile") "
  if [[ "$ctrls" == *" cpu "* ]]; then
    SCOPE_PROPS+=( -p "CPUQuota=${CPU_QUOTA}%" )
  else
    log "note: 'cpu' controller not delegated to the user manager; CPU headroom"
    log "      is still enforced by the thread cap (OMP_NUM_THREADS=$THREADS)."
  fi
  if [[ "$ctrls" == *" memory "* ]]; then
    SCOPE_PROPS+=( -p "MemoryHigh=${MEM_HIGH_MB}M" \
                   -p "MemoryMax=${MEM_MAX_MB}M" \
                   -p "MemorySwapMax=${SWAP_MAX}" )
  else
    log "note: 'memory' controller not delegated; memory hard-cap unavailable."
  fi
}

CAN_SCOPE=0
if [[ "$USE_CGROUP" -eq 1 ]] \
   && command -v systemd-run >/dev/null 2>&1 \
   && [[ -n "${XDG_RUNTIME_DIR:-}" ]] \
   && systemctl --user is-active --quiet default.target 2>/dev/null; then
  CAN_SCOPE=1
  build_scope_props
  [[ ${#SCOPE_PROPS[@]} -gt 0 ]] && HARDENED=1
fi

POLITE=( nice -n "$NICE_LEVEL" ionice -c3 )

# Run a command under the chosen limits. Echoes nothing; caller pipes to tee.
run_capped() {
  if [[ "$CAN_SCOPE" -eq 1 ]]; then
    systemd-run --user --scope --collect --quiet \
      --unit="${SCOPE_NAME}-$RANDOM" \
      "${SCOPE_PROPS[@]}" \
      "${POLITE[@]}" "$@"
  else
    "${POLITE[@]}" "$@"
  fi
}

# ---- banner -----------------------------------------------------------------
log "host=$(hostname)  nproc=$NPROC  reserve=$RESERVE_CORES  threads=$THREADS"
log "cpu_quota=${CPU_QUOTA}%  mem_high=${MEM_HIGH_MB}M  mem_max=${MEM_MAX_MB}M  swap_max=${SWAP_MAX}  nice=$NICE_LEVEL"
if [[ "$CAN_SCOPE" -eq 1 && "$HARDENED" -eq 1 ]]; then
  log "hardening: systemd --user scope (cgroup limits active)"
elif [[ "$CAN_SCOPE" -eq 1 ]]; then
  log "hardening: systemd scope present but no controllers delegated; thread cap + nice/ionice only"
else
  log "hardening: nice/ionice + thread cap only (systemd --user scope unavailable)"
fi
log "command: ${CMD[*]}"
log "log: $LOG_FILE"

# ---- run main command -------------------------------------------------------
t0=$(date +%s)
run_capped "${CMD[@]}" 2>&1 | tee -a "$LOG_FILE"
rc=${PIPESTATUS[0]}
dt=$(( $(date +%s) - t0 ))
log "main command finished rc=$rc in ${dt}s"

log "done rc=$rc"
exit "$rc"
