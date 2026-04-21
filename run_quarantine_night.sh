#!/usr/bin/env bash
set -u

TOTAL_HOURS=5
TOTAL_SECONDS=$((TOTAL_HOURS * 3600))

THREADS=15
TIMEOUT_PER_CALL=21600
RETRIES=2
N_SUBSETS=8
MAX_RUNS=1

STATUS_BASE="/mnt/c/Users/Asus/cern_db/dihiggs_lake/recomputed"
MASTER_LOG="${STATUS_BASE}/night_explore_5h.log"

START_EPOCH=$(date +%s)
END_EPOCH=$((START_EPOCH + TOTAL_SECONDS))

echo "[$(date '+%F %T')] START 5h exploration window" | tee -a "$MASTER_LOG"

mapfile -t CAMPAIGNS < <(
  python scripts/run_quarantine.py --list-campaigns | awk '{print $2}'
)

for CAMPAIGN in "${CAMPAIGNS[@]}"; do
  NOW=$(date +%s)
  if (( NOW >= END_EPOCH )); then
    echo "[$(date '+%F %T')] TIME LIMIT REACHED before campaign ${CAMPAIGN}" | tee -a "$MASTER_LOG"
    break
  fi

  echo "[$(date '+%F %T')] CAMPAIGN ${CAMPAIGN}" | tee -a "$MASTER_LOG"

  for ((k=0; k<N_SUBSETS; k++)); do
    NOW=$(date +%s)
    if (( NOW >= END_EPOCH )); then
      echo "[$(date '+%F %T')] TIME LIMIT REACHED during campaign ${CAMPAIGN}" | tee -a "$MASTER_LOG"
      break 2
    fi

    SUBSET="${k}/${N_SUBSETS}"
    LOG_FILE="${STATUS_BASE}/explore_${CAMPAIGN#campaign=}_subset${k}of${N_SUBSETS}.log"

    attempt=0
    success=0

    while (( attempt <= RETRIES )); do
      NOW=$(date +%s)
      if (( NOW >= END_EPOCH )); then
        echo "[$(date '+%F %T')] TIME LIMIT REACHED before subset ${SUBSET} of ${CAMPAIGN}" | tee -a "$MASTER_LOG"
        break 3
      fi

      REMAINING=$((END_EPOCH - NOW))
      CALL_TIMEOUT=$TIMEOUT_PER_CALL
      if (( REMAINING < CALL_TIMEOUT )); then
        CALL_TIMEOUT=$REMAINING
      fi
      if (( CALL_TIMEOUT <= 0 )); then
        echo "[$(date '+%F %T')] No time left for ${CAMPAIGN} subset ${SUBSET}" | tee -a "$MASTER_LOG"
        break 3
      fi

      echo "[$(date '+%F %T')] START campaign=${CAMPAIGN} subset=${SUBSET} attempt=${attempt} remaining_sec=${REMAINING}" | tee -a "$MASTER_LOG" | tee -a "$LOG_FILE"

      python scripts/run_quarantine.py \
        --campaign "$CAMPAIGN" \
        --subset "$SUBSET" \
        --max-runs "$MAX_RUNS" \
        --yukawas-type 1 \
        --omp-num-threads "$THREADS" \
        --process-timeout-sec "$CALL_TIMEOUT" \
        --status-csv "${STATUS_BASE}/_status_${CAMPAIGN#campaign=}_subset${k}of${N_SUBSETS}.csv" \
        --status-json "${STATUS_BASE}/_status_${CAMPAIGN#campaign=}_subset${k}of${N_SUBSETS}.json" \
        2>&1 | tee -a "$MASTER_LOG" | tee -a "$LOG_FILE"

      rc=${PIPESTATUS[0]}

      if [[ "$rc" -eq 0 ]]; then
        echo "[$(date '+%F %T')] OK campaign=${CAMPAIGN} subset=${SUBSET}" | tee -a "$MASTER_LOG" | tee -a "$LOG_FILE"
        success=1
        break
      fi

      echo "[$(date '+%F %T')] FAIL rc=${rc} campaign=${CAMPAIGN} subset=${SUBSET}" | tee -a "$MASTER_LOG" | tee -a "$LOG_FILE"
      attempt=$((attempt + 1))

      if (( attempt <= RETRIES )); then
        echo "[$(date '+%F %T')] RETRY in 30s..." | tee -a "$MASTER_LOG" | tee -a "$LOG_FILE"
        sleep 30
      fi
    done

    if [[ "$success" -eq 0 ]]; then
      echo "[$(date '+%F %T')] SKIP_TO_NEXT_SUBSET campaign=${CAMPAIGN} subset=${SUBSET}" | tee -a "$MASTER_LOG" | tee -a "$LOG_FILE"
    fi
  done
done

echo "[$(date '+%F %T')] FINISHED 5h exploration window" | tee -a "$MASTER_LOG"