#!/usr/bin/env bash

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$SCRIPT_DIR"

PYTHON_BIN="${PYTHON_BIN:-python3}"
VENV_DIR="${VENV_DIR:-$REPO_ROOT/.venv}"
STATE_DIR="${STATE_DIR:-$REPO_ROOT/.install-state}"
DATASET_ROOT="${DATASET_ROOT:-$REPO_ROOT/datasets}"
HIGGSTOOLS_SRC="${HIGGSTOOLS_SRC:-$REPO_ROOT/higgstools}"
HIGGSTOOLS_PREFIX="${HIGGSTOOLS_PREFIX:-$HIGGSTOOLS_SRC/installation}"
HB_DATASET="${HB_DATASET:-$DATASET_ROOT/HBDataset}"
HS_DATASET="${HS_DATASET:-$DATASET_ROOT/HSDataset}"

HIGGSTOOLS_URL="https://gitlab.com/higgsbounds/higgstools.git"
REQ_FILE="$REPO_ROOT/requirements.txt"
GET_DATASETS_SCRIPT="$REPO_ROOT/get_datasets.sh"

CURRENT_STAGE=""

log() {
  printf '[install] %s\n' "$*"
}

warn() {
  printf '[install][warn] %s\n' "$*" >&2
}

die() {
  printf '[install][error] %s\n' "$*" >&2
  exit 1
}

trap 'die "stage ${CURRENT_STAGE:-unknown} failed at line ${LINENO}"' ERR

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "required command not found: $1"
}

stage_done() {
  [[ -f "$STATE_DIR/$1.done" ]]
}

mark_stage_done() {
  touch "$STATE_DIR/$1.done"
}

run_stage() {
  local stage_name="$1"
  shift

  CURRENT_STAGE="$stage_name"
  if stage_done "$stage_name"; then
    log "skipping $stage_name (already completed)"
    return 0
  fi

  log "starting $stage_name"
  "$@"
  mark_stage_done "$stage_name"
  log "completed $stage_name"
}

write_activation_script() {
  cat > "$STATE_DIR/activate.sh" <<EOF
#!/usr/bin/env bash
export HB_DATASET="$HB_DATASET"
export HS_DATASET="$HS_DATASET"
export HIGGSTOOLS_PREFIX="$HIGGSTOOLS_PREFIX"
source "$VENV_DIR/bin/activate"
EOF
  chmod +x "$STATE_DIR/activate.sh"
}

stage_python_env() {
  require_cmd "$PYTHON_BIN"
  if [[ ! -d "$VENV_DIR" ]]; then
    log "creating virtual environment at $VENV_DIR"
    "$PYTHON_BIN" -m venv "$VENV_DIR"
  fi
  "$VENV_DIR/bin/python" -m pip install --upgrade pip setuptools wheel
  [[ -x "$VENV_DIR/bin/python" ]] || die "virtual environment python not available at $VENV_DIR/bin/python"
}

stage_python_deps() {
  [[ -f "$REQ_FILE" ]] || die "missing requirements file: $REQ_FILE"
  "$VENV_DIR/bin/python" -m pip install -r "$REQ_FILE"
}

stage_build_2hdmc() {
  require_cmd make
  require_cmd g++
  log "building 2HDMC static library"
  make -C "$REPO_ROOT/2hdmc" lib
  [[ -f "$REPO_ROOT/2hdmc/lib/lib2HDMC.a" ]] || die "2HDMC library was not created"
}

stage_build_higgstools() {
  require_cmd git
  require_cmd cmake

  if [[ ! -d "$HIGGSTOOLS_SRC" ]]; then
    log "cloning HiggsTools into $HIGGSTOOLS_SRC"
    git clone "$HIGGSTOOLS_URL" "$HIGGSTOOLS_SRC"
  fi

  mkdir -p "$HIGGSTOOLS_PREFIX"
  log "configuring HiggsTools"
  cmake -S "$HIGGSTOOLS_SRC" -B "$HIGGSTOOLS_SRC/build" \
    -DCMAKE_INSTALL_PREFIX="$HIGGSTOOLS_PREFIX" \
    -DHiggsTools_BUILD_TESTING=OFF \
    -DHiggsTools_BUILD_EXAMPLES=OFF

  local jobs
  jobs="$(getconf _NPROCESSORS_ONLN 2>/dev/null || nproc 2>/dev/null || echo 1)"
  log "building HiggsTools with $jobs parallel jobs"
  cmake --build "$HIGGSTOOLS_SRC/build" --parallel "$jobs"
  cmake --install "$HIGGSTOOLS_SRC/build"
  [[ -d "$HIGGSTOOLS_PREFIX/include/Higgs" ]] || die "HiggsTools headers were not installed"
  find "$HIGGSTOOLS_PREFIX/lib" -maxdepth 1 -name 'libHiggsTools*' | grep -q . || die "HiggsTools library was not installed"
}

stage_datasets() {
  require_cmd git
  [[ -x "$GET_DATASETS_SCRIPT" ]] || chmod +x "$GET_DATASETS_SCRIPT"
  log "downloading datasets into $DATASET_ROOT"
  "$GET_DATASETS_SCRIPT" "$DATASET_ROOT"
  [[ -d "$HB_DATASET" ]] || die "HB dataset directory missing: $HB_DATASET"
  [[ -d "$HS_DATASET" ]] || die "HS dataset directory missing: $HS_DATASET"
}

stage_build_dihiggs() {
  require_cmd make
  require_cmd g++
  log "building dihiggs core project"
  make -C "$REPO_ROOT/dihiggs" clean || true
  make -C "$REPO_ROOT/dihiggs"
  [[ -f "$REPO_ROOT/dihiggs/lib/libdihiggs.a" ]] || die "dihiggs library was not created"
  [[ -x "$REPO_ROOT/dihiggs/app/PhysLam1Scan" ]] || die "PhysLam1Scan executable was not created"
}

main() {
  mkdir -p "$STATE_DIR" "$DATASET_ROOT"

  require_cmd "$PYTHON_BIN"
  require_cmd git
  require_cmd make
  require_cmd g++
  require_cmd cmake

  run_stage python_env stage_python_env
  run_stage python_deps stage_python_deps
  run_stage build_2hdmc stage_build_2hdmc
  run_stage build_higgstools stage_build_higgstools
  run_stage datasets stage_datasets
  run_stage build_dihiggs stage_build_dihiggs

  write_activation_script

  log "installation finished"
  log "to activate the environment, run: source $STATE_DIR/activate.sh"
}

main "$@"