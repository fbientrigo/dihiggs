#!/usr/bin/env bash

set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
REPO_ROOT="$SCRIPT_DIR"

export PATH="$HOME/.local/bin:$PATH"

PYTHON_BIN="${PYTHON_BIN:-python3}"
PYTHON_SPEC="${PYTHON_SPEC:-3.11}"

VENV_DIR="${VENV_DIR:-$REPO_ROOT/.venv}"
STATE_DIR="${STATE_DIR:-$REPO_ROOT/.install-state}"
DATASET_ROOT="${DATASET_ROOT:-$REPO_ROOT/datasets}"
HIGGSTOOLS_SRC="${HIGGSTOOLS_SRC:-$REPO_ROOT/higgstools}"
HIGGSTOOLS_PREFIX="${HIGGSTOOLS_PREFIX:-$HIGGSTOOLS_SRC/installation}"
HB_DATASET="${HB_DATASET:-$DATASET_ROOT/HBDataset}"
HS_DATASET="${HS_DATASET:-$DATASET_ROOT/HSDataset}"

HIGGSTOOLS_URL="${HIGGSTOOLS_URL:-https://gitlab.com/higgsbounds/higgstools.git}"
REQ_FILE="${REQ_FILE:-$REPO_ROOT/requirements.txt}"
GET_DATASETS_SCRIPT="${GET_DATASETS_SCRIPT:-$REPO_ROOT/get_datasets.sh}"
PROJECT_CONFIG_FILE="${PROJECT_CONFIG_FILE:-$REPO_ROOT/project_config.json}"
DATA_LAKE_DIR="${DATA_LAKE_DIR:-$HOME/dihiggs_lake}"

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
  mkdir -p "$STATE_DIR"
  touch "$STATE_DIR/$1.done"
}

clear_stage_done() {
  mkdir -p "$STATE_DIR"
  local stage_name
  for stage_name in "$@"; do
    rm -f "$STATE_DIR/$stage_name.done"
  done
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

num_jobs() {
  getconf _NPROCESSORS_ONLN 2>/dev/null || nproc 2>/dev/null || echo 1
}

python_mm() {
  local py="$1"
  "$py" - <<'PY'
import sys
print(f"{sys.version_info.major}.{sys.version_info.minor}")
PY
}

python_is_supported_for_repo() {
  local py="$1"
  "$py" - <<'PY'
import sys
raise SystemExit(0 if (3, 8) <= sys.version_info[:2] < (3, 12) else 1)
PY
}

python_is_preferred() {
  local py="$1"
  "$py" - <<'PY'
import sys
raise SystemExit(0 if sys.version_info[:2] == (3, 11) else 1)
PY
}

bootstrap_uv() {
  if command -v uv >/dev/null 2>&1; then
    return 0
  fi

  log "uv not found; bootstrapping uv"

  if command -v curl >/dev/null 2>&1; then
    curl -LsSf https://astral.sh/uv/install.sh | sh
  elif command -v wget >/dev/null 2>&1; then
    wget -qO- https://astral.sh/uv/install.sh | sh
  else
    die "need curl or wget to bootstrap uv automatically"
  fi

  export PATH="$HOME/.local/bin:$PATH"
  require_cmd uv
}

reset_python_state_if_needed() {
  mkdir -p "$STATE_DIR"

  if [[ ! -x "$VENV_DIR/bin/python" ]]; then
    if stage_done python_env || stage_done python_deps; then
      warn "python stages marked done but venv python is missing; resetting python stages"
      clear_stage_done python_env python_deps
    fi
    return 0
  fi

  if ! python_is_supported_for_repo "$VENV_DIR/bin/python"; then
    warn "existing virtualenv uses unsupported Python: $("$VENV_DIR/bin/python" --version 2>/dev/null || echo unknown)"
    warn "removing incompatible virtualenv and resetting python stages"
    rm -rf "$VENV_DIR"
    clear_stage_done python_env python_deps
    return 0
  fi

  if ! python_is_preferred "$VENV_DIR/bin/python"; then
    warn "existing virtualenv is compatible but not preferred: $("$VENV_DIR/bin/python" --version)"
    warn "keeping it; preferred version is Python $PYTHON_SPEC"
  fi
}

write_activation_script() {
  mkdir -p "$STATE_DIR"
  cat > "$STATE_DIR/activate.sh" <<EOF
#!/usr/bin/env bash
export HB_DATASET="$HB_DATASET"
export HS_DATASET="$HS_DATASET"
export HIGGSTOOLS_PREFIX="$HIGGSTOOLS_PREFIX"
source "$VENV_DIR/bin/activate"
EOF
  chmod +x "$STATE_DIR/activate.sh"
}

stage_project_config() {
  if [[ -f "$PROJECT_CONFIG_FILE" ]]; then
    log "project config already exists: $PROJECT_CONFIG_FILE"
    return 0
  fi

  log "creating project config at $PROJECT_CONFIG_FILE"
  cat > "$PROJECT_CONFIG_FILE" <<EOF
{
  "data_lake_dir": "$DATA_LAKE_DIR"
}
EOF
}

stage_python_env() {
  local chosen_python=""
  local python_label=""

  reset_python_state_if_needed

  if command -v "$PYTHON_BIN" >/dev/null 2>&1; then
    chosen_python="$(command -v "$PYTHON_BIN")"
    if python_is_supported_for_repo "$chosen_python"; then
      python_label="$("$chosen_python" --version 2>/dev/null || echo "$chosen_python")"
      log "using system python for venv: $python_label"
    else
      warn "system python is unsupported for this repo: $("$chosen_python" --version 2>/dev/null || echo "$chosen_python")"
      chosen_python=""
    fi
  else
    warn "configured python binary not found: $PYTHON_BIN"
  fi

  if [[ -z "$chosen_python" ]]; then
    bootstrap_uv
  fi

  if [[ ! -d "$VENV_DIR" ]]; then
    if [[ -n "$chosen_python" ]]; then
      log "creating virtual environment at $VENV_DIR"
      "$chosen_python" -m venv "$VENV_DIR"
    else
      log "creating virtual environment at $VENV_DIR with Python $PYTHON_SPEC via uv"
      uv venv --python "$PYTHON_SPEC" "$VENV_DIR"
    fi
  fi

  [[ -x "$VENV_DIR/bin/python" ]] || die "virtual environment python not available at $VENV_DIR/bin/python"

  if ! python_is_supported_for_repo "$VENV_DIR/bin/python"; then
    die "virtual environment was created with unsupported Python: $("$VENV_DIR/bin/python" --version)"
  fi

  log "venv python: $("$VENV_DIR/bin/python" --version)"
  "$VENV_DIR/bin/python" -m pip install --upgrade pip setuptools wheel
}

stage_python_deps() {
  [[ -f "$REQ_FILE" ]] || die "missing requirements file: $REQ_FILE"

  log "installing python dependencies from $REQ_FILE"
  PIP_DISABLE_PIP_VERSION_CHECK=1 \
    "$VENV_DIR/bin/python" -m pip install --prefer-binary -r "$REQ_FILE"
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

  if [[ ! -d "$HIGGSTOOLS_SRC/.git" ]]; then
    if [[ -d "$HIGGSTOOLS_SRC" && ! -d "$HIGGSTOOLS_SRC/.git" ]]; then
      warn "existing $HIGGSTOOLS_SRC is not a git checkout; reusing directory as-is"
    else
      log "cloning HiggsTools into $HIGGSTOOLS_SRC"
      git clone "$HIGGSTOOLS_URL" "$HIGGSTOOLS_SRC"
    fi
  fi

  mkdir -p "$HIGGSTOOLS_PREFIX"

  log "configuring HiggsTools"
  cmake -S "$HIGGSTOOLS_SRC" -B "$HIGGSTOOLS_SRC/build" \
    -DCMAKE_INSTALL_PREFIX="$HIGGSTOOLS_PREFIX" \
    -DHiggsTools_BUILD_TESTING=OFF \
    -DHiggsTools_BUILD_EXAMPLES=OFF

  local jobs
  jobs="$(num_jobs)"
  log "building HiggsTools with $jobs parallel jobs"
  cmake --build "$HIGGSTOOLS_SRC/build" --parallel "$jobs"
  cmake --install "$HIGGSTOOLS_SRC/build"

  [[ -d "$HIGGSTOOLS_PREFIX/include/Higgs" ]] || die "HiggsTools headers were not installed"
  find "$HIGGSTOOLS_PREFIX/lib" -maxdepth 1 -name 'libHiggsTools*' | grep -q . || die "HiggsTools library was not installed"
}

stage_datasets() {
  require_cmd git
  [[ -f "$GET_DATASETS_SCRIPT" ]] || die "missing dataset helper script: $GET_DATASETS_SCRIPT"
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

  require_cmd git
  require_cmd make
  require_cmd g++
  require_cmd cmake

  reset_python_state_if_needed

  run_stage project_config stage_project_config
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