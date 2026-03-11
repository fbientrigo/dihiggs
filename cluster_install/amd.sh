#!/usr/bin/env bash
set -Eeuo pipefail

trap 'echo "==== ERROR ====" >&2; echo "Linea: $LINENO" >&2; exit 1' ERR

echo "==== DIHIGGS AMD INSTALLER ===="

# ------------------------------------------------------------------
# 0) ELEVAR PRIVILEGIOS SI HACE FALTA
# ------------------------------------------------------------------
if [ "${EUID:-$(id -u)}" -ne 0 ]; then
  exec sudo -E bash "$0" "$@"
fi

# ------------------------------------------------------------------
# 1) CONFIGURACION BASE
# ------------------------------------------------------------------
HOST_REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"

HOST_SHARED_BASE="${HOST_SHARED_BASE:-/shared-docker}"
HOST_RESULTS_ROOT="${HOST_RESULTS_ROOT:-$HOST_SHARED_BASE/results/dihiggs}"
HOST_DATA_ROOT="${HOST_DATA_ROOT:-$HOST_SHARED_BASE/data/dihiggs}"

IMAGE_NAME="${IMAGE_NAME:-dihiggs:local}"
CONTAINER_NAME="${CONTAINER_NAME:-dihiggs-dev}"

CONTAINER_REPO_ROOT="${CONTAINER_REPO_ROOT:-/workspaces/dihiggs}"
CONTAINER_RESULTS_ROOT="${CONTAINER_RESULTS_ROOT:-/shared-results}"
CONTAINER_DATA_ROOT="${CONTAINER_DATA_ROOT:-/shared-data}"

VSCODE_UID="${VSCODE_UID:-1000}"
VSCODE_GID="${VSCODE_GID:-1000}"

STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
RUN_NAME="${RUN_NAME:-${STAMP}_smoke}"
HOST_RUNROOT="$HOST_RESULTS_ROOT/runs/$RUN_NAME"
CONTAINER_RUNROOT="$CONTAINER_RESULTS_ROOT/runs/$RUN_NAME"

echo "==== CONFIG ===="
echo "HOST_REPO_ROOT=$HOST_REPO_ROOT"
echo "HOST_SHARED_BASE=$HOST_SHARED_BASE"
echo "HOST_RESULTS_ROOT=$HOST_RESULTS_ROOT"
echo "HOST_DATA_ROOT=$HOST_DATA_ROOT"
echo "IMAGE_NAME=$IMAGE_NAME"
echo "CONTAINER_NAME=$CONTAINER_NAME"
echo "CONTAINER_REPO_ROOT=$CONTAINER_REPO_ROOT"
echo "CONTAINER_RESULTS_ROOT=$CONTAINER_RESULTS_ROOT"
echo "CONTAINER_DATA_ROOT=$CONTAINER_DATA_ROOT"
echo "VSCODE_UID=$VSCODE_UID"
echo "VSCODE_GID=$VSCODE_GID"
echo "RUN_NAME=$RUN_NAME"

# ------------------------------------------------------------------
# 2) VALIDACIONES BASICAS
# ------------------------------------------------------------------
echo "==== VALIDACIONES BASICAS ===="

command -v docker >/dev/null
command -v git >/dev/null

docker --version
docker info >/dev/null

if [ ! -d "$HOST_REPO_ROOT/.git" ]; then
  echo "==== FALLO ===="
  echo "El script debe ejecutarse desde un repo clonado de dihiggs."
  echo "No se encontro .git en: $HOST_REPO_ROOT"
  exit 1
fi

if [ ! -f "$HOST_REPO_ROOT/.devcontainer/Dockerfile" ]; then
  echo "==== FALLO ===="
  echo "No se encontro .devcontainer/Dockerfile en: $HOST_REPO_ROOT"
  exit 1
fi

GIT_BRANCH="$(git -C "$HOST_REPO_ROOT" rev-parse --abbrev-ref HEAD)"
GIT_COMMIT="$(git -C "$HOST_REPO_ROOT" rev-parse HEAD)"

echo "==== REPO ===="
echo "Branch: $GIT_BRANCH"
echo "Commit: $GIT_COMMIT"

# ------------------------------------------------------------------
# 3) PREPARAR DIRECTORIOS PERSISTENTES
# ------------------------------------------------------------------
echo "==== PREPARAR DIRECTORIOS PERSISTENTES ===="

mkdir -p "$HOST_RESULTS_ROOT/runs"
mkdir -p "$HOST_RESULTS_ROOT/logs"
mkdir -p "$HOST_DATA_ROOT"
mkdir -p "$HOST_SHARED_BASE/admin/runbooks"

# Importante: el repo ya existe; los mounts deben ser escribibles por vscode:1000
chown -R "$VSCODE_UID:$VSCODE_GID" "$HOST_REPO_ROOT"
chown -R "$VSCODE_UID:$VSCODE_GID" "$HOST_RESULTS_ROOT"
chown -R "$VSCODE_UID:$VSCODE_GID" "$HOST_DATA_ROOT"

echo "==== OWNERSHIP ===="
ls -ld "$HOST_REPO_ROOT"
ls -ld "$HOST_RESULTS_ROOT"
ls -ld "$HOST_DATA_ROOT"

# ------------------------------------------------------------------
# 4) BUILD DE IMAGEN DESDE DEVCONTAINER
# ------------------------------------------------------------------
echo "==== BUILD DE IMAGEN LOCAL ===="

docker build \
  -t "$IMAGE_NAME" \
  -f "$HOST_REPO_ROOT/.devcontainer/Dockerfile" \
  "$HOST_REPO_ROOT"

IMAGE_ID="$(docker image inspect "$IMAGE_NAME" --format '{{.Id}}')"

echo "==== IMAGEN CONSTRUIDA ===="
echo "IMAGE_NAME=$IMAGE_NAME"
echo "IMAGE_ID=$IMAGE_ID"

# ------------------------------------------------------------------
# 5) RECREAR CONTENEDOR PERSISTENTE
# ------------------------------------------------------------------
echo "==== RECREAR CONTENEDOR PERSISTENTE ===="

docker rm -f "$CONTAINER_NAME" 2>/dev/null || true

docker run -d \
  --name "$CONTAINER_NAME" \
  --restart unless-stopped \
  --user "$VSCODE_UID:$VSCODE_GID" \
  -v "$HOST_REPO_ROOT:$CONTAINER_REPO_ROOT" \
  -v "$HOST_RESULTS_ROOT:$CONTAINER_RESULTS_ROOT" \
  -v "$HOST_DATA_ROOT:$CONTAINER_DATA_ROOT" \
  -w "$CONTAINER_REPO_ROOT" \
  --entrypoint bash \
  "$IMAGE_NAME" \
  -lc "sleep infinity"

echo "==== CONTENEDOR ACTIVO ===="
docker ps --filter "name=$CONTAINER_NAME" --format 'table {{.Names}}\t{{.Status}}\t{{.Image}}\t{{.RunningFor}}'

echo "==== RESTART POLICY ===="
docker inspect -f '{{.HostConfig.RestartPolicy.Name}}' "$CONTAINER_NAME"

# ------------------------------------------------------------------
# 6) SANITY CHECK DE ESCRITURA
# ------------------------------------------------------------------
echo "==== SANITY CHECK DE ESCRITURA ===="

docker exec "$CONTAINER_NAME" bash -lc "
set -Eeuo pipefail
whoami
id
mkdir -p '$CONTAINER_RESULTS_ROOT/runs/perm_test'
touch '$CONTAINER_RESULTS_ROOT/runs/perm_test/ok.txt'
mkdir -p '$CONTAINER_DATA_ROOT/perm_test'
touch '$CONTAINER_DATA_ROOT/perm_test/ok.txt'
touch '$CONTAINER_REPO_ROOT/.perm_test'
rm -f '$CONTAINER_REPO_ROOT/.perm_test'
ls -lah '$CONTAINER_RESULTS_ROOT/runs/perm_test'
ls -lah '$CONTAINER_DATA_ROOT/perm_test'
"

# ------------------------------------------------------------------
# 7) MANIFIESTO DEL RUN
# ------------------------------------------------------------------
echo "==== CREAR MANIFIESTO DEL RUN ===="

mkdir -p "$HOST_RUNROOT"

cat > "$HOST_RUNROOT/00_manifest_host.txt" <<EOF
RUN_NAME=$RUN_NAME
TIMESTAMP_UTC=$STAMP
GIT_BRANCH=$GIT_BRANCH
GIT_COMMIT=$GIT_COMMIT
HOST_REPO_ROOT=$HOST_REPO_ROOT
HOST_RESULTS_ROOT=$HOST_RESULTS_ROOT
HOST_DATA_ROOT=$HOST_DATA_ROOT
IMAGE_NAME=$IMAGE_NAME
IMAGE_ID=$IMAGE_ID
CONTAINER_NAME=$CONTAINER_NAME
CONTAINER_REPO_ROOT=$CONTAINER_REPO_ROOT
CONTAINER_RESULTS_ROOT=$CONTAINER_RESULTS_ROOT
CONTAINER_DATA_ROOT=$CONTAINER_DATA_ROOT
VSCODE_UID=$VSCODE_UID
VSCODE_GID=$VSCODE_GID
EOF

# ------------------------------------------------------------------
# 8) BUILD Y TESTS DENTRO DEL CONTENEDOR
# ------------------------------------------------------------------
echo "==== BUILD Y TESTS DENTRO DEL CONTENEDOR ===="

docker exec \
  -e REPO_ROOT="$CONTAINER_REPO_ROOT" \
  -e DATA_ROOT="$CONTAINER_DATA_ROOT" \
  -e HB_DATASET="$CONTAINER_DATA_ROOT/HBDataset" \
  -e HS_DATASET="$CONTAINER_DATA_ROOT/HSDataset" \
  -e HIGGSTOOLS_PREFIX="$CONTAINER_REPO_ROOT/higgstools/installation" \
  -e RUNROOT="$CONTAINER_RUNROOT" \
  "$CONTAINER_NAME" \
  bash -lc '
set -Eeuo pipefail

echo "==== PATHS ACTIVOS ====" | tee "$RUNROOT/00_paths.log"
echo "REPO_ROOT=$REPO_ROOT" | tee -a "$RUNROOT/00_paths.log"
echo "DATA_ROOT=$DATA_ROOT" | tee -a "$RUNROOT/00_paths.log"
echo "HB_DATASET=$HB_DATASET" | tee -a "$RUNROOT/00_paths.log"
echo "HS_DATASET=$HS_DATASET" | tee -a "$RUNROOT/00_paths.log"
echo "HIGGSTOOLS_PREFIX=$HIGGSTOOLS_PREFIX" | tee -a "$RUNROOT/00_paths.log"

echo "==== 1) BUILD 2HDMC ===="
(
  cd "$REPO_ROOT/2hdmc"
  make clean || true
  make -j"$(nproc)"
) 2>&1 | tee "$RUNROOT/01_build_2hdmc.log"

echo "==== 2) BUILD HIGGSTOOLS ===="
rm -rf "$REPO_ROOT/higgstools/build" "$REPO_ROOT/higgstools/installation"
mkdir -p "$REPO_ROOT/higgstools/build"

(
  cd "$REPO_ROOT/higgstools/build"
  cmake .. -DCMAKE_INSTALL_PREFIX="$HIGGSTOOLS_PREFIX"
  make -j"$(nproc)"
  make install
) 2>&1 | tee "$RUNROOT/02_build_higgstools.log"

echo "==== 3) DATASETS ===="

sync_dataset() {
  local dir="$1"
  if [ -d "$dir/.git" ]; then
    echo "Actualizando dataset en $dir"
    git -C "$dir" fetch --all --tags || true
    git -C "$dir" pull --ff-only || true
    return 0
  fi
  return 1
}

if sync_dataset "$HB_DATASET" && sync_dataset "$HS_DATASET"; then
  echo "Datasets ya existentes: actualizados o verificados." | tee "$RUNROOT/03_get_datasets.log"
else
  (
    cd "$REPO_ROOT"
    bash ./get_datasets.sh "$DATA_ROOT"
  ) 2>&1 | tee "$RUNROOT/03_get_datasets.log"
fi

echo "==== 4) VALIDAR DATASETS ===="
ls -lah "$DATA_ROOT" | tee "$RUNROOT/04_datasets_ls.log"
ls -lah "$HB_DATASET" | tee -a "$RUNROOT/04_datasets_ls.log"
ls -lah "$HS_DATASET" | tee -a "$RUNROOT/04_datasets_ls.log"

echo "==== 5) BUILD CORE DIHIGGS ===="
(
  cd "$REPO_ROOT/dihiggs"
  make clean || true
  make HIGGSTOOLS_PREFIX="$HIGGSTOOLS_PREFIX" lib/libdihiggs.a
  make HIGGSTOOLS_PREFIX="$HIGGSTOOLS_PREFIX" all
) 2>&1 | tee "$RUNROOT/05_build_core.log"

echo "==== 6) BUILD TESTS ===="
(
  cd "$REPO_ROOT/tests/unit"
  make clean || true
  make dihiggs_unit_tests
) 2>&1 | tee "$RUNROOT/06_build_tests.log"

echo "==== 7) RUN TESTS ===="
(
  cd "$REPO_ROOT/tests/unit"
  ./dihiggs_unit_tests -s
) 2>&1 | tee "$RUNROOT/07_run_tests.log"

echo "==== 8) INVENTARIO FINAL ===="
find "$REPO_ROOT/dihiggs" -maxdepth 3 \( -type f -perm -111 -o -name "*.a" \) | sort | tee "$RUNROOT/08_binaries.txt"
'

# ------------------------------------------------------------------
# 9) HELPERS EN EL HOST
# ------------------------------------------------------------------
echo "==== INSTALAR HELPERS EN HOST ===="

cat > /usr/local/bin/dihiggs-status <<EOF
#!/usr/bin/env bash
docker ps -a --filter name=$CONTAINER_NAME --format 'table {{.Names}}\t{{.Status}}\t{{.Image}}\t{{.RunningFor}}'
EOF

cat > /usr/local/bin/dihiggs-enter <<EOF
#!/usr/bin/env bash
exec docker exec -it $CONTAINER_NAME bash
EOF

cat > /usr/local/bin/dihiggs-start <<EOF
#!/usr/bin/env bash
docker start $CONTAINER_NAME
docker ps --filter name=$CONTAINER_NAME
EOF

cat > /usr/local/bin/dihiggs-stop <<EOF
#!/usr/bin/env bash
docker stop $CONTAINER_NAME
EOF

cat > /usr/local/bin/dihiggs-logs <<EOF
#!/usr/bin/env bash
docker logs --tail 200 $CONTAINER_NAME
EOF

chmod +x /usr/local/bin/dihiggs-status /usr/local/bin/dihiggs-enter /usr/local/bin/dihiggs-start /usr/local/bin/dihiggs-stop /usr/local/bin/dihiggs-logs

# ------------------------------------------------------------------
# 10) MARCAS DE EXITO
# ------------------------------------------------------------------
echo "==== MARCAS DE EXITO ===="

ln -sfn "$HOST_RUNROOT" "$HOST_RESULTS_ROOT/runs/latest_smoke"
echo "$HOST_RUNROOT" > "$HOST_RESULTS_ROOT/LAST_SUCCESS_RUN.txt"

cat > "$HOST_RESULTS_ROOT/INSTALL_SUMMARY.txt" <<EOF
Estado: OK
Ultimo run exitoso: $HOST_RUNROOT
Contenedor: $CONTAINER_NAME
Imagen: $IMAGE_NAME
Commit: $GIT_COMMIT
Branch: $GIT_BRANCH
EOF

echo "==== INSTALACION COMPLETA ===="
echo "Run root: $HOST_RUNROOT"
echo "Contenedor persistente: $CONTAINER_NAME"
echo "Entrar luego con: dihiggs-enter"
echo "Ver estado con: dihiggs-status"