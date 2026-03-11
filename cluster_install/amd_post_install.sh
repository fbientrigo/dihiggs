#!/usr/bin/env bash
set -Eeuo pipefail

trap 'echo "==== ERROR ====" >&2; echo "Linea: $LINENO" >&2; exit 1' ERR

echo "==== DIHIGGS AMD POST-INSTALL ===="

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
HOST_ADMIN_RUNBOOKS="${HOST_ADMIN_RUNBOOKS:-$HOST_SHARED_BASE/admin/runbooks}"

CONTAINER_NAME="${CONTAINER_NAME:-dihiggs-dev}"
IMAGE_NAME="${IMAGE_NAME:-dihiggs:local}"

CONTAINER_REPO_ROOT="${CONTAINER_REPO_ROOT:-/workspaces/dihiggs}"
CONTAINER_RESULTS_ROOT="${CONTAINER_RESULTS_ROOT:-/shared-results}"
CONTAINER_DATA_ROOT="${CONTAINER_DATA_ROOT:-/shared-data}"

SSH_CLIENT_ALIVE_INTERVAL="${SSH_CLIENT_ALIVE_INTERVAL:-120}"
SSH_CLIENT_ALIVE_COUNT_MAX="${SSH_CLIENT_ALIVE_COUNT_MAX:-10}"

SSH_DROPIN="${SSH_DROPIN:-/etc/ssh/sshd_config.d/99-dihiggs-keepalive.conf}"

STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
POST_INSTALL_LOG="$HOST_SHARED_BASE/admin/logs/${STAMP}_amd_post_install.log"

mkdir -p "$HOST_SHARED_BASE/admin/logs"
mkdir -p "$HOST_ADMIN_RUNBOOKS"

exec > >(tee -a "$POST_INSTALL_LOG") 2>&1

echo "==== CONFIG ===="
echo "HOST_REPO_ROOT=$HOST_REPO_ROOT"
echo "HOST_SHARED_BASE=$HOST_SHARED_BASE"
echo "HOST_ADMIN_RUNBOOKS=$HOST_ADMIN_RUNBOOKS"
echo "CONTAINER_NAME=$CONTAINER_NAME"
echo "IMAGE_NAME=$IMAGE_NAME"
echo "CONTAINER_REPO_ROOT=$CONTAINER_REPO_ROOT"
echo "CONTAINER_RESULTS_ROOT=$CONTAINER_RESULTS_ROOT"
echo "CONTAINER_DATA_ROOT=$CONTAINER_DATA_ROOT"
echo "SSH_DROPIN=$SSH_DROPIN"
echo "SSH_CLIENT_ALIVE_INTERVAL=$SSH_CLIENT_ALIVE_INTERVAL"
echo "SSH_CLIENT_ALIVE_COUNT_MAX=$SSH_CLIENT_ALIVE_COUNT_MAX"
echo "POST_INSTALL_LOG=$POST_INSTALL_LOG"

# ------------------------------------------------------------------
# 2) VALIDACIONES BASICAS
# ------------------------------------------------------------------
echo "==== VALIDACIONES BASICAS ===="

command -v docker >/dev/null
command -v systemctl >/dev/null
command -v apt-get >/dev/null

docker --version
docker info >/dev/null

if ! docker image inspect "$IMAGE_NAME" >/dev/null 2>&1; then
  echo "==== FALLO ===="
  echo "No existe la imagen $IMAGE_NAME."
  echo "Ejecuta primero: bash cluster_install/amd.sh"
  exit 1
fi

if ! docker container inspect "$CONTAINER_NAME" >/dev/null 2>&1; then
  echo "==== FALLO ===="
  echo "No existe el contenedor $CONTAINER_NAME."
  echo "Ejecuta primero: bash cluster_install/amd.sh"
  exit 1
fi

echo "==== CONTENEDOR ACTUAL ===="
docker ps -a --filter "name=$CONTAINER_NAME" --format 'table {{.Names}}\t{{.Status}}\t{{.Image}}\t{{.RunningFor}}'
docker inspect -f 'RestartPolicy={{.HostConfig.RestartPolicy.Name}}' "$CONTAINER_NAME"

# ------------------------------------------------------------------
# 3) PAQUETES DE UTILIDAD EN EL HOST
# ------------------------------------------------------------------
echo "==== INSTALAR HERRAMIENTAS DE HOST ===="

apt-get update
apt-get install -y \
  python3.12-venv \
  tmux \
  htop \
  tree \
  jq \
  rsync \
  ripgrep \
  less \
  lsof \
  netcat-openbsd \
  unzip \
  zip

# ------------------------------------------------------------------
# 4) ENDURECER COMODIDAD SSH (KEEPALIVE)
# ------------------------------------------------------------------
echo "==== CONFIGURAR SSH KEEPALIVE ===="

mkdir -p "$(dirname "$SSH_DROPIN")"

cat > "$SSH_DROPIN" <<EOF
ClientAliveInterval $SSH_CLIENT_ALIVE_INTERVAL
ClientAliveCountMax $SSH_CLIENT_ALIVE_COUNT_MAX
TCPKeepAlive yes
EOF

if systemctl is-active --quiet ssh; then
  systemctl reload ssh || systemctl restart ssh
elif systemctl is-active --quiet sshd; then
  systemctl reload sshd || systemctl restart sshd
else
  echo "==== AVISO ===="
  echo "No se encontro servicio ssh/sshd activo para recargar."
fi

echo "==== SSH DROP-IN ACTIVO ===="
cat "$SSH_DROPIN"

# ------------------------------------------------------------------
# 5) ASEGURAR CONTENEDOR LEVANTADO Y PERSISTENTE
# ------------------------------------------------------------------
echo "==== ASEGURAR CONTENEDOR LEVANTADO ===="

docker update --restart unless-stopped "$CONTAINER_NAME" >/dev/null
docker start "$CONTAINER_NAME" >/dev/null || true

docker ps --filter "name=$CONTAINER_NAME" --format 'table {{.Names}}\t{{.Status}}\t{{.Image}}\t{{.RunningFor}}'

# ------------------------------------------------------------------
# 6) HELPERS DE OPERACION EN EL HOST
# ------------------------------------------------------------------
echo "==== INSTALAR HELPERS ===="

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

cat > /usr/local/bin/dihiggs-paths <<EOF
#!/usr/bin/env bash
cat <<TXT
CONTAINER_NAME=$CONTAINER_NAME
HOST_REPO_ROOT=$HOST_REPO_ROOT
HOST_SHARED_BASE=$HOST_SHARED_BASE
CONTAINER_REPO_ROOT=$CONTAINER_REPO_ROOT
CONTAINER_RESULTS_ROOT=$CONTAINER_RESULTS_ROOT
CONTAINER_DATA_ROOT=$CONTAINER_DATA_ROOT
TXT
EOF

cat > /usr/local/bin/dihiggs-tmux <<EOF
#!/usr/bin/env bash
set -euo pipefail
SESSION="dihiggs-dev"

if ! docker ps --format '{{.Names}}' | grep -qx '$CONTAINER_NAME'; then
  echo "Contenedor $CONTAINER_NAME no esta corriendo. Intentando iniciarlo..."
  docker start '$CONTAINER_NAME' >/dev/null
fi

if tmux has-session -t "\$SESSION" 2>/dev/null; then
  exec tmux attach -t "\$SESSION"
else
  exec tmux new-session -s "\$SESSION" "docker exec -it $CONTAINER_NAME bash"
fi
EOF

cat > /usr/local/bin/dihiggs-vscode-info <<EOF
#!/usr/bin/env bash
set -euo pipefail

PUBLIC_IP=\$(hostname -I | awk '{print \$1}')

cat <<TXT
=== Flujo recomendado para VS Code ===

1) Desde tu maquina local, agrega algo como esto a ~/.ssh/config:

Host dihiggs-do
    HostName \$PUBLIC_IP
    User root
    Port 22
    ServerAliveInterval 30
    ServerAliveCountMax 6
    TCPKeepAlive yes

2) En VS Code local:
   - instala "Remote - SSH"
   - instala "Dev Containers"

3) Conectate por:
   Remote-SSH: Connect to Host... -> dihiggs-do

4) Ya dentro de la ventana remota:
   Dev Containers: Attach to Running Container...

5) Elige el contenedor:
   $CONTAINER_NAME

6) Abre la carpeta:
   $CONTAINER_REPO_ROOT

Rutas utiles:
- codigo:      $CONTAINER_REPO_ROOT
- resultados:  $CONTAINER_RESULTS_ROOT
- datos:       $CONTAINER_DATA_ROOT
TXT
EOF

chmod +x \
  /usr/local/bin/dihiggs-status \
  /usr/local/bin/dihiggs-enter \
  /usr/local/bin/dihiggs-start \
  /usr/local/bin/dihiggs-stop \
  /usr/local/bin/dihiggs-logs \
  /usr/local/bin/dihiggs-paths \
  /usr/local/bin/dihiggs-tmux \
  /usr/local/bin/dihiggs-vscode-info

echo "==== PRUEBA HELPERS ===="
dihiggs-status
dihiggs-paths

# ------------------------------------------------------------------
# 7) ALIASES MINIMOS PARA ROOT
# ------------------------------------------------------------------
echo "==== CONFIGURAR ALIASES EN /root/.bashrc ===="

if ! grep -q "DIHIGGS_HELPERS_BEGIN" /root/.bashrc 2>/dev/null; then
  cat >> /root/.bashrc <<'EOF'

# DIHIGGS_HELPERS_BEGIN
alias ds='dihiggs-status'
alias de='dihiggs-enter'
alias dt='dihiggs-tmux'
alias dl='dihiggs-logs'
alias dp='dihiggs-paths'
# DIHIGGS_HELPERS_END
EOF
fi

# ------------------------------------------------------------------
# 8) RUNBOOKS Y ARCHIVOS DE AYUDA
# ------------------------------------------------------------------
echo "==== CREAR RUNBOOKS ===="

PUBLIC_IP="$(hostname -I | awk '{print $1}')"

cat > "$HOST_ADMIN_RUNBOOKS/LOCAL_SSH_CONFIG.example" <<EOF
Host dihiggs-do
    HostName $PUBLIC_IP
    User root
    Port 22
    ServerAliveInterval 30
    ServerAliveCountMax 6
    TCPKeepAlive yes
EOF

cat > "$HOST_ADMIN_RUNBOOKS/VSCODE_ATTACH.md" <<EOF
# VS Code Attach

## Requisitos en la maquina local

- Extension: Remote - SSH
- Extension: Dev Containers

## Flujo

1. Conectarse por Remote-SSH al host
2. En la ventana remota: Attach to Running Container
3. Elegir: $CONTAINER_NAME
4. Abrir: $CONTAINER_REPO_ROOT

## Rutas utiles

- codigo: $CONTAINER_REPO_ROOT
- resultados: $CONTAINER_RESULTS_ROOT
- datos: $CONTAINER_DATA_ROOT
EOF

cat > "$HOST_ADMIN_RUNBOOKS/HOST_OPERATIONS.md" <<EOF
# Host operations

## Helpers

- dihiggs-status
- dihiggs-enter
- dihiggs-start
- dihiggs-stop
- dihiggs-logs
- dihiggs-paths
- dihiggs-tmux
- dihiggs-vscode-info

## Contenedor oficial

- nombre: $CONTAINER_NAME
- imagen: $IMAGE_NAME
EOF

# ------------------------------------------------------------------
# 9) RESUMEN FINAL
# ------------------------------------------------------------------
echo "==== RESUMEN FINAL ===="

cat > "$HOST_SHARED_BASE/admin/POST_INSTALL_SUMMARY.txt" <<EOF
Estado: OK
Fecha UTC: $STAMP
Contenedor: $CONTAINER_NAME
Imagen: $IMAGE_NAME
Keepalive SSH: ClientAliveInterval=$SSH_CLIENT_ALIVE_INTERVAL ClientAliveCountMax=$SSH_CLIENT_ALIVE_COUNT_MAX
Runbooks: $HOST_ADMIN_RUNBOOKS
Log post-install: $POST_INSTALL_LOG
EOF

cat "$HOST_SHARED_BASE/admin/POST_INSTALL_SUMMARY.txt"

echo "==== POST-INSTALL COMPLETO ===="
echo "Entrar al contenedor: dihiggs-enter"
echo "Entrar via tmux:      dihiggs-tmux"
echo "Ver pasos VS Code:    dihiggs-vscode-info"
echo "Runbooks en:          $HOST_ADMIN_RUNBOOKS"
echo "Log en:               $POST_INSTALL_LOG"