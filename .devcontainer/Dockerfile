# Usamos Ubuntu 22.04 (Jammy) para tener Python 3.10 nativo y compatibilidad con Numba
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# 1. Instalación de dependencias del sistema
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    wget \
    libgsl-dev \
    python3 \
    python3-pip \
    python3-venv \
    sudo \
    && rm -rf /var/lib/apt/lists/*

# 2. Creación del usuario 'vscode' (Estándar para DevContainers)
# Como en Ubuntu 22.04 el UID 1000 está libre, lo creamos explícitamente.
ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=$USER_UID

RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME \
    && echo $USERNAME ALL=\(ALL\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

# 3. Configuramos el entorno
USER $USERNAME
WORKDIR /workspaces/dihiggs

# 4. Instalación de librerías Python
COPY --chown=$USERNAME:$USERNAME requirements.txt /tmp/pip-tmp/
# En Ubuntu 22.04 pip es más permisivo, no necesitamos --break-system-packages
RUN pip3 install --no-cache-dir -r /tmp/pip-tmp/requirements.txt \
    && rm -rf /tmp/pip-tmp