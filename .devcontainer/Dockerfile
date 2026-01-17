FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

# 1. Instalación de Dependencias del Sistema
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

# 2. Configuración del usuario 'ubuntu' (Ya existe en la imagen base)
# Le damos permisos de sudo sin contraseña para comodidad en desarrollo
RUN echo "ubuntu ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/ubuntu \
    && chmod 0440 /etc/sudoers.d/ubuntu

# 3. Establecemos el usuario por defecto
USER ubuntu
WORKDIR /workspaces/dihiggs

# 4. Preparamos el entorno Python
# Copiamos requirements.txt y cambiamos el dueño al usuario 'ubuntu'
COPY --chown=ubuntu:ubuntu requirements.txt /tmp/pip-tmp/

# Instalamos las dependencias
# Nota: --break-system-packages es necesario en Ubuntu 24.04 para pip global
# o usamos un venv. Aquí, por simplicidad en devcontainer, usamos la flag.
RUN pip3 install --no-cache-dir --break-system-packages -r /tmp/pip-tmp/requirements.txt \
    && rm -rf /tmp/pip-tmp