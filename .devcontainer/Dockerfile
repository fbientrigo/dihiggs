FROM ubuntu:latest

# Evita que la instalación nos haga preguntas interactivas (como zona horaria)
ENV DEBIAN_FRONTEND=noninteractive

# 2. Instalación de Dependencias del Sistema (Requisitos del README)
# - build-essential: Contiene gcc, g++, make (Vital para paso 3 y 5 del README)
# - cmake: Necesario para HiggsTools (Paso 4)
# - libgsl-dev: La GNU Scientific Library (Requisito explícito)
# - git, wget: Herramientas auxiliares para descargar datasets
# - python3-full, pip: Para el stack de Python
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    wget \
    libgsl-dev \
    python3 \
    python3-pip \
    python3-venv \
    && rm -rf /var/lib/apt/lists/*

# 3. Creamos un usuario no-root (Best Practice)
# Es buena práctica no desarrollar como 'root' para evitar problemas de permisos de archivos
ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=$USER_UID

RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME \
    && apt-get update \
    && apt-get install -y sudo \
    && echo $USERNAME ALL=\(ALL\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

# 4. Establecemos el usuario por defecto
USER $USERNAME
WORKDIR /workspaces/dihiggs

# 5. Preparamos el entorno Python (Paso 1 y 2 del README)
# Instalamos las dependencias globales para el usuario para simplificar
# (Nota: En devcontainers, a veces es más fácil instalar pip user-wide que usar venv por proyecto)
COPY --chown=$USERNAME:$USERNAME requirements.txt /tmp/pip-tmp/
RUN pip3 install --no-cache-dir -r /tmp/pip-tmp/requirements.txt \
    && rm -rf /tmp/pip-tmp