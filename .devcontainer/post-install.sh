#!/bin/bash
set -e

echo "--- Iniciando configuración de entorno (Modo: $(whoami)) ---"

# DETECCIÓN DE ENTORNO:
# Si somos root (EUID 0), no usamos sudo. Si somos usuario, usamos sudo.
if [ "$EUID" -eq 0 ]; then
    SUDO_CMD=""
    PIP_FLAGS="--no-cache-dir" # Instalación global
    echo ">> Ejecutando como ROOT (Build de Apptainer/Docker)"
else
    SUDO_CMD="sudo"
    PIP_FLAGS="--no-cache-dir --user" # Instalación local de usuario
    echo ">> Ejecutando como USUARIO (VS Code DevContainer)"
fi

# 1. Compilar 2HDMC
echo ">> Compilando 2HDMC..."
if [ -d "2hdmc" ]; then
    cd 2hdmc
    make clean && make
    cd ..
else
    echo "!! Carpeta 2hdmc no encontrada."
fi

# 2. Instalar HiggsTools
echo ">> Configurando HiggsTools..."
if [ -d "higgstools" ]; then
    cd higgstools
    rm -rf build
    mkdir -p build && cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local
    make -j$(nproc)
    
    # Aquí usamos la variable dinámica SUDO_CMD
    echo "Instalando en /usr/local..."
    $SUDO_CMD make install
    cd ../..
else
    echo "!! Carpeta higgstools no encontrada."
fi

# 3. Descargar Datasets (Opcional, a veces es mejor no meter datos en la imagen)
if [ -f "get_datasets.sh" ]; then
    chmod +x get_datasets.sh
    ./get_datasets.sh
fi

# 4. Compilación del proyecto principal
echo ">> Compilando DiHiggs..."
if [ -d "dihiggs" ]; then
    cd dihiggs
    make clean && make all
    cd ..
fi

echo "--- Configuración de Python ---"

# Instalación de librerías
# Buscamos el requirements en la ruta estándar del repo o en /tmp si Docker lo puso ahí
REQ_FILE="requirements.txt"
[ -f ".devcontainer/requirements.txt" ] && REQ_FILE=".devcontainer/requirements.txt"
[ -f "/tmp/requirements.txt" ] && REQ_FILE="/tmp/requirements.txt"

if [ -f "$REQ_FILE" ]; then
    echo "Instalando dependencias desde $REQ_FILE..."
    # pip3 install respeta los flags definidos arriba (global vs user)
    pip3 install $PIP_FLAGS -r "$REQ_FILE"
fi

# Validación final
if python3 -c "import ROOT" &> /dev/null; then
    echo "✅ ROOT importado correctamente."
else
    echo "⚠️  Advertencia: No se pudo importar ROOT (¿Está configurado PYTHONPATH?)"
fi

echo "✅✅ Configuración completada exitosamente."