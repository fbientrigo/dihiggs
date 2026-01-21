#!/bin/bash
set -e # Detener si hay error

echo "Iniciando configuración post-creación..."

# 1. Compilar 2HDMC (Paso 3 del README)
echo "Compilando 2HDMC..."
if [ -d "2hdmc" ]; then
    cd 2hdmc
    make clean
    make
    cd ..
else
    echo "Carpeta 2hdmc no encontrada."
fi

# 2. Instalar HiggsTools (Paso 4 del README)
echo "Configurando HiggsTools..."
# Asumimos que la carpeta higgstools ya viene en tu repo según lo que subiste
if [ -d "higgstools" ]; then
    cd higgstools
    # Limpiamos builds previos por si acaso
    rm -rf build
    mkdir -p build && cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=../installation
    make -j$(nproc)
    sudo make install
    cd ../..
else
    echo "!!  Carpeta higgstools no encontrada. Clonando..."
    git clone https://gitlab.com/higgstools/higgstools.git
    # (Aquí podrías repetir los pasos de compilación si se clona de cero)
fi

# 3. Descargar Datasets (Paso 5 del README)
echo "Descargando Datasets..."
if [ -f "get_datasets.sh" ]; then
    chmod +x get_datasets.sh
    ./get_datasets.sh
fi

# 4. Compilación final del proyecto (Paso 5 del README)
echo "Compilando proyecto principal..."
if [ -d "dihiggs" ]; then     # Verificación defensiva
    cd dihiggs                # <--- FIX: Entrar al directorio
    make clean                # Buena práctica para limpiar builds previos
    make all                     # Ahora sí encontrará el Makefile
    cd ..                     # Volver a la raíz (opcional, pero limpio)
else
    echo "Error: No se encuentra la carpeta 'dihiggs' para compilar."
    exit 1
fi

echo "Proyectos C++ compilados exitosamente."

echo "--- Iniciando configuración de Python ---"

# 1. Instalar requerimientos si el archivo existe
if [ -f "/tmp/requirements.txt" ]; then
    echo "Instalando dependencias de Python..."
    pip3 install --no-cache-dir --user -r /tmp/requirements.txt
elif [ -f ".devcontainer/requirements.txt" ]; then
    pip3 install --no-cache-dir --user -r .devcontainer/requirements.txt
fi

# 2. Verificar ROOT (Grado científico: validamos el entorno)
if python3 -c "import ROOT; print('ROOT import exitoso')" &> /dev/null; then
    echo "Entorno de ROOT validado correctamente."
else
    echo "ADVERTENCIA: ROOT no se puede importar en Python. Revisa el PYTHONPATH."
fi

echo "--- Configuración completada ---"
echo "¡Entorno listo para la física de partículas!"