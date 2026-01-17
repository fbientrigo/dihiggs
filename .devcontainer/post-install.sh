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
# make clean # Opcional
make

echo "¡Entorno listo para la física de partículas!"