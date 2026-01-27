Vamos a implementar la **Opción A (Docker Archive)** para solucionar tu error de red/DNS, y simultáneamente haremos un **Refactor DRY (Don't Repeat Yourself)** para que `post-install.sh` sea la única fuente de verdad.

el plan de limpieza y ejecución.

### 1. Limpieza del Área de Trabajo

Ejecuta esto en tu terminal (WSL) para borrar intentos fallidos y liberar espacio:

```bash
# Borrar contenedores SIF fallidos o parciales
rm -f *.sif

# Borrar caché de construcción de Apptainer (opcional pero recomendado si hubo errores raros)
apptainer cache clean --force

```

---

### 2. El Nuevo `post-install.sh` (Inteligente)

Vamos a modificar este script para que sea "polimórfico": detectará si está corriendo como **root** (dentro del build de Apptainer) o como **usuario normal** (en VS Code) y ajustará los comandos `sudo` y `pip` automáticamente.

Sobrescribe `.devcontainer/post-install.sh` con esto:

```bash
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

```

*No olvides hacerlo ejecutable:* `chmod +x .devcontainer/post-install.sh`

---

### 3. El Nuevo `apptainer.def` (Simplificado)

Ahora tu archivo de definición es minúsculo. Solo copia los archivos y llama al script maestro. Usamos `docker-archive` para solucionar tu error de red.

Sobrescribe `apptainer.def` con esto:

```singularity
Bootstrap: docker-archive
From: dihiggs-base.tar

%files
    # Copiamos TODO el repo al interior del contenedor en /workspaces/dihiggs
    # Esto incluye .devcontainer, 2hdmc, higgstools, etc.
    . /workspaces/dihiggs

%post
    # Preparar entorno
    export DEBIAN_FRONTEND=noninteractive
    
    # Ir al directorio de trabajo
    cd /workspaces/dihiggs
    
    # Ejecutar el Script Maestro (reutilización de lógica)
    # Al correr en %post, somos ROOT, así que el script instalará globalmente.
    bash .devcontainer/post-install.sh
    
    # Limpieza final
    apt-get clean
    rm -rf /var/lib/apt/lists/*

%environment
    export LC_ALL=C
    export PYTHONPATH=/workspaces/dihiggs:$PYTHONPATH
    # HiggsTools se instala en /usr/local/lib, aseguramos que el linker lo vea
    export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

%runscript
    exec "$@"

%labels
    Author DiHiggs-Team
    Version 2.0-DRY

```

---

### 4. Instrucciones de Ejecución (La Solución "Opción A")

Sigue estos tres pasos exactos en tu terminal WSL (fuera de VS Code):

**Paso 1: Construir la base de Docker**
(Esto ya lo hiciste, pero asegúrate de que esté fresca).

```bash
docker build -f .devcontainer/Dockerfile -t dihiggs-base .

```

**Paso 2: Exportar a archivo (El truco mágico)**
Esto crea un archivo `.tar` físico. Apptainer leerá este archivo en lugar de intentar conectarse al demonio de Docker, eliminando el error de DNS/Timeouts.

```bash
docker save dihiggs-base -o dihiggs-base.tar

```

**Paso 3: Construir el SIF final**
Ahora Apptainer lee el archivo local y ejecuta tu `post-install.sh` dentro.

```bash
# Nota: Ya no necesitamos argumentos complejos porque el .def apunta al archivo tar
apptainer build dihiggs.sif apptainer.def

```

### ¿Por qué esto es mejor?

1. **Cero Redundancia:** Si cambias cómo se compila `2hdmc`, solo editas `post-install.sh`. Apptainer y VS Code se actualizan automáticamente.
2. **Sin Errores de Red:** `docker-archive` funciona 100% offline una vez tienes el tar.
3. **Grado Científico:** El script ahora verifica explícitamente el entorno (ROOT/User) antes de actuar.