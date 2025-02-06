# Nombre del compilador
CXX = g++

# Flags de compilación
CXXFLAGS = -I../2hdmc/src -fopenmp -std=c++11 -O2 -Wall

# Flags de enlazado. Ajusta si necesitas más libs (lapack, etc.) si fuera el caso.
LDFLAGS = -L../2hdmc/lib -l2HDMC -lgsl -lgslcblas -lm

# Directorios
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = app

# ============
# Archivos fuente con main
# Para cada uno, generaremos un ejecutable distinto.
# ============
MAIN_SRCS = \
  ParamScanPhys.cpp \
  ParamScanGen_pll.cpp \
  ParamScanHiggs.cpp \
  testing.cpp

# ============
# Archivos fuente comunes (sin main):
# Por ejemplo ParamUtils.cpp.
# ============
COMMON_SRCS = \
  ParamUtils.cpp

# Generamos la lista completa de fuentes
SRCS = $(MAIN_SRCS:%=$(SRC_DIR)/%) $(COMMON_SRCS:%=$(SRC_DIR)/%)

# Objetos correspondientes
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))

# ============
# Nombre de los binarios resultantes
# Se asume que cada "main" produce uno de estos binarios:
# Si el nombre del .cpp es 'ParamScanGen.cpp', el binario se llamará 'ParamScanGen'
# ============
BINARIES = \
  $(BIN_DIR)/ParamScanPhys \
  $(BIN_DIR)/ParamScanGen_pll \
  $(BIN_DIR)/ParamScanHiggs \
  $(BIN_DIR)/testing

# Regla principal
all: $(BIN_DIR) $(OBJ_DIR) $(BINARIES)

# ----
# Reglas de enlace de cada ejecutable
# Notar que cada ejecutable depende de su .o (con main) + todos los .o de la parte común.
# ----

$(BIN_DIR)/ParamScanPhys: $(OBJ_DIR)/ParamScanPhys.o $(OBJ_DIR)/ParamUtils.o
	$(CXX) -fopenmp $^ -o $@ $(LDFLAGS)

$(BIN_DIR)/ParamScanGen_pll: $(OBJ_DIR)/ParamScanGen_pll.o $(OBJ_DIR)/ParamUtils.o
	$(CXX) -fopenmp $^ -o $@ $(LDFLAGS)

$(BIN_DIR)/ParamScanHiggs: $(OBJ_DIR)/ParamScanHiggs.o $(OBJ_DIR)/ParamUtils.o
	$(CXX) $^ -o $@ $(LDFLAGS)

$(BIN_DIR)/testing: $(OBJ_DIR)/testing.o $(OBJ_DIR)/ParamUtils.o
	$(CXX) $^ -o $@ $(LDFLAGS)

# Regla genérica para compilar .cpp -> .o
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Crear directorios si no existen
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Limpieza
clean:
	rm -rf $(OBJ_DIR)/*.o $(BINARIES)

cleanall: clean
	rm -rf $(BIN_DIR)/*
	rm -rf data/*

.PHONY: all clean cleanall
