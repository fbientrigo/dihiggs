# Nombre del compilador
CXX = g++

# Prefijo de instalación local de HiggsTools (puedes sobreescribir por CLI)
HIGGSTOOLS_PREFIX ?= $(HOME)/dihiggs/higgstools/installation
HIGGSTOOLS_INC    := $(HIGGSTOOLS_PREFIX)/include
HIGGSTOOLS_LIB    := $(HIGGSTOOLS_PREFIX)/lib

# Flags de compilación
CXXFLAGS     = -I../2hdmc/src -I$(HIGGSTOOLS_INC) -fopenmp -std=gnu++17 -O2 -Wall
SPEED_FLAGS  = -DSPEED_TEST

# Flags de enlazado (ajusta si necesitas más librerías)
LDFLAGS = -L../2hdmc/lib -l2HDMC \
          -L$(HIGGSTOOLS_LIB) -lHiggsTools -lgsl -lgslcblas -lm \
          -Wl,-rpath,$(HIGGSTOOLS_LIB)

# Directorios
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = app

# ============ 
# Archivos fuente con main (cada uno genera un ejecutable).
# ============ 
MAIN_SRCS = \
  GenParamScan.cpp \
  HiggsParamScan.cpp \
  PhysParamScan.cpp \
  testing.cpp \
  Oracle.cpp \
  ModularParamScan.cpp \
  PhysScanWithFixings.cpp \
  PhysScanPolynomial.cpp

# ============ 
# Archivos fuente comunes (sin main).
# ============ 
COMMON_SRCS = \
  ParamUtils.cpp

# Lista completa de fuentes
SRCS = $(MAIN_SRCS:%=$(SRC_DIR)/%) $(COMMON_SRCS:%=$(SRC_DIR)/%)

# Objetos normales
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))

# Objeto especializado para speed test
OBJ_SPEED = $(OBJ_DIR)/PhysScanWithFixings_speed_test.o

# ============ 
# Binarios resultantes
# ============ 
BINARIES = \
  $(BIN_DIR)/GenParamScan \
  $(BIN_DIR)/HiggsParamScan \
  $(BIN_DIR)/PhysParamScan \
  $(BIN_DIR)/testing \
  $(BIN_DIR)/Oracle \
  $(BIN_DIR)/ModularParamScan \
  $(BIN_DIR)/PhysScanWithFixings \
  $(BIN_DIR)/PhysScanPolynomial

# Binario adicional de speed test
SPEED_BINARIES = \
  $(BIN_DIR)/PhysScanWithFixings_speed_test

# Regla principal
all: $(BIN_DIR) $(OBJ_DIR) $(BINARIES) $(SPEED_BINARIES)

# ---- Regla de enlace de cada ejecutable normal ----
$(BIN_DIR)/%: $(OBJ_DIR)/%.o $(OBJ_DIR)/ParamUtils.o
	$(CXX) -fopenmp $^ -o $@ $(LDFLAGS)

# ---- Regla de enlace del speed test ----
$(BIN_DIR)/PhysScanWithFixings_speed_test: $(OBJ_SPEED) $(OBJ_DIR)/ParamUtils.o
	$(CXX) -fopenmp $(OBJ_SPEED) $(OBJ_DIR)/ParamUtils.o -o $@ $(LDFLAGS)

# ---- Compilación genérica de .cpp -> .o ----
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ---- Compilación específica para speed test (mismo .cpp, distinta bandera) ----
$(OBJ_DIR)/PhysScanWithFixings_speed_test.o: $(SRC_DIR)/PhysScanWithFixings.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(SPEED_FLAGS) -c $< -o $@

# Crear directorios si no existen
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Limpieza
clean:
	rm -rf $(OBJ_DIR)/*.o $(BINARIES) $(SPEED_BINARIES)

cleanall: clean
	rm -rf $(BIN_DIR)/*
	rm -rf data/*

.PHONY: all clean cleanall
