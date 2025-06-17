# Nombre del compilador
CXX = g++

# Flags de compilación
CXXFLAGS = -I../2hdmc/src -fopenmp -std=c++11 -O2 -Wall

# Flags de enlazado (ajusta si necesitas más librerías)
LDFLAGS = -L../2hdmc/lib -l2HDMC -lgsl -lgslcblas -lm

# Directorios
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = app

# ============
# Archivos fuente con main (cada uno genera un ejecutable).
# Incluye el nuevo 'ModularParamScan.cpp' y 'PhysScanWithFixings.cpp'
# ============
MAIN_SRCS = \
  GenParamScan.cpp \
  HiggsParamScan.cpp \
  PhysParamScan.cpp \
  testing.cpp \
  Oracle.cpp \
  ModularParamScan.cpp \
  PhysScanWithFixings.cpp

# ============
# Archivos fuente comunes (sin main).
# ============
COMMON_SRCS = \
  ParamUtils.cpp

# Generamos la lista completa de fuentes
SRCS = $(MAIN_SRCS:%=$(SRC_DIR)/%) $(COMMON_SRCS:%=$(SRC_DIR)/%)

# Objetos correspondientes
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))

# ============
# Nombre de los binarios resultantes
# ============
BINARIES = \
  $(BIN_DIR)/GenParamScan \
  $(BIN_DIR)/HiggsParamScan \
  $(BIN_DIR)/PhysParamScan \
  $(BIN_DIR)/testing \
  $(BIN_DIR)/Oracle \
  $(BIN_DIR)/ModularParamScan \
  $(BIN_DIR)/PhysScanWithFixings

# Regla principal
all: $(BIN_DIR) $(OBJ_DIR) $(BINARIES)

# ----
# Reglas de enlace de cada ejecutable
# ----
$(BIN_DIR)/GenParamScan: $(OBJ_DIR)/GenParamScan.o $(OBJ_DIR)/ParamUtils.o
	$(CXX) -fopenmp $^ -o $@ $(LDFLAGS)

$(BIN_DIR)/HiggsParamScan: $(OBJ_DIR)/HiggsParamScan.o $(OBJ_DIR)/ParamUtils.o
	$(CXX) -fopenmp $^ -o $@ $(LDFLAGS)

$(BIN_DIR)/PhysParamScan: $(OBJ_DIR)/PhysParamScan.o $(OBJ_DIR)/ParamUtils.o
	$(CXX) -fopenmp $^ -o $@ $(LDFLAGS)

$(BIN_DIR)/testing: $(OBJ_DIR)/testing.o $(OBJ_DIR)/ParamUtils.o
	$(CXX) -fopenmp $^ -o $@ $(LDFLAGS)

$(BIN_DIR)/Oracle: $(OBJ_DIR)/Oracle.o $(OBJ_DIR)/ParamUtils.o
	$(CXX) -fopenmp $^ -o $@ $(LDFLAGS)

$(BIN_DIR)/ModularParamScan: $(OBJ_DIR)/ModularParamScan.o $(OBJ_DIR)/ParamUtils.o
	$(CXX) -fopenmp $^ -o $@ $(LDFLAGS)

$(BIN_DIR)/PhysScanWithFixings: $(OBJ_DIR)/PhysScanWithFixings.o $(OBJ_DIR)/ParamUtils.o
	$(CXX) -fopenmp $^ -o $@ $(LDFLAGS)

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