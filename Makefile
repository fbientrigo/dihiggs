# Nombre del compilador
CXX = g++

# Flags de compilación (añade aquí -O2, -Wall, etc. según desees)
CXXFLAGS = -I../2hdmc/src -fopenmp -std=c++11

# Flags de enlazado. Ajusta si necesitas más libs.
LDFLAGS = -L../2hdmc/lib -l2HDMC -lgsl -lgslcblas -lm

# Directorios
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = app

# Archivos fuente principales (tus escaneos)
# Puedes enumerar manualmente o usar wildcards.
SRCS = $(SRC_DIR)/ParamScanGen.cpp \
       $(SRC_DIR)/ParamScanGen_pll.cpp \
       $(SRC_DIR)/ParamScanHiggs.cpp \
       $(SRC_DIR)/ParamUtils.cpp

# Genera la lista de .o sustituyendo src/ por obj/
OBJS = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRCS))

# Aquí defines los nombres de los ejecutables que quieres generar.
# Puedes hacer uno por cada .cpp principal de tu "main".
# Por ejemplo, si ParamScanGen.cpp y ParamScanGen_pll.cpp tienen su "main",
# generamos 2 binarios. ParamScanHiggs.cpp generaría un tercero, etc.
# A continuación, un ejemplo con 3 ejecutables:
BINARIES = $(BIN_DIR)/ParamScanGen \
           $(BIN_DIR)/ParamScanGen_pll \
           $(BIN_DIR)/ParamScanHiggs

# Regla "all" para compilar todo
all: $(BIN_DIR) $(OBJ_DIR) $(BINARIES)

# Regla genérica para compilar cada ejecutable enlazando los .o que necesita
# Aquí hacemos un pequeño truco: cada ejecutable se asocia con su fuente .cpp principal.
# De este modo, compila con todos los .o (incluidos ParamUtils.o, etc.).
$(BIN_DIR)/ParamScanGen: $(OBJ_DIR)/ParamScanGen.o $(OBJ_DIR)/ParamUtils.o
	$(CXX) $^ -o $@ $(LDFLAGS)

$(BIN_DIR)/ParamScanGen_pll: $(OBJ_DIR)/ParamScanGen_pll.o $(OBJ_DIR)/ParamUtils.o
	$(CXX) $^ -o $@ $(LDFLAGS)

$(BIN_DIR)/ParamScanHiggs: $(OBJ_DIR)/ParamScanHiggs.o $(OBJ_DIR)/ParamUtils.o
	$(CXX) $^ -o $@ $(LDFLAGS)

# Regla genérica para compilar .cpp -> .o
# Crea el directorio obj/ si no existe y compila.
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Crea directorios si no existen
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Limpieza
clean:
	rm -rf $(OBJ_DIR)/*.o $(BINARIES)

# "make cleanall" si también quisieras borrar app/ y data/ (opcional)
cleanall: clean
	rm -rf $(BIN_DIR)/*
	rm -rf data/*

.PHONY: all clean cleanall
