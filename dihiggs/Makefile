# =============================================================================
# DIHIGGS CORE MAKEFILE - VERSION ROBUSTA CI/CD
# =============================================================================

CXX      := g++
SHELL    := /bin/bash

MAKEFILE_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
REPO_ROOT    := $(abspath $(MAKEFILE_DIR)/..)

# Directorios
SRC_DIR      := src
OBJ_DIR      := obj
BIN_DIR      := app
LIB_DIR      := lib

# Dependencias (HiggsTools & 2HDMC)
HIGGSTOOLS_PREFIX ?= $(REPO_ROOT)/higgstools/installation
HIGGSTOOLS_INC    := $(HIGGSTOOLS_PREFIX)/include
HIGGSTOOLS_LIB    := $(HIGGSTOOLS_PREFIX)/lib
THDMC_SRC         := ../2hdmc/src
THDMC_LIB         := ../2hdmc/lib

# Flags
CXXFLAGS := -I$(THDMC_SRC) -I$(HIGGSTOOLS_INC) -I$(SRC_DIR)
CXXFLAGS += -fopenmp -std=gnu++17 -Wall -Wextra -Wpedantic -O0 -g
LDFLAGS  := -L$(LIB_DIR) -ldihiggs \
            -L$(THDMC_LIB) -l2HDMC \
            -L$(HIGGSTOOLS_LIB) -lHiggsTools \
            -lgsl -lgslcblas -lm \
            -Wl,-rpath,$(HIGGSTOOLS_LIB)

# Fuentes y Objetos
MAIN_SRCS   := $(wildcard $(SRC_DIR)/PhysParamScan.cpp $(SRC_DIR)/Oracle.cpp \
               $(SRC_DIR)/HBHS_testing.cpp $(SRC_DIR)/PostProcessHBHS.cpp \
               $(SRC_DIR)/PhysScanWithFixings.cpp $(SRC_DIR)/PhysLam1Scan.cpp)
COMMON_SRCS := $(SRC_DIR)/ParamUtils.cpp

# Generación de rutas de objetos
MAIN_OBJS   := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(MAIN_SRCS))
COMMON_OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(COMMON_SRCS))
BINARIES    := $(patsubst $(SRC_DIR)/%.cpp,$(BIN_DIR)/%,$(MAIN_SRCS))
TARGET_LIB  := $(LIB_DIR)/libdihiggs.a

# =============================================================================
# REGLAS DE CONSTRUCCIÓN
# =============================================================================

.PHONY: all clean cleanall directories

# Forzamos que 'directories' ocurra antes de cualquier otra cosa
all: directories
	@$(MAKE) $(TARGET_LIB)
	@$(MAKE) $(BINARIES)
	@echo "[INFO] Build complete."

# Regla de la Librería Estática
$(TARGET_LIB): $(COMMON_OBJS)
	@echo "[AR] Creating static library $@"
	@mkdir -p $(LIB_DIR)
	@ar rcs $@ $^

# Regla de enlace de Ejecutables
$(BIN_DIR)/%: $(OBJ_DIR)/%.o $(TARGET_LIB)
	@echo "[LD] Linking executable $@"
	@mkdir -p $(BIN_DIR)
	@$(CXX) -fopenmp $< -o $@ $(LDFLAGS)

# REGLA CRÍTICA: Compilación de .cpp a .o
# Se elimina la dependencia de pipe '|' y se añade 'mkdir -p' interno
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@echo "[CXX] Compiling $<"
	@if [ ! -f $< ]; then echo "[ERROR] Source file $< not found"; exit 1; fi
	@mkdir -p $(dir $@)
	@$(CXX) $(CXXFLAGS) -c $< -o $@

directories:
	@mkdir -p $(BIN_DIR) $(OBJ_DIR) $(LIB_DIR)

clean:
	@echo "[CLEAN] Removing objects and binaries"
	@rm -rf $(OBJ_DIR) $(BIN_DIR) $(LIB_DIR)

cleanall: clean
	@rm -rf data/*