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

# Dependencias (2HDMC)
THDMC_SRC         := ../2hdmc/src
THDMC_LIB         := ../2hdmc/lib

# Flags
CXXFLAGS := -I$(THDMC_SRC) -I$(SRC_DIR) -Iinclude
CXXFLAGS += -fopenmp -std=gnu++17 -Wall -Wextra -Wpedantic -O0 -g
LDFLAGS  := -L$(LIB_DIR) -ldihiggs \
            -L$(THDMC_LIB) -l2HDMC \
            -lgsl -lgslcblas -lm

# Fuentes y Objetos
# PhysLam1Scan is not in MAIN_SRCS: it is the frozen legacy oracle behind
# tests/test_golden_lambda1.py, built by scripts/build_lambda1_characterization.sh
# (or `make app/PhysLam1Scan` for tests/unit).
MAIN_SRCS   := $(SRC_DIR)/Lambda1EvaluatorV2.cpp \
               $(SRC_DIR)/DihiggsPointV2Evaluator.cpp \
               $(SRC_DIR)/Phys_M2BandTracker.cpp

COMMON_SRCS := $(SRC_DIR)/ParamUtils.cpp \
               $(SRC_DIR)/ReplaySafeOutput.cpp \
               $(SRC_DIR)/M2PointEvaluator.cpp \
               $(SRC_DIR)/M2BatchEvaluator.cpp \
               $(SRC_DIR)/M2IntervalDetector.cpp \
               $(SRC_DIR)/M2ContinuationPredictor.cpp \
               $(SRC_DIR)/M2BoundaryRefiner.cpp \
               $(SRC_DIR)/M2FallbackPolicy.cpp

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

$(BIN_DIR)/Lambda1EvaluatorV2: $(OBJ_DIR)/Lambda1EvaluatorV2.o $(OBJ_DIR)/ReplaySafeOutput.o
	@echo "[LD] Linking executable $@ (2HDMC only)"
	@mkdir -p $(BIN_DIR)
	@$(CXX) -fopenmp $^ -o $@ \
	    -L$(THDMC_LIB) -l2HDMC \
	    -lgsl -lgslcblas -lm

$(BIN_DIR)/DihiggsPointV2Evaluator: $(OBJ_DIR)/DihiggsPointV2Evaluator.o $(OBJ_DIR)/ReplaySafeOutput.o
	@echo "[LD] Linking executable $@ (2HDMC only)"
	@mkdir -p $(BIN_DIR)
	@$(CXX) -fopenmp $^ -o $@ \
	    -L$(THDMC_LIB) -l2HDMC \
	    -lgsl -lgslcblas -lm

# REGLA CRÍTICA: Compilación de .cpp a .o
# Se elimina la dependencia de pipe '|' y se añade 'mkdir -p' interno
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@echo "[CXX] Compiling $<"
	@if [ ! -f $< ]; then echo "[ERROR] Source file $< not found"; exit 1; fi
	@mkdir -p $(dir $@)
	@$(CXX) $(CXXFLAGS) -c $< -o $@

directories:
	@mkdir -p $(BIN_DIR) $(OBJ_DIR) $(LIB_DIR) include

clean:
	@echo "[CLEAN] Removing objects and binaries"
	@rm -rf $(OBJ_DIR) $(LIB_DIR)
	@rm -f $(BINARIES) $(BIN_DIR)/PhysLam1Scan

cleanall: clean
	@rm -rf data/*
