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
CXXFLAGS := -I$(THDMC_SRC) -I$(HIGGSTOOLS_INC) -I$(SRC_DIR) -Iinclude
CXXFLAGS += -fopenmp -std=gnu++17 -Wall -Wextra -Wpedantic -O0 -g
LDFLAGS  := -L$(LIB_DIR) -ldihiggs \
            -L$(THDMC_LIB) -l2HDMC \
            -L$(HIGGSTOOLS_LIB) -lHiggsTools \
            -lgsl -lgslcblas -lm \
            -Wl,-rpath,$(HIGGSTOOLS_LIB)

# Fuentes y Objetos
MAIN_SRCS   := $(wildcard $(SRC_DIR)/PhysParamScan.cpp $(SRC_DIR)/Oracle.cpp \
               $(SRC_DIR)/HBHS_testing.cpp $(SRC_DIR)/PostProcessHBHS.cpp \
               $(SRC_DIR)/PhysScanWithFixings.cpp $(SRC_DIR)/PhysLam1Scan.cpp \
			   $(SRC_DIR)/ColabRecomputeValidation.cpp \
               $(SRC_DIR)/Phys_M2BandTracker.cpp \
               $(SRC_DIR)/GenScanWithFixings.cpp)

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

.PHONY: all clean cleanall directories ci-calcphys-gate

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

# GenScanWithFixings is built standalone, bypassing libdihiggs.a: COMMON_SRCS
# (M2PointEvaluator.cpp etc.) currently fails to compile because
# M2PointEvaluator.hpp does not exist (pre-existing breakage, unrelated to
# this binary). Once that's fixed, this could be folded back into the
# default $(BIN_DIR)/% pattern rule.
# Note: deliberately does NOT link -lHiggsTools (HiggsTools is not built in
# this checkout, and GenScanPointEvaluator/GenScanWithFixings only use
# 2HDMC's THDM/Constraints/DecayTable, never HiggsTools).
$(BIN_DIR)/GenScanWithFixings: $(OBJ_DIR)/GenScanWithFixings.o $(OBJ_DIR)/GenScanPointEvaluator.o $(OBJ_DIR)/ReplaySafeOutput.o
	@echo "[LD] Linking executable $@ (standalone, bypasses libdihiggs.a)"
	@mkdir -p $(BIN_DIR)
	@$(CXX) -fopenmp $^ -o $@ \
	    -L$(THDMC_LIB) -l2HDMC \
	    -lgsl -lgslcblas -lm

directories:
	@mkdir -p $(BIN_DIR) $(OBJ_DIR) $(LIB_DIR) include

clean:
	@echo "[CLEAN] Removing objects and binaries"
	@rm -rf $(OBJ_DIR) $(LIB_DIR)
	@rm -f $(BINARIES)

cleanall: clean
	@rm -rf data/*

GATE_THREADS_HI ?= 4
GATE_STEPS ?= 25
GATE_OUT_DIR ?= /tmp/calcphys_ci_gate

# Paths are anchored to the checkout (MAKEFILE_DIR / REPO_ROOT) so they resolve
# correctly regardless of the cwd used to invoke the target (e.g.
# `make -C dihiggs ci-calcphys-gate`). The scan and CalcPhys binaries are built
# from committed sources at these locations.
GATE_SCRIPT   ?= $(REPO_ROOT)/scripts/out/calcphys_ci_gate.py
GATE_SCAN_EXE ?= $(MAKEFILE_DIR)/app/PhysScanWithFixings
GATE_CALCPHYS ?= $(REPO_ROOT)/2hdmc/CalcPhys

ci-calcphys-gate: $(BIN_DIR)/PhysScanWithFixings
	@echo "[CI] Running CalcPhys parity + thread gate (1 thread vs >=4 threads)"
	@if [ ! -f "$(GATE_SCRIPT)" ]; then \
		echo "[ERROR] CalcPhys CI gate script not found: $(GATE_SCRIPT)"; \
		echo "[ERROR] This harness is not committed to the repository, so the"; \
		echo "[ERROR] target cannot run. Provide it via GATE_SCRIPT=/path/to/calcphys_ci_gate.py"; \
		echo "[ERROR] or restore it into scripts/out/ before invoking this target."; \
		exit 1; \
	fi
	@python3 "$(GATE_SCRIPT)" --repo-root "$(REPO_ROOT)" --scan-exe "$(GATE_SCAN_EXE)" --calcphys "$(GATE_CALCPHYS)" --threads-hi $(GATE_THREADS_HI) --steps $(GATE_STEPS) --out-dir $(GATE_OUT_DIR)
