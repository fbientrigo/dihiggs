# =============================================================================
# DIHIGGS CORE MAKEFILE
# =============================================================================

# Compiler and Shell setup
CXX      := g++
SHELL    := /bin/bash

# Directory layout
MAKEFILE_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
REPO_ROOT    := $(abspath $(MAKEFILE_DIR)/..)
SRC_DIR      := src
OBJ_DIR      := obj
BIN_DIR      := app
LIB_DIR      := lib

# Dependencies: HiggsTools
HIGGSTOOLS_PREFIX ?= $(REPO_ROOT)/higgstools/installation
HIGGSTOOLS_INC    := $(HIGGSTOOLS_PREFIX)/include
HIGGSTOOLS_LIB    := $(HIGGSTOOLS_PREFIX)/lib

# Dependencies: 2HDMC
THDMC_SRC := ../2hdmc/src
THDMC_LIB := ../2hdmc/lib

# Compilation Flags
# Use -O0 -g for debugging/testing, -O3 for production
CXXFLAGS := -I$(THDMC_SRC) -I$(HIGGSTOOLS_INC) -I$(SRC_DIR)
CXXFLAGS += -fopenmp -std=gnu++17 -Wall -Wextra -Wpedantic -O0 -g

# Linker Flags (Order matters: Local libs -> 2HDMC -> HiggsTools -> System)
LDFLAGS  := -L$(LIB_DIR) -ldihiggs \
            -L$(THDMC_LIB) -l2HDMC \
            -L$(HIGGSTOOLS_LIB) -lHiggsTools \
            -lgsl -lgslcblas -lm \
            -Wl,-rpath,$(HIGGSTOOLS_LIB)

# =============================================================================
# SOURCES & OBJECTS
# =============================================================================

# Archivos con 'main' (ahora con el prefijo del directorio)
MAIN_SRCS := \
  $(SRC_DIR)/PhysParamScan.cpp \
  $(SRC_DIR)/Oracle.cpp \
  $(SRC_DIR)/HBHS_testing.cpp \
  $(SRC_DIR)/PostProcessHBHS.cpp \
  $(SRC_DIR)/PhysScanWithFixings.cpp

# Archivos comunes para la librería (ahora con el prefijo del directorio)
COMMON_SRCS := \
  $(SRC_DIR)/ParamUtils.cpp

# Ajustamos la generación de OBJS para que mantengan la estructura de nombres
MAIN_OBJS   := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(MAIN_SRCS))
COMMON_OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(COMMON_SRCS))

# Los binarios resultantes
BINARIES    := $(patsubst $(SRC_DIR)/%.cpp,$(BIN_DIR)/%,$(MAIN_SRCS))

# The Main Static Library
TARGET_LIB  := $(LIB_DIR)/libdihiggs.a

# =============================================================================
# RULES
# =============================================================================

.PHONY: all clean cleanall info directories

all: directories $(TARGET_LIB) $(BINARIES)
	@echo "[INFO] Build complete."

# --- 1. Compile Static Library ---
$(TARGET_LIB): $(COMMON_OBJS)
	@echo "[AR] Creating static library $@"
	@ar rcs $@ $^

# --- 2. Link Executables ---
# Depends on the .o file AND the static library
$(BIN_DIR)/%: $(OBJ_DIR)/%.o $(TARGET_LIB)
	@echo "[LD] Linking executable $@"
	@$(CXX) -fopenmp $< -o $@ $(LDFLAGS)

# --- 3. Compile Source Files ---
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	@echo "[CXX] Compiling $<"
	@$(CXX) $(CXXFLAGS) -c $< -o $@

# --- Helper Rules ---
directories:
	@mkdir -p $(BIN_DIR) $(OBJ_DIR) $(LIB_DIR)

clean:
	@echo "[CLEAN] Removing objects and binaries"
	@rm -rf $(OBJ_DIR)/*.o $(BINARIES) $(TARGET_LIB)

cleanall: clean
	@echo "[CLEAN] Removing all generated data"
	@rm -rf $(BIN_DIR)/*
	@rm -rf data/*