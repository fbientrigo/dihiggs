# =============================================================================
# DIHIGGS UNIT TESTS MAKEFILE
# =============================================================================

CXX      := g++
REPO_ROOT := ../..

# Absolute project root (resolved at make time so the binary path baked into
# the test binary matches the container path at test-run time)
REPO_ROOT_ABS := $(abspath $(REPO_ROOT))

# Include paths (Core Src + Dependencies headers)
INCLUDES := -I$(REPO_ROOT)/dihiggs/src \
            -I$(REPO_ROOT)/2hdmc/src \
            -I$(REPO_ROOT)/higgstools/installation/include

# Compilation Flags (Must match Core flags essentially)
CXXFLAGS := -std=c++17 -Wall -Wextra -O0 -g $(INCLUDES)

# Bake the absolute path to PhysLam1Scan into the test binary so that oracle
# tests can locate the executable regardless of the working directory.
CXXFLAGS += -DPHYSLAM1SCAN_BINARY_PATH=\"$(REPO_ROOT_ABS)/dihiggs/app/PhysLam1Scan\"

# Library Paths
LIB_CORE      := $(REPO_ROOT)/dihiggs/lib
LIB_2HDMC     := $(REPO_ROOT)/2hdmc/lib
LIB_HIGGSTOOLS:= $(REPO_ROOT)/higgstools/installation/lib

# Link against the Core Static Library
LDFLAGS  := -L$(LIB_CORE) -ldihiggs \
            -L$(LIB_2HDMC) -l2HDMC \
            -L$(LIB_HIGGSTOOLS) -lHiggsTools \
            -lgsl -lgslcblas -lm \
            -Wl,-rpath,$(LIB_HIGGSTOOLS)

# Sources
TEST_SRCS := test_main.cpp test_ParamUtils.cpp test_PhysLam1Scan.cpp
TARGET    := dihiggs_unit_tests

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(TEST_SRCS)
	@echo "[TEST-BUILD] Compiling Unit Tests..."
	@# Note: We depend on the main lib being built first. 
	@# Ideally, the CI script calls 'make' in root first.
	@$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)
	@echo "[TEST-BUILD] Success: ./$(TARGET)"

run: $(TARGET)
	@echo "[TEST-RUN] Executing..."
	@./$(TARGET) -s

clean:
	@rm -f $(TARGET)