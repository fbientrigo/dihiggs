# Unit Tests Configuration
# Configures the binary path at CMake time and processes the test template

# Find the PhysLam1Scan executable (adjust path based on your build structure)
# This assumes the binary is built in the dihiggs build tree
set(PHYSLAM1SCAN_BINARY_PATH "${CMAKE_BINARY_DIR}/dihiggs/app/PhysLam1Scan")

# Configure the test file from the template
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/test_PhysLam1Scan.cpp.in"
    "${CMAKE_CURRENT_BINARY_DIR}/test_PhysLam1Scan.cpp"
    @ONLY
)

# Note: If you're using the Makefile approach, you may need to:
# 1. Run this CMake configuration as part of your build
# 2. Or manually pass the binary path to your build system

message(STATUS "PhysLam1Scan binary path: ${PHYSLAM1SCAN_BINARY_PATH}")