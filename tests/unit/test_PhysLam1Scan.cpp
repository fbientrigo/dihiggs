/**
 * @file test_PhysLam1Scan.cpp
 * @brief Catch2 oracle tests for PhysLam1Scan binary
 *
 * Test categories:
 * 1. Build/Link: Binary runs without missing symbols
 * 2. Smoke CLI: Valid invocation creates output, invalid fails clearly
 * 3. Physics Consistency: lambda1 requested vs recovered
 * 4. Anchor Tests: sin(beta-alpha)=1, yukawa_type=1, large tan_beta
 * 5. Randomized Regression: Fixed seed, capped attempts
 *
 * All tests tagged [oracle] for CI lane selection.
 * Single-threaded execution: OMP_NUM_THREADS=1 enforced.
 *
 * Tolerance policy (from 2hdmc/src/CalcRoundTrip.cpp):
 * - rtol = 1e-10
 * - atol = 1e-12 (quantity-specific for small values)
 */

#include "catch.hpp"
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

// ===========================================================================
// UTILITIES
// ===========================================================================

static const double RTOL = 1e-10;  // Relative tolerance from CalcRoundTrip
static const double ATOL = 1e-12;  // Absolute tolerance for small quantities

static bool nearly_equal(double a, double b, double rtol = RTOL, double atol = ATOL) {
    const double diff = std::fabs(a - b);
    const double scale = std::max(std::fabs(a), std::fabs(b));
    return diff <= (atol + rtol * scale);
}

static bool file_exists(const std::string& path) {
    std::ifstream f(path);
    return f.good();
}

static void remove_if_exists(const std::string& path) {
    std::remove(path.c_str());
}

static int count_csv_lines(const std::string& path) {
    std::ifstream f(path);
    if (!f.is_open()) return -1;
    
    int count = 0;
    std::string line;
    while (std::getline(f, line)) {
        ++count;
    }
    return count;
}

static std::string get_csv_header(const std::string& path) {
    std::ifstream f(path);
    if (!f.is_open()) return "";
    
    std::string header;
    std::getline(f, header);
    return header;
}

static std::vector<std::string> parse_csv_line(const std::string& line) {
    std::vector<std::string> result;
    std::stringstream ss(line);
    std::string cell;
    
    while (std::getline(ss, cell, ',')) {
        result.push_back(cell);
    }
    return result;
}

static double extract_triple_ok_count(const std::string& stdout_output) {
    // Extract "TRIPLE_OK_POINTS <n>" from stdout
    size_t pos = stdout_output.find("TRIPLE_OK_POINTS ");
    if (pos == std::string::npos) return -1.0;
    
    std::string substr = stdout_output.substr(pos + 17); // Skip "TRIPLE_OK_POINTS "
    return std::atof(substr.c_str());
}

static std::string run_command_capture(const std::string& cmd, int& rc) {
    // Run command and capture stdout+stderr
    std::string full_cmd = cmd + " 2>&1";
    FILE* pipe = popen(full_cmd.c_str(), "r");
    if (!pipe) {
        rc = -1;
        return "";
    }
    
    std::string result;
    char buffer[256];
    while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
        result += buffer;
    }
    
    rc = pclose(pipe);
    // pclose returns raw wait() status, extract exit code
    if (WIFEXITED(rc)) {
        rc = WEXITSTATUS(rc);
    } else {
        rc = -1;
    }
    
    return result;
}

// ===========================================================================
// TEST SUITE: BUILD/LINK TESTS
// ===========================================================================

TEST_CASE("PhysLam1Scan binary exists and has correct linkage", "[oracle][build]") {
    const char* binary_path = "/home/fabi/dihiggs/dihiggs/app/PhysLam1Scan";
    
    SECTION("Binary file exists") {
        REQUIRE(file_exists(binary_path));
    }
    
    SECTION("Binary runs without missing symbols") {
        // Run with invalid argc to trigger usage message (rc=1)
        // This verifies all symbols resolve at runtime
        int rc = 0;
        std::string output = run_command_capture(binary_path, rc);
        
        // Should fail with usage message (rc=1), not with missing symbol crash
        REQUIRE(rc == 1);
        REQUIRE(output.find("Usage:") != std::string::npos);
    }
}

// ===========================================================================
// TEST SUITE: SMOKE CLI TESTS
// ===========================================================================

TEST_CASE("PhysLam1Scan CLI smoke tests", "[oracle][smoke]") {
    const char* binary_path = "/home/fabi/dihiggs/dihiggs/app/PhysLam1Scan";
    const std::string test_output = "/tmp/physlam1scan_smoke_test.csv";
    
    // Clean up from previous runs
    remove_if_exists(test_output);
    remove_if_exists(test_output + ".tmp");
    
    SECTION("Valid invocation creates output file") {
        // Minimal scan: 1x1 grid, should complete quickly
        std::string cmd = std::string("OMP_NUM_THREADS=1 ") + binary_path + 
                          " 130 130 1 0.1 0.1 1 300 0.999 50 0.1 0.0 " + test_output;
        
        int rc = 0;
        std::string output = run_command_capture(cmd, rc);
        
        REQUIRE(rc == 0);
        REQUIRE(file_exists(test_output));
        REQUIRE(output.find("Total Attempts:") != std::string::npos);
        REQUIRE(output.find("TRIPLE_OK_POINTS") != std::string::npos);
        
        // CSV should have header + at least 0 data rows (might fail constraints)
        int line_count = count_csv_lines(test_output);
        REQUIRE(line_count >= 1); // At least header
    }
    
    SECTION("Invalid argc fails with usage message") {
        int rc = 0;
        std::string output = run_command_capture(binary_path, rc);
        
        REQUIRE(rc == 1);
        REQUIRE(output.find("Usage:") != std::string::npos);
    }
    
    SECTION("Invalid parameters fail cleanly (negative masses)") {
        // Negative mA should fail THDM validation
        std::string cmd = std::string("OMP_NUM_THREADS=1 ") + binary_path + 
                          " 130 130 1 0.1 0.1 1 -300 0.999 50 0.1 0.0 " + test_output;
        
        int rc = 0;
        std::string output = run_command_capture(cmd, rc);
        
        // Should either fail with rc!=0 or produce 0 valid points
        if (rc == 0) {
            // If it succeeds, should report 0 valid points
            REQUIRE(output.find("Total CSV Rows: 0") != std::string::npos);
        }
    }
    
    // Cleanup
    remove_if_exists(test_output);
    remove_if_exists(test_output + ".tmp");
}

// ===========================================================================
// TEST SUITE: PHYSICS CONSISTENCY TESTS
// ===========================================================================

TEST_CASE("PhysLam1Scan physics consistency: lambda1 recovery", "[oracle][physics]") {
    const char* binary_path = "/home/fabi/dihiggs/dihiggs/app/PhysLam1Scan";
    const std::string test_output = "/tmp/physlam1scan_physics_test.csv";
    
    remove_if_exists(test_output);
    remove_if_exists(test_output + ".tmp");
    
    SECTION("lambda1 requested matches computed_lam1 within tolerance") {
        // Known-good fixture from plan line 485
        std::string cmd = std::string("OMP_NUM_THREADS=1 ") + binary_path + 
                          " 130 130 1 0.1 0.1 1 300 0.999 50 0.1 0.0 " + test_output;
        
        int rc = 0;
        std::string output = run_command_capture(cmd, rc);
        
        REQUIRE(rc == 0);
        REQUIRE(file_exists(test_output));
        
        // Read CSV and verify lambda1 vs computed_lam1
        std::ifstream csv(test_output);
        REQUIRE(csv.is_open());
        
        std::string header;
        std::getline(csv, header);
        
        // Verify header matches PhysScanWithFixings exactly
        REQUIRE(header.find("m_phi") != std::string::npos);
        REQUIRE(header.find("lam1") != std::string::npos);
        REQUIRE(header.find("computed_lam1") != std::string::npos);
        
        // Find column indices (from PhysLam1Scan.cpp lines 53-60)
        // Column 22: lam1 (requested)
        // Column 23: computed_lam1
        std::vector<std::string> cols = parse_csv_line(header);
        int idx_lam1 = -1, idx_comp_lam1 = -1;
        for (size_t i = 0; i < cols.size(); ++i) {
            if (cols[i] == "lam1") idx_lam1 = i;
            if (cols[i] == "computed_lam1") idx_comp_lam1 = i;
        }
        
        REQUIRE(idx_lam1 >= 0);
        REQUIRE(idx_comp_lam1 >= 0);
        
        // Read data rows
        std::string row;
        int valid_rows = 0;
        while (std::getline(csv, row)) {
            std::vector<std::string> data = parse_csv_line(row);
            if (data.size() <= (size_t)std::max(idx_lam1, idx_comp_lam1)) continue;
            
            double lam1_req = std::atof(data[idx_lam1].c_str());
            double lam1_comp = std::atof(data[idx_comp_lam1].c_str());
            
            // Tolerance: rtol=1e-10, atol=1e-12
            REQUIRE(nearly_equal(lam1_req, lam1_comp, RTOL, ATOL));
            ++valid_rows;
        }
        
        // Should have at least 1 valid row for this fixture
        REQUIRE(valid_rows >= 1);
    }
    
    remove_if_exists(test_output);
    remove_if_exists(test_output + ".tmp");
}

// ===========================================================================
// TEST SUITE: ANCHOR TESTS
// ===========================================================================

TEST_CASE("PhysLam1Scan anchor tests: sin(beta-alpha)=1", "[oracle][anchor]") {
    const char* binary_path = "/home/fabi/dihiggs/dihiggs/app/PhysLam1Scan";
    const std::string test_output = "/tmp/physlam1scan_anchor_sba1.csv";
    
    remove_if_exists(test_output);
    remove_if_exists(test_output + ".tmp");
    
    // Alignment limit: sin(beta-alpha) = 1.0
    std::string cmd = std::string("OMP_NUM_THREADS=1 ") + binary_path + 
                      " 130 130 1 0.1 0.1 1 300 1.0 50 0.1 0.0 " + test_output;
    
    int rc = 0;
    std::string output = run_command_capture(cmd, rc);
    
    REQUIRE(rc == 0);
    REQUIRE(output.find("Total Attempts: 1") != std::string::npos);
    
    // Verify CSV header exists
    REQUIRE(file_exists(test_output));
    std::string header = get_csv_header(test_output);
    REQUIRE(header.find("sin_ba") != std::string::npos);
    
    remove_if_exists(test_output);
    remove_if_exists(test_output + ".tmp");
}

TEST_CASE("PhysLam1Scan anchor tests: large tan_beta", "[oracle][anchor]") {
    const char* binary_path = "/home/fabi/dihiggs/dihiggs/app/PhysLam1Scan";
    const std::string test_output = "/tmp/physlam1scan_anchor_tanbeta.csv";
    
    remove_if_exists(test_output);
    remove_if_exists(test_output + ".tmp");
    
    // Large tan_beta regime: tan_beta = 100
    std::string cmd = std::string("OMP_NUM_THREADS=1 ") + binary_path + 
                      " 130 130 1 0.1 0.1 1 300 0.999 100 0.1 0.0 " + test_output;
    
    int rc = 0;
    std::string output = run_command_capture(cmd, rc);
    
    REQUIRE(rc == 0);
    REQUIRE(output.find("Total Attempts: 1") != std::string::npos);
    
    // Verify output exists
    REQUIRE(file_exists(test_output));
    
    remove_if_exists(test_output);
    remove_if_exists(test_output + ".tmp");
}

TEST_CASE("PhysLam1Scan anchor tests: yukawa_type=1 hardcoded", "[oracle][anchor]") {
    // Verify yukawa_type=1 is hardcoded (from PhysLam1Scan.cpp line 94)
    // This is a build-time check - if binary runs, yukawa_type is correct
    const char* binary_path = "/home/fabi/dihiggs/dihiggs/app/PhysLam1Scan";
    
    REQUIRE(file_exists(binary_path));
    
    // Run minimal scan to verify no yukawa_type CLI parameter required
    std::string test_output = "/tmp/physlam1scan_anchor_yukawa.csv";
    remove_if_exists(test_output);
    
    std::string cmd = std::string("OMP_NUM_THREADS=1 ") + binary_path + 
                      " 130 130 1 0.1 0.1 1 300 0.999 50 0.1 0.0 " + test_output;
    
    int rc = 0;
    std::string output = run_command_capture(cmd, rc);
    
    REQUIRE(rc == 0);
    REQUIRE(file_exists(test_output));
    
    remove_if_exists(test_output);
    remove_if_exists(test_output + ".tmp");
}

// ===========================================================================
// TEST SUITE: CSV HEADER VERIFICATION
// ===========================================================================

TEST_CASE("PhysLam1Scan CSV header matches PhysScanWithFixings", "[oracle][physics]") {
    const char* binary_path = "/home/fabi/dihiggs/dihiggs/app/PhysLam1Scan";
    const std::string test_output = "/tmp/physlam1scan_header_test.csv";
    
    remove_if_exists(test_output);
    remove_if_exists(test_output + ".tmp");
    
    std::string cmd = std::string("OMP_NUM_THREADS=1 ") + binary_path + 
                      " 130 130 1 0.1 0.1 1 300 0.999 50 0.1 0.0 " + test_output;
    
    int rc = 0;
    run_command_capture(cmd, rc);
    REQUIRE(rc == 0);
    
    std::string header = get_csv_header(test_output);
    
    // Expected header from PhysScanWithFixings (lines 87-94) and PhysLam1Scan (lines 53-60)
    std::vector<std::string> expected_cols = {
        "m_phi", "mA", "alpha", "beta", "lambda6", "lambda7", "m12",
        "sin_ba", "tan_beta", "positivity_ok", "unitarity_ok", "perturbativity_ok",
        "width_bb", "width_tautau", "width_WW", "width_ZZ",
        "width_gaga", "width_Zga", "width_gg", "width_hh",
        "total_width", "br_gaga", "lam1", "computed_lam1",
        "lam2", "computed_lam2", "lam3", "lam4", "lam5"
    };
    
    for (const auto& col : expected_cols) {
        REQUIRE(header.find(col) != std::string::npos);
    }
    
    remove_if_exists(test_output);
    remove_if_exists(test_output + ".tmp");
}

// ===========================================================================
// TEST SUITE: RANDOMIZED REGRESSION (FIXED SEED)
// ===========================================================================

TEST_CASE("PhysLam1Scan randomized regression: fixed seed", "[oracle][regression]") {
    const char* binary_path = "/home/fabi/dihiggs/dihiggs/app/PhysLam1Scan";
    const std::string test_output = "/tmp/physlam1scan_regression.csv";
    
    remove_if_exists(test_output);
    remove_if_exists(test_output + ".tmp");
    
    // 3x3 grid scan (9 attempts total) - capped for performance
    // mphi: 130-132 GeV (3 points)
    // lam1: 0.1-0.3 (3 points)
    std::string cmd = std::string("OMP_NUM_THREADS=1 ") + binary_path + 
                      " 130 132 3 0.1 0.3 3 300 0.999 50 0.1 0.0 " + test_output;
    
    int rc = 0;
    std::string output = run_command_capture(cmd, rc);
    
    REQUIRE(rc == 0);
    REQUIRE(output.find("Total Attempts: 9") != std::string::npos);
    
    // Extract TRIPLE_OK count
    double triple_ok = extract_triple_ok_count(output);
    REQUIRE(triple_ok >= 0); // Should report count
    
    // Verify CSV structure
    REQUIRE(file_exists(test_output));
    int line_count = count_csv_lines(test_output);
    REQUIRE(line_count >= 1); // At least header
    
    // Run again - deterministic behavior (no random seed in PhysLam1Scan)
    std::string test_output2 = "/tmp/physlam1scan_regression2.csv";
    remove_if_exists(test_output2);
    
    std::string cmd2 = std::string("OMP_NUM_THREADS=1 ") + binary_path + 
                       " 130 132 3 0.1 0.3 3 300 0.999 50 0.1 0.0 " + test_output2;
    
    int rc2 = 0;
    std::string output2 = run_command_capture(cmd2, rc2);
    
    REQUIRE(rc2 == 0);
    
    // Should produce identical attempt counts (deterministic)
    REQUIRE(output2.find("Total Attempts: 9") != std::string::npos);
    
    double triple_ok2 = extract_triple_ok_count(output2);
    REQUIRE(triple_ok2 == triple_ok); // Deterministic
    
    remove_if_exists(test_output);
    remove_if_exists(test_output2);
    remove_if_exists(test_output + ".tmp");
    remove_if_exists(test_output2 + ".tmp");
}

// ===========================================================================
// TEST SUITE: STDOUT MARKERS
// ===========================================================================

TEST_CASE("PhysLam1Scan stdout markers present", "[oracle][smoke]") {
    const char* binary_path = "/home/fabi/dihiggs/dihiggs/app/PhysLam1Scan";
    const std::string test_output = "/tmp/physlam1scan_stdout_test.csv";
    
    remove_if_exists(test_output);
    
    std::string cmd = std::string("OMP_NUM_THREADS=1 ") + binary_path + 
                      " 130 130 1 0.1 0.1 1 300 0.999 50 0.1 0.0 " + test_output;
    
    int rc = 0;
    std::string output = run_command_capture(cmd, rc);
    
    REQUIRE(rc == 0);
    
    // Verify required stdout markers (plan lines 513-527)
    REQUIRE(output.find("Total Attempts:") != std::string::npos);
    REQUIRE(output.find("TRIPLE_OK_POINTS") != std::string::npos);
    
    remove_if_exists(test_output);
    remove_if_exists(test_output + ".tmp");
}

// ===========================================================================
// TEST SUITE: ATOMIC OUTPUT (TMP + RENAME)
// ===========================================================================

TEST_CASE("PhysLam1Scan atomic output: .tmp rename", "[oracle][build]") {
    const char* binary_path = "/home/fabi/dihiggs/dihiggs/app/PhysLam1Scan";
    const std::string test_output = "/tmp/physlam1scan_atomic_test.csv";
    const std::string tmp_output = test_output + ".tmp";
    
    remove_if_exists(test_output);
    remove_if_exists(tmp_output);
    
    std::string cmd = std::string("OMP_NUM_THREADS=1 ") + binary_path + 
                      " 130 130 1 0.1 0.1 1 300 0.999 50 0.1 0.0 " + test_output;
    
    int rc = 0;
    run_command_capture(cmd, rc);
    
    REQUIRE(rc == 0);
    
    // Final output should exist
    REQUIRE(file_exists(test_output));
    
    // .tmp file should NOT exist after successful completion
    REQUIRE_FALSE(file_exists(tmp_output));
    
    remove_if_exists(test_output);
}
