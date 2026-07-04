#!/bin/bash
# ============================================================================
# Pipeline Runner for Parquet Analysis
# ============================================================================
# Usage:
#   ./run_pipeline.sh all                  # Run full workflow
#   ./run_pipeline.sh eda                  # Phase 1 only
#   ./run_pipeline.sh ctau                 # Phase 2 only
#   ./run_pipeline.sh mphi                 # Phase 3 only
#   ./run_pipeline.sh parallel             # Phase 4 only
# ============================================================================

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Configuration
PARQUET_FILE="./temp_subspace_l6_tb_high.parquet"
WORK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Optional: point HIGGS_ENV_ACTIVATE at your virtualenv's activate script
ENV_PATH="${HIGGS_ENV_ACTIVATE:-}"

# ============================================================================
# Helper Functions
# ============================================================================
print_header() {
    echo -e "\n${BLUE}════════════════════════════════════════════════════════════${NC}"
    echo -e "${BLUE}  $1${NC}"
    echo -e "${BLUE}════════════════════════════════════════════════════════════${NC}\n"
}

print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

print_error() {
    echo -e "${RED}✗ $1${NC}"
}

print_info() {
    echo -e "${YELLOW}ℹ $1${NC}"
}

check_environment() {
    print_header "Step 1: Checking Environment"
    
    cd "$WORK_DIR"
    print_success "Working directory: $WORK_DIR"

    if [ -n "$ENV_PATH" ]; then
        if [ ! -f "$ENV_PATH" ]; then
            print_error "Python environment not found at: $ENV_PATH"
            exit 1
        fi
        source "$ENV_PATH"
        print_success "Environment activated"
    else
        print_info "HIGGS_ENV_ACTIVATE not set; using current Python environment"
    fi

    if [ ! -f "$PARQUET_FILE" ]; then
        print_error "Parquet file not found: $PARQUET_FILE"
        exit 1
    fi
    print_success "Parquet file exists: $PARQUET_FILE"
}

verify_dependencies() {
    print_header "Step 2: Verifying Dependencies"
    
    python -c "import polars; print(f'Polars: {polars.__version__}')" && print_success "Polars available" || {
        print_error "Polars not installed"
        echo "Install with: pip install polars --upgrade"
        exit 1
    }
    
    python -c "import pandas; print(f'Pandas: {pandas.__version__}')" && print_success "Pandas available" || {
        print_error "Pandas not installed"
        exit 1
    }
    
    python -c "import matplotlib; print(f'Matplotlib: {matplotlib.__version__}')" && print_success "Matplotlib available" || {
        print_error "Matplotlib not installed"
        exit 1
    }
    
    python -c "import numpy; print(f'NumPy: {numpy.__version__}')" && print_success "NumPy available" || {
        print_error "NumPy not installed"
        exit 1
    }
}

verify_parquet() {
    print_header "Step 3: Verifying Parquet Structure"
    
    python << 'EOF'
import polars as pl
lf = pl.scan_parquet('./temp_subspace_l6_tb_high.parquet')
schema = lf.collect_schema()
row_count = lf.select(pl.len()).collect().item()

required = ['mA', 'lambda6', 'tan_beta', 'br_gaga', 'total_width']
missing = [c for c in required if c not in schema.names()]

print(f"Total rows: {row_count:,}")
print(f"Columns: {len(schema.names())}")
print(f"First 5: {schema.names()[:5]}")

if missing:
    print(f"⚠ Missing columns: {missing}")
else:
    print("✓ All required columns present")
EOF
}

run_eda() {
    print_header "Phase 1: Exploratory Data Analysis (EDA)"
    print_info "Opening Jupyter notebook: EDA_subSpace.ipynb"
    print_info "Run cells manually to generate statistics and phase-space plots"
    
    jupyter notebook EDA_subSpace.ipynb
}

run_ctau() {
    print_header "Phase 2: Lifetime vs Branching Ratio (ctau vs BR)"
    print_info "Generating plots in: ctau_plots/"
    print_info "This may take 5-15 minutes..."
    
    python ctau_vs_br_multislice_patched.py \
        --input ./temp_subspace_l6_tb_high.parquet \
        --color-by both \
        --br-col br_gaga
    
    print_success "ctau plots generated!"
    
    if [ -f "ctau_plots/ctau_summary.json" ]; then
        print_info "Summary: ctau_plots/ctau_summary.json"
    fi
}

run_mphi() {
    print_header "Phase 3: Mass vs Branching Ratio (m_phi vs BR)"
    print_info "Generating publication-quality plots in: paper_plots_mphi_br/"
    print_info "This may take 5-15 minutes..."
    
    python paper_like_mphi_vs_br_patched.py \
        --input ./temp_subspace_l6_tb_high.parquet \
        --mphi-col m_phi \
        --br-col br_gaga
    
    print_success "m_phi plots generated!"
    
    if [ -f "paper_plots_mphi_br/summary.json" ]; then
        print_info "Summary: paper_plots_mphi_br/summary.json"
    fi
}

run_parallel() {
    print_header "Phase 4: Parallel Coordinates (Global + ROI)"
    print_info "Generating plots in: parallel_plots/"
    print_info "ROI defaults: ctau < q20 and br_gaga > q80"

    python parallel_coordinates_roi.py \
        --input ./temp_subspace_l6_tb_high.parquet \
        --background-alpha 0.03 \
        --roi-alpha 0.8

    print_success "Parallel coordinates plots generated!"

    if [ -f "parallel_plots/parallel_summary.json" ]; then
        print_info "Summary: parallel_plots/parallel_summary.json"
    fi
}

show_summary() {
    print_header "📊 Analysis Summary"
    
    echo "Output directories created:"
    
    if [ -d "ctau_plots" ]; then
        count=$(ls -1 ctau_plots/*.png 2>/dev/null | wc -l)
        echo -e "  ${GREEN}✓ ctau_plots/${NC}           ($count PNG files)"
    fi
    
    if [ -d "paper_plots_mphi_br" ]; then
        count=$(ls -1 paper_plots_mphi_br/*.png 2>/dev/null | wc -l)
        echo -e "  ${GREEN}✓ paper_plots_mphi_br/${NC}  ($count PNG files)"
    fi
    
    if [ -d "subspace_comparisons_logs" ]; then
        count=$(ls -1 subspace_comparisons_logs/*.png 2>/dev/null | wc -l)
        echo -e "  ${GREEN}✓ subspace_comparisons_logs/${NC} ($count PNG files)"
    fi

    if [ -d "parallel_plots" ]; then
        count=$(ls -1 parallel_plots/*.png 2>/dev/null | wc -l)
        echo -e "  ${GREEN}✓ parallel_plots/${NC}      ($count PNG files)"
    fi
    
    echo ""
    print_success "All analyses complete!"
}

show_usage() {
    cat << 'EOF'
Usage: ./run_pipeline.sh [STAGE]

Stages:
    all         Run full workflow (all phases)
  check       Environment & dependency check only
  eda         Phase 1: Exploratory data analysis (Jupyter notebook)
  ctau        Phase 2: ctau vs BR plots
  mphi        Phase 3: m_phi vs BR plots
  parallel    Phase 4: Parallel coordinates (global + ROI)
  help        Show this message

Examples:
  ./run_pipeline.sh all           # Full pipeline
  ./run_pipeline.sh ctau          # Generate ctau plots only
  ./run_pipeline.sh check         # Verify setup

EOF
}

# ============================================================================
# Main Execution
# ============================================================================

if [ "$#" -eq 0 ] || [ "$1" == "help" ]; then
    show_usage
    exit 0
fi

STAGE="${1,,}"  # Convert to lowercase

case "$STAGE" in
    check)
        check_environment
        verify_dependencies
        verify_parquet
        print_success "All checks passed!"
        ;;
    eda)
        check_environment
        verify_dependencies
        verify_parquet
        run_eda
        ;;
    ctau)
        check_environment
        verify_dependencies
        verify_parquet
        run_ctau
        show_summary
        ;;
    mphi)
        check_environment
        verify_dependencies
        verify_parquet
        run_mphi
        show_summary
        ;;
    parallel)
        check_environment
        verify_dependencies
        verify_parquet
        run_parallel
        show_summary
        ;;
    all)
        check_environment
        verify_dependencies
        verify_parquet
        
        read -p "Start full pipeline? [y/N] " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            # Note: EDA is skipped in automated run; user should run manually if needed
            print_info "Skipping EDA phase (notebook). Run manually: jupyter notebook EDA_subSpace.ipynb"
            run_ctau
            run_mphi
            run_parallel
            show_summary
        else
            print_info "Pipeline cancelled"
        fi
        ;;
    *)
        print_error "Unknown stage: $STAGE"
        show_usage
        exit 1
        ;;
esac

echo ""
