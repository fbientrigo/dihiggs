#!/bin/bash
# Quick launcher for parquet analysis pipeline documentation & tasks
# Usage: ./launch.sh [doc|run|test|help]

DOCS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VIEWER="less"  # or 'more', 'cat -P', etc.

show_menu() {
    cat << 'EOF'

╔════════════════════════════════════════════════════════════════════╗
║              📊 Parquet Analysis Pipeline Launcher               ║
╚════════════════════════════════════════════════════════════════════╝

📚 DOCUMENTATION:
  1. readme        - Main documentation index (START HERE)
  2. quick-ref     - Quick reference card (fastest lookup)
  3. quick-run     - Complete run guide (detailed walkthrough)
  4. schema        - Parquet schema reference (understand data)

▶ RUN ANALYSIS:
  5. test          - Quick test run (2-3 minutes)
  6. run-all       - Full pipeline (30-40 minutes)
  7. run-ctau      - Just ctau plots (10 minutes)
  8. run-mphi      - Just m_phi plots (10 minutes)
  9. run-compare   - Just comparison plots (5 minutes)

⚙  UTILITIES:
  10. check        - Verify environment
  11. inspect      - Inspect parquet structure
  12. help         - Show this menu

EOF
}

show_menu() {
    cat << 'EOF'

╔════════════════════════════════════════════════════════════════════╗
║              📊 Parquet Analysis Pipeline Launcher               ║
╚════════════════════════════════════════════════════════════════════╝

🎯 Quick Commands:
  readme              View main documentation
  quick-ref           Quick reference card  
  quick-run           Complete guide
  schema              Data schema reference
  test                Run quick test (2 min)
  run-all             Run full pipeline (35 min)
  check               Verify setup
  help                Show full menu

EOF
}

case "${1,,}" in
    readme|index)
        echo "📖 Opening README.md..."
        cat README.md | $VIEWER
        ;;
    quick-ref|quickref|qr)
        echo "📋 Opening QUICK_REFERENCE.md..."
        cat QUICK_REFERENCE.md | $VIEWER
        ;;
    quick-run|quickrun|run-guide)
        echo "📘 Opening QUICK_RUN_GUIDE.md..."
        cat QUICK_RUN_GUIDE.md | $VIEWER
        ;;
    schema|parquet)
        echo "🔍 Opening PARQUET_SCHEMA.md..."
        cat PARQUET_SCHEMA.md | $VIEWER
        ;;
    test|quick-test)
        echo "⚡ Running quick test (max-combinations=3)..."
        python ctau_vs_br_multislice_patched.py \
            --input ./temp_subspace_l6_tb_high.parquet \
            --color-by lambda6 \
            --max-combinations 3
        ;;
    run-all|all)
        echo "🚀 Running full pipeline..."
        ./run_pipeline.sh all
        ;;
    run-ctau|ctau)
        echo "📈 Running ctau plots..."
        ./run_pipeline.sh ctau
        ;;
    run-mphi|mphi)
        echo "📊 Running m_phi plots..."
        ./run_pipeline.sh mphi
        ;;
    run-compare|compare)
        echo "🔄 Running comparison plots..."
        ./run_pipeline.sh compare
        ;;
    check|verify)
        echo "✓ Checking environment..."
        ./run_pipeline.sh check
        ;;
    inspect)
        echo "🔍 Inspecting parquet structure..."
        python << 'INSPECT'
import polars as pl
lf = pl.scan_parquet('./temp_subspace_l6_tb_high.parquet')
schema = lf.collect_schema()
count = lf.select(pl.len()).collect().item()
print(f"Rows: {count:,}")
print(f"Columns: {len(schema.names())}")
print("\nColumn listing:")
for i, col in enumerate(schema.names(), 1):
    print(f"  {i:2d}. {col:<30s} {str(schema[col])}")
INSPECT
        ;;
    help|menu|"")
        show_menu
        ;;
    *)
        echo "❌ Unknown command: $1"
        show_menu
        exit 1
        ;;
esac

