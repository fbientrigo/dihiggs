#!/usr/bin/env python3
"""
Aggregate and visualize physics results from completed autoresearch campaigns.

Usage:
    python aggregate_results.py /path/to/campaign/output [--output summary.csv] [--plots]
    
Example:
    python aggregate_results.py runs/lam1-tb1000-mphi130-debug-no-convergence --plots
"""

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd


def find_results_files(campaign_dir: Path) -> List[Path]:
    """Find all results.csv files in campaign checkpoint directories."""
    checkpoint_dir = campaign_dir / "checkpoints"
    if not checkpoint_dir.exists():
        print(f"ERROR: No checkpoints directory found at {checkpoint_dir}", file=sys.stderr)
        return []
    
    results_files = list(checkpoint_dir.glob("**/results.csv"))
    return sorted(results_files)


def aggregate_results(results_files: List[Path]) -> pd.DataFrame:
    """Load and concatenate all results.csv files into a single DataFrame."""
    if not results_files:
        return pd.DataFrame()
    
    dfs = []
    for i, filepath in enumerate(results_files):
        try:
            df = pd.read_csv(filepath)
            # Add metadata columns
            parts = filepath.parts
            iter_idx = next((i for i, p in enumerate(parts) if p.startswith("iter_")), None)
            run_idx = next((i for i, p in enumerate(parts) if p.startswith("run_")), None)
            
            if iter_idx is not None:
                df['iteration'] = parts[iter_idx]
            if run_idx is not None:
                df['run_id'] = parts[run_idx]
            
            df['file_path'] = str(filepath)
            dfs.append(df)
        except Exception as e:
            print(f"WARNING: Failed to read {filepath}: {e}", file=sys.stderr)
    
    if not dfs:
        return pd.DataFrame()
    
    return pd.concat(dfs, ignore_index=True)


def summarize_results(df: pd.DataFrame) -> Dict:
    """Generate summary statistics from aggregated results."""
    if df.empty:
        return {}
    
    summary = {
        'total_evaluations': len(df),
        'unique_iterations': df['iteration'].nunique() if 'iteration' in df.columns else 0,
    }
    
    # Constraint satisfaction
    if 'positivity_ok' in df.columns:
        summary['positivity_pass_rate'] = df['positivity_ok'].mean()
    if 'unitarity_ok' in df.columns:
        summary['unitarity_pass_rate'] = df['unitarity_ok'].mean()
    if 'perturbativity_ok' in df.columns:
        summary['perturbativity_pass_rate'] = df['perturbativity_ok'].mean()
    
    # All constraints passed
    constraint_cols = ['positivity_ok', 'unitarity_ok', 'perturbativity_ok']
    if all(c in df.columns for c in constraint_cols):
        df['all_constraints_ok'] = df[constraint_cols].all(axis=1).astype(int)
        summary['all_constraints_pass_rate'] = df['all_constraints_ok'].mean()
        summary['viable_points'] = int(df['all_constraints_ok'].sum())
    
    # Parameter ranges
    param_cols = ['m_phi', 'lambda6', 'lambda7', 'tan_beta', 'lam1', 'lam2', 'lam3', 'lam4', 'lam5']
    for col in param_cols:
        if col in df.columns:
            summary[f'{col}_min'] = df[col].min()
            summary[f'{col}_max'] = df[col].max()
            summary[f'{col}_mean'] = df[col].mean()
    
    # Branching ratio stats
    if 'br_gaga' in df.columns:
        summary['br_gaga_mean'] = df['br_gaga'].mean()
        summary['br_gaga_max'] = df['br_gaga'].max()
    
    return summary


def print_summary(summary: Dict):
    """Pretty-print summary statistics."""
    print("\n" + "="*60)
    print("CAMPAIGN RESULTS SUMMARY")
    print("="*60)
    
    print(f"\n📊 Evaluation Statistics:")
    print(f"  Total evaluations: {summary.get('total_evaluations', 0):,}")
    print(f"  Unique iterations: {summary.get('unique_iterations', 0)}")
    print(f"  Viable points (all constraints): {summary.get('viable_points', 0):,}")
    
    print(f"\n✓ Constraint Pass Rates:")
    print(f"  Positivity:      {summary.get('positivity_pass_rate', 0)*100:.1f}%")
    print(f"  Unitarity:       {summary.get('unitarity_pass_rate', 0)*100:.1f}%")
    print(f"  Perturbativity:  {summary.get('perturbativity_pass_rate', 0)*100:.1f}%")
    print(f"  ALL constraints: {summary.get('all_constraints_pass_rate', 0)*100:.1f}%")
    
    print(f"\n📈 Parameter Coverage:")
    param_names = ['lam1', 'm_phi', 'tan_beta', 'lambda6', 'lambda7']
    for p in param_names:
        if f'{p}_min' in summary:
            print(f"  {p:12s}: [{summary[f'{p}_min']:9.3f}, {summary[f'{p}_max']:9.3f}]  (mean: {summary[f'{p}_mean']:9.3f})")
    
    if 'br_gaga' in summary:
        print(f"\n🔬 Physics Results:")
        print(f"  BR(h→γγ) mean: {summary.get('br_gaga_mean', 0):.6f}")
        print(f"  BR(h→γγ) max:  {summary.get('br_gaga_max', 0):.6f}")
    
    print("\n" + "="*60 + "\n")


def plot_results(df: pd.DataFrame, campaign_dir: Path):
    """Generate visualization plots of exploration results."""
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        print("WARNING: matplotlib/seaborn not installed. Skipping plots.", file=sys.stderr)
        return
    
    if df.empty:
        print("WARNING: No data to plot.", file=sys.stderr)
        return
    
    output_dir = campaign_dir / "analysis_plots"
    output_dir.mkdir(exist_ok=True)
    
    # Set style
    sns.set_style("whitegrid")
    
    # 1. Lambda1 distribution over iterations
    if 'lam1' in df.columns and 'iteration' in df.columns:
        fig, ax = plt.subplots(figsize=(12, 6))
        df_plot = df.sort_values('iteration')
        ax.scatter(df_plot.index, df_plot['lam1'], alpha=0.3, s=10)
        ax.set_xlabel('Evaluation Index')
        ax.set_ylabel('λ₁')
        ax.set_title('λ₁ Coverage Over Campaign')
        fig.tight_layout()
        fig.savefig(output_dir / "lambda1_coverage.png", dpi=150)
        plt.close(fig)
        print(f"  ✓ Saved: {output_dir / 'lambda1_coverage.png'}")
    
    # 2. Constraint satisfaction over iterations
    constraint_cols = ['positivity_ok', 'unitarity_ok', 'perturbativity_ok']
    if all(c in df.columns for c in constraint_cols) and 'iteration' in df.columns:
        fig, ax = plt.subplots(figsize=(12, 6))
        
        iter_order = sorted(df['iteration'].unique())
        iter_stats = []
        for it in iter_order:
            iter_df = df[df['iteration'] == it]
            iter_stats.append({
                'iteration': it,
                'positivity': iter_df['positivity_ok'].mean() * 100,
                'unitarity': iter_df['unitarity_ok'].mean() * 100,
                'perturbativity': iter_df['perturbativity_ok'].mean() * 100,
                'all': iter_df[constraint_cols].all(axis=1).mean() * 100
            })
        
        stats_df = pd.DataFrame(iter_stats)
        x = range(len(stats_df))
        
        ax.plot(x, stats_df['positivity'], label='Positivity', marker='o')
        ax.plot(x, stats_df['unitarity'], label='Unitarity', marker='s')
        ax.plot(x, stats_df['perturbativity'], label='Perturbativity', marker='^')
        ax.plot(x, stats_df['all'], label='All Constraints', marker='D', linewidth=2)
        
        ax.set_xlabel('Iteration')
        ax.set_ylabel('Pass Rate (%)')
        ax.set_title('Constraint Satisfaction Evolution')
        ax.legend()
        ax.set_xticks(x[::2])
        ax.set_xticklabels([stats_df.iloc[i]['iteration'] for i in x[::2]], rotation=45)
        fig.tight_layout()
        fig.savefig(output_dir / "constraint_evolution.png", dpi=150)
        plt.close(fig)
        print(f"  ✓ Saved: {output_dir / 'constraint_evolution.png'}")
    
    # 3. Parameter space coverage (lambda1 vs m_phi)
    if 'lam1' in df.columns and 'm_phi' in df.columns:
        df['all_constraints_ok'] = df[constraint_cols].all(axis=1) if all(c in df.columns for c in constraint_cols) else True
        
        fig, ax = plt.subplots(figsize=(10, 8))
        viable = df[df['all_constraints_ok'] == True]
        invalid = df[df['all_constraints_ok'] == False]
        
        ax.scatter(invalid['lam1'], invalid['m_phi'], alpha=0.3, s=10, c='red', label='Invalid')
        ax.scatter(viable['lam1'], viable['m_phi'], alpha=0.5, s=15, c='green', label='Viable')
        
        ax.set_xlabel('λ₁')
        ax.set_ylabel('m_φ (GeV)')
        ax.set_title('Parameter Space Coverage (λ₁ vs m_φ)')
        ax.legend()
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(output_dir / "parameter_space.png", dpi=150)
        plt.close(fig)
        print(f"  ✓ Saved: {output_dir / 'parameter_space.png'}")
    
    # 4. Lambda1 histogram
    if 'lam1' in df.columns:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(df['lam1'], bins=50, alpha=0.7, edgecolor='black')
        ax.set_xlabel('λ₁')
        ax.set_ylabel('Frequency')
        ax.set_title('λ₁ Distribution')
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(output_dir / "lambda1_histogram.png", dpi=150)
        plt.close(fig)
        print(f"  ✓ Saved: {output_dir / 'lambda1_histogram.png'}")
    
    print(f"\n📊 All plots saved to: {output_dir}/")


def main():
    parser = argparse.ArgumentParser(
        description="Aggregate and visualize autoresearch campaign results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic aggregation and summary
  python aggregate_results.py runs/my-campaign
  
  # Generate plots
  python aggregate_results.py runs/my-campaign --plots
  
  # Save aggregated CSV
  python aggregate_results.py runs/my-campaign --output all_results.csv
  
  # Full analysis
  python aggregate_results.py runs/my-campaign --plots --output all_results.csv
        """
    )
    
    parser.add_argument(
        'campaign_dir',
        type=Path,
        help='Path to campaign output directory (contains checkpoints/)'
    )
    parser.add_argument(
        '--output', '-o',
        type=Path,
        help='Save aggregated results to CSV file'
    )
    parser.add_argument(
        '--plots', '-p',
        action='store_true',
        help='Generate visualization plots'
    )
    
    args = parser.parse_args()
    
    if not args.campaign_dir.exists():
        print(f"ERROR: Campaign directory not found: {args.campaign_dir}", file=sys.stderr)
        sys.exit(1)
    
    print(f"🔍 Scanning campaign: {args.campaign_dir}")
    
    # Find all results files
    results_files = find_results_files(args.campaign_dir)
    print(f"   Found {len(results_files)} results.csv files")
    
    if not results_files:
        print("ERROR: No results.csv files found in campaign checkpoints", file=sys.stderr)
        sys.exit(1)
    
    # Aggregate
    print(f"📦 Aggregating results...")
    df = aggregate_results(results_files)
    
    if df.empty:
        print("ERROR: Failed to load any results", file=sys.stderr)
        sys.exit(1)
    
    print(f"   Loaded {len(df)} evaluations")
    
    # Generate summary
    summary = summarize_results(df)
    print_summary(summary)
    
    # Save aggregated CSV
    if args.output:
        df.to_csv(args.output, index=False)
        print(f"💾 Saved aggregated results to: {args.output}")
    
    # Generate plots
    if args.plots:
        print(f"📊 Generating plots...")
        plot_results(df, args.campaign_dir)
    
    print("✅ Analysis complete!")


if __name__ == '__main__':
    main()
