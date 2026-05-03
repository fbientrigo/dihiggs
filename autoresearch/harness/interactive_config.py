#!/usr/bin/env python3
"""
Interactive CLI wizard for generating autoresearch campaign configurations.

Usage:
    python interactive_config.py --wizard
    python interactive_config.py --template lambda1
    python interactive_config.py --edit existing_config.json
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, Optional


# Template configurations
TEMPLATES = {
    "lambda1": {
        "name": "Lambda1 Exploration",
        "description": "Explore lambda1 parameter space with fixed m_phi, tan_beta",
        "config": {
            "supervisor": {
                "max_rounds": 0,
                "enable_convergence": True,
                "checkpoint_interval": 1
            },
            "runtime": {
                "timeout_sec": 300,
                "threads": 4
            },
            "limits": {
                "max_new_run_dirs_per_round": 100,
                "max_new_run_dirs_per_arm_call": 100
            },
            "convergence": {
                "plateau_threshold": 0.01,
                "plateau_window": 3,
                "max_empty_rounds": 2
            },
            "autoscaling": {
                "enabled": True,
                "scale_threads": True,
                "scale_timeout": True,
                "min_threads": 2,
                "max_threads": 8,
                "min_timeout_sec": 60,
                "max_timeout_sec": 600
            },
            "bandit_config": {
                "ucb_c": 1.414
            },
            "orchestrator": {
                "arms": [
                    {
                        "name": "adaptive-smoke",
                        "command": [
                            "{python}",
                            "{repo_root}/dihiggs/app/adaptive_explorer_lam1.py",
                            "--m-phi", "130",
                            "--m-A", "300",
                            "--m-Hp", "300",
                            "--tan-beta", "1000",
                            "--lambda6", "0.001",
                            "--lambda7", "0",
                            "--sin-ba", "1.0",
                            "--lam1-min", "0",
                            "--lam1-max", "12",
                            "--lam1-bins", "2000",
                            "--hb-dataset", "{hb_dataset}",
                            "--hs-dataset", "{hs_dataset}",
                            "--n-proposals", "100",
                            "--output-dir", "{checkpoint_root}/{arm_name}/iter_{iter:04d}"
                        ],
                        "env": {}
                    }
                ]
            }
        }
    },
    "mphi_scan": {
        "name": "m_phi Grid Scan",
        "description": "Scan m_phi parameter space with fixed lambda parameters",
        "config": {
            "supervisor": {
                "max_rounds": 0,
                "enable_convergence": True,
                "checkpoint_interval": 1
            },
            "runtime": {
                "timeout_sec": 300,
                "threads": 4
            },
            "limits": {
                "max_new_run_dirs_per_round": 50,
                "max_new_run_dirs_per_arm_call": 50
            },
            "convergence": {
                "plateau_threshold": 0.01,
                "plateau_window": 3,
                "max_empty_rounds": 2
            },
            "autoscaling": {
                "enabled": True,
                "scale_threads": True,
                "scale_timeout": True,
                "min_threads": 2,
                "max_threads": 8,
                "min_timeout_sec": 60,
                "max_timeout_sec": 600
            },
            "bandit_config": {
                "ucb_c": 1.414
            },
            "orchestrator": {
                "arms": [
                    {
                        "name": "mphi-scan",
                        "command": [
                            "{python}",
                            "{repo_root}/dihiggs/app/adaptive_explorer.py",
                            "--m-phi-min", "130",
                            "--m-phi-max", "300",
                            "--m-A", "300",
                            "--m-Hp", "300",
                            "--tan-beta", "1000",
                            "--lambda6", "0.001",
                            "--lambda7", "0",
                            "--sin-ba", "1.0",
                            "--lam1", "6.0",
                            "--hb-dataset", "{hb_dataset}",
                            "--hs-dataset", "{hs_dataset}",
                            "--n-proposals", "50",
                            "--output-dir", "{checkpoint_root}/{arm_name}/iter_{iter:04d}"
                        ],
                        "env": {}
                    }
                ]
            }
        }
    }
}


def prompt_value(prompt: str, default: Any = None, value_type: type = str) -> Any:
    """Prompt user for a value with optional default."""
    if default is not None:
        prompt_str = f"{prompt} [{default}]: "
    else:
        prompt_str = f"{prompt}: "
    
    while True:
        response = input(prompt_str).strip()
        
        if not response and default is not None:
            return default
        
        if not response:
            print("  ⚠️  Value required. Please try again.")
            continue
        
        try:
            if value_type == bool:
                return response.lower() in ('y', 'yes', 'true', '1')
            elif value_type == int:
                return int(response)
            elif value_type == float:
                return float(response)
            else:
                return response
        except ValueError:
            print(f"  ⚠️  Invalid {value_type.__name__}. Please try again.")


def wizard_mode() -> Dict:
    """Interactive wizard to build a campaign configuration."""
    print("\n" + "="*60)
    print("AUTORESEARCH CAMPAIGN CONFIG WIZARD")
    print("="*60 + "\n")
    
    # Choose template
    print("Available templates:")
    for i, (key, tmpl) in enumerate(TEMPLATES.items(), 1):
        print(f"  {i}. {tmpl['name']}")
        print(f"     {tmpl['description']}")
    
    print(f"  {len(TEMPLATES) + 1}. Custom (start from scratch)")
    
    tmpl_choice = prompt_value(f"\nSelect template (1-{len(TEMPLATES) + 1})", default=1, value_type=int)
    
    if 1 <= tmpl_choice <= len(TEMPLATES):
        template_key = list(TEMPLATES.keys())[tmpl_choice - 1]
        config = TEMPLATES[template_key]["config"].copy()
        print(f"\n✓ Using template: {TEMPLATES[template_key]['name']}\n")
    else:
        config = {
            "supervisor": {},
            "runtime": {},
            "limits": {},
            "convergence": {},
            "autoscaling": {},
            "bandit_config": {},
            "orchestrator": {"arms": []}
        }
        print("\n✓ Starting custom configuration\n")
    
    # Campaign basics
    print("=== Campaign Settings ===\n")
    
    max_rounds = prompt_value("Max rounds (0 = run until convergence)", 
                              default=config.get("supervisor", {}).get("max_rounds", 0), 
                              value_type=int)
    config["supervisor"]["max_rounds"] = max_rounds
    
    enable_conv = prompt_value("Enable convergence detection? (y/n)", 
                                default="y" if config.get("supervisor", {}).get("enable_convergence", True) else "n",
                                value_type=bool)
    config["supervisor"]["enable_convergence"] = enable_conv
    
    # Time budget options
    use_time_budget = prompt_value("Use time budget (instead of/in addition to max_rounds)? (y/n)", default="n", value_type=bool)
    
    if use_time_budget:
        max_duration_hours = prompt_value("Maximum campaign duration (hours)", default=8.0, value_type=float)
        config["supervisor"]["max_duration_hours"] = max_duration_hours
        print(f"  ✓ Campaign will stop after {max_duration_hours} hours")
    else:
        # Remove time budget if present
        config["supervisor"].pop("max_duration_hours", None)
        config["supervisor"].pop("stop_at_timestamp", None)
    
    # Runtime settings
    print("\n=== Runtime Settings ===\n")
    
    threads = prompt_value("Number of parallel threads", 
                           default=config.get("runtime", {}).get("threads", 4), 
                           value_type=int)
    config["runtime"]["threads"] = threads
    
    timeout = prompt_value("Timeout per subprocess (seconds)", 
                           default=config.get("runtime", {}).get("timeout_sec", 300), 
                           value_type=int)
    config["runtime"]["timeout_sec"] = timeout
    
    # Limits
    print("\n=== Resource Limits ===\n")
    
    max_per_round = prompt_value("Max evaluations per round", 
                                  default=config.get("limits", {}).get("max_new_run_dirs_per_round", 100), 
                                  value_type=int)
    config["limits"]["max_new_run_dirs_per_round"] = max_per_round
    config["limits"]["max_new_run_dirs_per_arm_call"] = max_per_round
    
    # Physics parameters (if using template with arms)
    if config.get("orchestrator", {}).get("arms"):
        print("\n=== Physics Parameters ===\n")
        print("(Modifying first arm configuration)\n")
        
        arm = config["orchestrator"]["arms"][0]
        cmd = arm["command"]
        
        # Helper to update command arg
        def update_arg(flag: str, prompt_text: str, default_val: Any, val_type: type = float):
            try:
                idx = cmd.index(flag)
                current = cmd[idx + 1]
                new_val = prompt_value(prompt_text, default=current, value_type=val_type)
                cmd[idx + 1] = str(new_val)
            except (ValueError, IndexError):
                pass
        
        # Lambda1 range (if present)
        if "--lam1-min" in cmd:
            update_arg("--lam1-min", "Lambda1 minimum", 0, float)
            update_arg("--lam1-max", "Lambda1 maximum", 12, float)
            update_arg("--lam1-bins", "Lambda1 bins (resolution)", 2000, int)
        
        # m_phi (fixed or range)
        if "--m-phi" in cmd:
            update_arg("--m-phi", "m_phi (GeV)", 130, float)
        elif "--m-phi-min" in cmd:
            update_arg("--m-phi-min", "m_phi minimum (GeV)", 130, float)
            update_arg("--m-phi-max", "m_phi maximum (GeV)", 300, float)
        
        # Other parameters
        update_arg("--tan-beta", "tan(beta)", 1000, float)
        update_arg("--lambda6", "Lambda6", 0.001, float)
        update_arg("--lambda7", "Lambda7", 0, float)
        
        # Proposals per round
        if "--n-proposals" in cmd:
            update_arg("--n-proposals", "Proposals per round", 100, int)
    
    # Output
    print("\n=== Output Settings ===\n")
    
    output_name = prompt_value("Campaign output directory name", default="my-campaign")
    
    return config, output_name


def edit_mode(config_path: Path) -> Dict:
    """Load and interactively edit an existing configuration."""
    print(f"\n📝 Loading: {config_path}\n")
    
    with open(config_path) as f:
        config = json.load(f)
    
    print("Current configuration loaded. You can now modify key parameters.\n")
    print("=== Quick Edit Mode ===\n")
    
    # Show current values and allow edits
    max_rounds = config.get("supervisor", {}).get("max_rounds", 0)
    new_max_rounds = prompt_value(f"Max rounds", default=max_rounds, value_type=int)
    config.setdefault("supervisor", {})["max_rounds"] = new_max_rounds
    
    threads = config.get("runtime", {}).get("threads", 4)
    new_threads = prompt_value(f"Threads", default=threads, value_type=int)
    config.setdefault("runtime", {})["threads"] = new_threads
    
    max_per_round = config.get("limits", {}).get("max_new_run_dirs_per_round", 100)
    new_max = prompt_value(f"Max evaluations per round", default=max_per_round, value_type=int)
    config.setdefault("limits", {})["max_new_run_dirs_per_round"] = new_max
    config.setdefault("limits", {})["max_new_run_dirs_per_arm_call"] = new_max
    
    # Edit first arm if present
    if config.get("orchestrator", {}).get("arms"):
        print("\n=== Arm Parameters ===\n")
        arm = config["orchestrator"]["arms"][0]
        cmd = arm["command"]
        
        def update_arg(flag: str, prompt_text: str, val_type: type = str):
            try:
                idx = cmd.index(flag)
                current = cmd[idx + 1]
                new_val = prompt_value(f"{prompt_text}", default=current, value_type=val_type)
                cmd[idx + 1] = str(new_val)
            except (ValueError, IndexError):
                pass
        
        if "--n-proposals" in cmd:
            update_arg("--n-proposals", "Proposals per round", int)
        
        if "--lam1-min" in cmd:
            update_arg("--lam1-min", "Lambda1 min", float)
            update_arg("--lam1-max", "Lambda1 max", float)
    
    return config, config_path.stem


def save_config(config: Dict, output_path: Path):
    """Save configuration to JSON file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w') as f:
        json.dump(config, f, indent=2)
    
    print(f"\n✅ Configuration saved to: {output_path}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Interactive campaign configuration builder",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Interactive wizard (recommended for new campaigns)
  python interactive_config.py --wizard
  
  # Start from a template
  python interactive_config.py --template lambda1
  
  # Edit existing config
  python interactive_config.py --edit configs/my_campaign.json
  
  # List available templates
  python interactive_config.py --list-templates
        """
    )
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--wizard', '-w', action='store_true',
                       help='Launch interactive wizard')
    group.add_argument('--template', '-t', choices=list(TEMPLATES.keys()),
                       help='Use a template as starting point')
    group.add_argument('--edit', '-e', type=Path,
                       help='Edit existing configuration file')
    group.add_argument('--list-templates', '-l', action='store_true',
                       help='List available templates')
    
    parser.add_argument('--output', '-o', type=Path,
                        help='Output file path (default: auto-generated in configs/)')
    
    args = parser.parse_args()
    
    # List templates
    if args.list_templates:
        print("\nAvailable templates:\n")
        for key, tmpl in TEMPLATES.items():
            print(f"  {key}")
            print(f"    Name: {tmpl['name']}")
            print(f"    Description: {tmpl['description']}")
            print()
        return
    
    # Wizard mode
    if args.wizard:
        config, output_name = wizard_mode()
    
    # Template mode
    elif args.template:
        template = TEMPLATES[args.template]
        print(f"\n✓ Using template: {template['name']}")
        print(f"  {template['description']}\n")
        config = template["config"]
        output_name = args.template
    
    # Edit mode
    elif args.edit:
        if not args.edit.exists():
            print(f"ERROR: File not found: {args.edit}", file=sys.stderr)
            sys.exit(1)
        config, output_name = edit_mode(args.edit)
    
    # Determine output path
    if args.output:
        output_path = args.output
    else:
        configs_dir = Path(__file__).parent.parent / "configs"
        output_path = configs_dir / f"{output_name}.json"
    
    # Save
    save_config(config, output_path)
    
    # Show next steps
    print("Next steps:")
    print(f"  1. Review: cat {output_path}")
    print(f"  2. Run campaign:")
    print(f"     python autoresearch/run_supervisor.py {output_path}")
    print(f"  3. Monitor:")
    print(f"     python autoresearch/harness/live_monitor.py --watch runs/{output_name}")
    print()


if __name__ == '__main__':
    main()
