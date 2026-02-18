#!/usr/bin/env python3
"""
postgwas-pleio — Unified CLI for Pleiotropy and Multi-trait analysis.
"""

import sys
import logging
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich.text import Text

# Initialize Rich Console
console = Console()

# 1. Import Category Dispatchers 
# Using a try-except here prevents the whole tool from failing if one module is broken
from postgwas_pleio.meta_analysis.meta_analysis_cli import main as meta_analysis_main

# -------------------------------------------------------------------
# CATEGORY REGISTRY
# -------------------------------------------------------------------
# Format: "command": (function_to_call, "description")
CATEGORIES = {
    "meta-analysis": (meta_analysis_main, "Multi-trait analysis (MTAG, ASSET, FastASSET, Pleio)"),
    "colocalisation": (None, "Genetic colocalisation analysis (e.g., coloc, eCAVIAR)"),
    "twas": (None, "Transcriptome-wide association studies"),
    "conditional-analysis": (None, "Conditional & joint association analysis (GCTA-COJO)"),
    "mr": (None, "Mendelian Randomisation analysis"),
}

# Sort categories alphabetically for a consistent, professional help menu
CATEGORIES = dict(sorted(CATEGORIES.items(), key=lambda x: x[0].lower()))

def print_global_help() -> None:
    """Render the high-level Rich help menu."""
    
    header = Text("PostGWAS-Pleio: Unified Pleiotropy Toolkit", style="bold white on blue")
    console.print(Panel(header, expand=False))

    console.print("[bold yellow]Available Analysis Categories:[/bold yellow]")
    
    table = Table(box=None, padding=(0, 2))
    table.add_column("Category", style="cyan", no_wrap=True)
    table.add_column("Description", style="white")

    for name, (_, desc) in CATEGORIES.items():
        table.add_row(name, desc)

    console.print(table)
    console.print("\n[bold yellow]Usage:[/bold yellow]  postgwas-pleio <category> <tool> [options]")
    console.print("Example: [green]postgwas-pleio meta-analysis mtag --help[/green]\n")

def main():
    argv = sys.argv
    prog = "postgwas-pleio"

    # 1. Handle help or empty execution at the root level
    if len(argv) == 1 or argv[1] in ("-h", "--help"):
        print_global_help()
        sys.exit(0)

    category_cmd = argv[1]

    # 2. Validate Category
    if category_cmd not in CATEGORIES:
        console.print(f"[bold red]Error:[/] '{category_cmd}' is not a valid category.", style="red")
        print_global_help()
        sys.exit(1)

    # 3. Dispatch to Category CLI
    category_func, _ = CATEGORIES[category_cmd]

    if category_func is None:
        console.print(f"[bold yellow]Module '{category_cmd}' is currently under development.[/]")
        sys.exit(0)

    # 4. Reconstruct sys.argv for the child CLI
    # This ensures the child CLI's 'prog' name is 'postgwas-pleio meta-analysis'
    # argv[0] is usually the script path; we replace it with the composite command
    sys.argv = [f"{prog} {category_cmd}"] + argv[2:]

    try:
        # Transfer control to the category-specific main()
        return category_func()
    
    except KeyboardInterrupt:
        console.print("\n[yellow]Execution interrupted by user.[/]")
        sys.exit(130)
    
    except Exception as e:
        # This catch-all ensures that even if a sub-module crashes, 
        # the user gets a semi-helpful error message instead of a raw traceback.
        console.print(f"[bold red]Critical Error in {category_cmd}:[/] {e}")
        # Only show the full traceback if the user might be a developer (check for a flag)
        if "--verbose" in argv:
            raise e
        sys.exit(1)

if __name__ == "__main__":
    main()