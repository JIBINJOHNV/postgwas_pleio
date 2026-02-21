import sys
from rich.console import Console
from rich.table import Table

# Import tool entrypoints
# Ensure these paths match your directory structure
from postgwas_pleio.meta_analysis.mtag.main import main as mtag_main
from postgwas_pleio.meta_analysis.pleio.main import main as pleio_main
from postgwas_pleio.meta_analysis.asset.main import main as asset_main
console = Console()

# Registry of tools available under the 'meta-analysis' category
TOOLS = {
    "mtag": (mtag_main, "Multi-trait Analysis of GWAS (Turley et al. 2018)"),
    "asset": (asset_main, "Association Analysis based on Subsets (ASSET)"),
    "pleio": (pleio_main, "Pleiotropy-informed meta-analysis"),
}

def main():
    argv = sys.argv
    
    # 1. Handle help or empty execution for this category
    if len(argv) == 1 or argv[1] in ("-h", "--help"):
        table = Table(title="Meta-Analysis Tools", title_style="bold cyan")
        table.add_column("Tool", style="green", no_wrap=True)
        table.add_column("Description", style="white")
        
        for name, (_, desc) in TOOLS.items():
            table.add_row(name, desc)
            
        console.print(table)
        console.print("\n[bold yellow]Usage:[/bold yellow] postgwas-pleio meta-analysis <tool> [options]")
        return

    tool_cmd = argv[1]

    # 2. Validate Tool Selection
    if tool_cmd not in TOOLS:
        console.print(f"[bold red]Error:[/bold red] '{tool_cmd}' is not a valid tool in meta-analysis.", style="red")
        sys.exit(1)

    # 3. Dispatch to the specific Tool's main function
    tool_func, _ = TOOLS[tool_cmd]

    if tool_func is None:
        console.print(f"[bold yellow]Note: Tool '{tool_cmd}' is currently under development.[/bold yellow]")
        sys.exit(0)

    # 4. Reconstruct sys.argv for the tool-specific CLI
    # This shifts the arguments so the tool thinks it's the root (e.g., 'mtag' becomes argv[0])
    sys.argv = [f"postgwas-pleio meta-analysis {tool_cmd}"] + argv[2:]

    try:
        return tool_func()
    except KeyboardInterrupt:
        console.print("\n[yellow]Execution interrupted by user.[/yellow]")
        sys.exit(130)
    except Exception as e:
        console.print(f"[bold red]Critical Error in {tool_cmd}:[/bold red] {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()