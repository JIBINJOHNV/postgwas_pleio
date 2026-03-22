import sys
import argparse
from postgwas_pleio.clis.metanalysis_cli import parse_metal_pipeline_args
from postgwas_pleio.meta_analysis.metal.workflow import metal_workflow



def main():
    parser = argparse.ArgumentParser(
        prog="postgwas-pleio meta-analysis metal",
        description="METAL: Meta-Analysis of GWAS (Metal et al.)",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # -------------------------------------------------
    # Subparsers
    # -------------------------------------------------
    subparsers = parser.add_subparsers(dest="mode", help="Execution modes")

    parse_metal_pipeline_args(subparsers)
    # parse_metal_direct_args(subparsers)

    # -------------------------------------------------
    # Show help if no arguments
    # -------------------------------------------------
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    # -------------------------------------------------
    # Parse arguments
    # -------------------------------------------------
    args = parser.parse_args()

    # -------------------------------------------------
    # Validation
    # -------------------------------------------------
    if hasattr(args, "scheme") and args.scheme.upper() == "STDERR":
        if args.sample_size_approach.lower() != "totalnef":
            parser.error("--scheme STDERR requires --sample_size_approach totalnef")

    # -------------------------------------------------
    # Route execution
    # -------------------------------------------------
    if args.mode == "pipeline":
        metal_workflow(args)

    elif args.mode == "direct":
        metal_workflow(args)

    else:
        parser.print_help()


if __name__ == "__main__":
    main()