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
    
    # Define subparsers
    subparsers = parser.add_subparsers(dest="mode", help="MTAG execution modes")

    # Call your two definition functions
    parse_metal_pipeline_args(subparsers)
    #parse_metal_direct_args(subparsers)

    # Handle the help logic
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    # Parse everything at once
    args = parser.parse_args()

    # Route based on the 'mode' (pipeline or direct)
    if args.mode == "pipeline":
        print("Running Pipeline...")
        metal_workflow(args)
        # call mtag_workflow(args)
    elif args.mode == "direct":
        print("Running Direct Mode...")
        metal_workflow(args)
        # call mtag_workflow(args)

if __name__ == "__main__":
    main()