import sys
import argparse
from postgwas_pleio.clis.metanalysis_cli import parse_mtag_pipeline_args, parse_mtag_direct_args
from postgwas_pleio.meta_analysis.mtag.workflow import mtag_workflow


def main():
    parser = argparse.ArgumentParser(
        prog="postgwas-pleio meta-analysis mtag",
        description="MTAG: Multi-Trait Analysis of GWAS (Turley et al. 2018)",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Define subparsers
    subparsers = parser.add_subparsers(dest="mode", help="MTAG execution modes")

    # Call your two definition functions
    parse_mtag_pipeline_args(subparsers)
    parse_mtag_direct_args(subparsers)

    # Handle the help logic
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    # Parse everything at once
    args = parser.parse_args()

    # Route based on the 'mode' (pipeline or direct)
    if args.mode == "pipeline":
        print("Running Pipeline...")
        mtag_workflow(args)
        # call mtag_workflow(args)
    elif args.mode == "direct":
        print("Running Direct Mode...")
        mtag_workflow(args)
        # call mtag_workflow(args)

if __name__ == "__main__":
    main()