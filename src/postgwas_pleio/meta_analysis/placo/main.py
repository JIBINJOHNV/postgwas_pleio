import sys
import argparse
from postgwas_pleio.clis.metanalysis_cli import parse_placo_pipeline_args, parse_placo_direct_args
from postgwas_pleio.meta_analysis.placo.workflow import placo_workflow


def main():
    parser = argparse.ArgumentParser(
        prog="postgwas-pleio meta-analysis placo",
        description="PLACO: Polygenic Linkage Disequilibrium Adjusted Covariance",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Define subparsers
    subparsers = parser.add_subparsers(dest="mode", help="Execution modes", required=True)

    # Call your two definition functions
    parse_placo_pipeline_args(subparsers)
    parse_placo_direct_args(subparsers)

    # Handle the help logic
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    # Parse everything at once
    args = parser.parse_args()

    # Route based on the 'mode' (pipeline or direct)
    if args.mode == "pipeline":
        print("Running Pipeline...")
        placo_workflow(args)
        # call mtag_workflow(args)
    elif args.mode == "direct":
        print("Running Direct Mode...")
        placo_workflow(args)
        # call mtag_workflow(args)

if __name__ == "__main__":
    main()