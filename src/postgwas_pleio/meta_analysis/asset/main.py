import sys
import argparse
from postgwas_pleio.clis.metanalysis_cli import parse_asset_pipeline_args, parse_asset_direct_args
from postgwas_pleio.meta_analysis.asset.workflow import asset_workflow


def main():
    parser = argparse.ArgumentParser(
        prog="postgwas-pleio meta-analysis asset",
        description="ASSET: Multi-Trait Analysis of GWAS (Turley et al. 2018)",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Define subparsers
    subparsers = parser.add_subparsers(dest="mode", help="MTAG execution modes")

    # Call your two definition functions
    parse_asset_pipeline_args(subparsers)
    parse_asset_direct_args(subparsers)

    # Handle the help logic
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    # Parse everything at once
    args = parser.parse_args()

    # Route based on the 'mode' (pipeline or direct)
    if args.mode == "pipeline":
        print("Running Pipeline...")
        asset_workflow(args)
        # call mtag_workflow(args)
    elif args.mode == "direct":
        print("Running Direct Mode...")
        asset_workflow(args)
        # call mtag_workflow(args)

if __name__ == "__main__":
    main()