import sys
import argparse

from postgwas_pleio.clis.metanalysis_cli import (
    parse_genomicpca_pipeline_args,
)

from postgwas_pleio.meta_analysis.genomicpca.workflow import genomicpca_workflow


def main():

    parser = argparse.ArgumentParser(
        prog="postgwas-pleio meta-analysis genomicpca",
        description="GenomicPCA: Multivariate GWAS Meta-Analysis using genetic correlation structure",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # -----------------------------------------------------
    # Subparsers (pipeline / direct)
    # -----------------------------------------------------
    subparsers = parser.add_subparsers(
        dest="mode",
        help="GenomicPCA execution modes"
    )

    # CLI definitions
    parse_genomicpca_pipeline_args(subparsers)

    # -----------------------------------------------------
    # Show help if no arguments
    # -----------------------------------------------------
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    # -----------------------------------------------------
    # Parse arguments
    # -----------------------------------------------------
    args = parser.parse_args()

    # -----------------------------------------------------
    # Route execution mode
    # -----------------------------------------------------
    if args.mode == "pipeline":
        print("Running GenomicPCA pipeline...")
        genomicpca_workflow(args)

    elif args.mode == "direct":
        print("Running GenomicPCA direct mode...")
        genomicpca_workflow(args)

    else:
        parser.print_help()


if __name__ == "__main__":
    main()