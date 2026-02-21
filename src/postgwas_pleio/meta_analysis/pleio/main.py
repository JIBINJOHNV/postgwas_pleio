import sys
import argparse
from postgwas_pleio.clis.metanalysis_cli import parse_pleio_pipeline_args, parse_pleio_direct_args
from postgwas_pleio.meta_analysis.pleio.workflow import pleio_workflow
import logging

logger = logging.getLogger(__name__)

def main():
    parser = argparse.ArgumentParser(
        prog="postgwas-pleio",
        description="PLEIO: Pleiotropy-Informed LD Score Regression",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    subparsers = parser.add_subparsers(dest="mode", help="Execution modes", required=True)

    parse_pleio_pipeline_args(subparsers)
    parse_pleio_direct_args(subparsers)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    # ⭐ PASS RAW TOKENS
    remainder = sys.argv[2:]

    pleio_workflow(args, remainder)

if __name__ == "__main__":
    main()