from postgwas_pleio.meta_analysis.asset.asset_direct import asset_direct_runner
from postgwas_pleio.meta_analysis.asset.asset_pipeline import asset_pipeline_runner
import logging

logger = logging.getLogger(__name__)


def asset_workflow(args, remainder=None):
    """
    remainder will be None for pipeline mode, 
    but will contain a list of strings for direct mode.
    """
    if args.mode == "pipeline":
        logger.info("Executing Asset Pipeline Runner...")
        print("Running Pipeline...")
        asset_pipeline_runner(args)

    elif args.mode == "direct":
        logger.info("Executing Asset Direct Passthrough...")
        # Pass the remainder list to the runner
        print("Running Direct...")
        asset_direct_runner(args)