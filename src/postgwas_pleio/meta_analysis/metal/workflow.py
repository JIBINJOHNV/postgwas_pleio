#from postgwas_pleio.meta_analysis.metal.metal_direct import metal_direct_runner
from postgwas_pleio.meta_analysis.metal.metal_pipeline import metal_pipeline_runner
from postgwas_pleio.meta_analysis.metal.metal_direct import metal_direct_runner
import logging

logger = logging.getLogger(__name__)


def metal_workflow(args, remainder=None):
    """
    remainder will be None for pipeline mode, 
    but will contain a list of strings for direct mode.
    """
    if args.mode == "pipeline":
        logger.info("Executing Metal Pipeline Runner...")
        print("Running Pipeline...")
        metal_pipeline_runner(args)

    elif args.mode == "direct":
        logger.info("Executing Metal Direct Passthrough...")
        # Pass the remainder list to the runner
        print("Running Direct...")
        metal_direct_runner(args)