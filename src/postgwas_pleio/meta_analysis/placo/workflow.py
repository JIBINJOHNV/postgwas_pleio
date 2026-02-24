from postgwas_pleio.meta_analysis.placo.placo_direct import placo_direct_runner
from postgwas_pleio.meta_analysis.placo.placo_pipeline import placo_pipeline_runner
import logging

logger = logging.getLogger(__name__)


def placo_workflow(args, remainder=None):
    """
    remainder will be None for pipeline mode, 
    but will contain a list of strings for direct mode.
    """
    if args.mode == "pipeline":
        logger.info("Executing Placo Pipeline Runner...")
        print("Running Pipeline...")
        placo_pipeline_runner(args)

    elif args.mode == "direct":
        logger.info("Executing Placo Direct Passthrough...")
        # Pass the remainder list to the runner
        print("Running Direct...")
        placo_direct_runner(args)