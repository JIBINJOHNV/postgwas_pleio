from postgwas_pleio.meta_analysis.pleio.pleio_direct import pleio_direct_runner
from postgwas_pleio.meta_analysis.pleio.pleio_pipeline import pleio_pipeline_runner
import logging

logger = logging.getLogger(__name__)


def pleio_workflow(args, remainder=None):
    """
    remainder will be None for pipeline mode, 
    but will contain a list of strings for direct mode.
    """
    if args.mode == "pipeline":
        logger.info("Executing Pleio Pipeline Runner...")
        print("Running Pipeline...")
        pleio_pipeline_runner(args)

    elif args.mode == "direct":
        logger.info("Executing Pleio Direct Passthrough...")
        # Pass the remainder list to the runner
        print("Running Direct...")
        pleio_direct_runner(args)