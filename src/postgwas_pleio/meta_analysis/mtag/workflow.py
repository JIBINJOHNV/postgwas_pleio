import logging
from postgwas_pleio.meta_analysis.mtag.mtag_pipeline import mtag_pipeline_runner
from postgwas_pleio.meta_analysis.mtag.mtag_direct import mtag_direct_runner

logger = logging.getLogger(__name__)

def mtag_workflow(args, remainder=None):
    """
    remainder will be None for pipeline mode, 
    but will contain a list of strings for direct mode.
    """
    if args.mode == "pipeline":
        logger.info("Executing MTAG Pipeline Runner...")
        print("Running Pipeline...")
        mtag_pipeline_runner(args)

    elif args.mode == "direct":
        logger.info("Executing MTAG Direct Passthrough...")
        # Pass the remainder list to the runner
        print("Running Direct...")
        mtag_direct_runner(remainder)