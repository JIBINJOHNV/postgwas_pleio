#from postgwas_pleio.meta_analysis.genomicpca.genomicpca_direct import genomicpca_direct_runner
from postgwas_pleio.meta_analysis.genomicpca.genomicpca_pipeline import genomicpca_pipeline_runner
import logging

logger = logging.getLogger(__name__)


def genomicpca_workflow(args, remainder=None):
    """
    remainder will be None for pipeline mode, 
    but will contain a list of strings for direct mode.
    """
    if args.mode == "pipeline":
        logger.info("Executing Asset Pipeline Runner...")
        print("Running Pipeline...")
        genomicpca_pipeline_runner(args)
    
    elif args.mode == "direct":
        print("direct mode not yet implimented")
    #     logger.info("Executing Asset Direct Passthrough...")
    #     # Pass the remainder list to the runner
    #     print("Running Direct...")
    #     genomicpca_direct_runner(args)