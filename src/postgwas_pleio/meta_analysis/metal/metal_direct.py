
import subprocess
import logging
from pathlib import Path

from postgwas_pleio.meta_analysis.metal.script_generator import generate_metal_script



def metal_direct_runner(args):
    """
    Validate manifest → merge VCFs → create METAL per-study inputs → generate METAL script → run METAL.

    Notes:
      - METAL produces ONE meta-analysis output (not one output per study).
      - We therefore write sample_order.tsv to preserve input order / provenance.
    """
    pass