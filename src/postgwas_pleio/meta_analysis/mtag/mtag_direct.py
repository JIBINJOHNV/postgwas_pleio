import subprocess
import logging
import os
import shutil
from pathlib import Path
from typing import List, Optional

logger = logging.getLogger(__name__)

def validate_paths(paths: List[Optional[str]], description: str):
    """Checks if files exist before starting MTAG."""
    for p in paths:
        if p:
            # Handle comma-separated strings (common in MTAG)
            individual_paths = p.split(",") if "," in p else [p]
            for path_str in individual_paths:
                path_obj = Path(path_str.strip())
                if not path_obj.exists():
                    logger.error(f"Validation Failed: {description} not found at {path_obj}")
                    raise FileNotFoundError(f"Missing {description}: {path_obj}")
    logger.info(f"Validation Passed: {description}")


def mtag_direct_runner(args):
    """
    Takes the parsed 'args' Namespace and rebuilds the command 
    to call the external mtag.py script.
    """
    cmd = ["micromamba", "run", "-n", "mtag", "python", "/opt/mtag/mtag.py"]
    
    # Convert the 'args' object back into a list of strings for subprocess
    # We skip 'mode' because mtag.py doesn't know what that is
    for key, value in vars(args).items():
        if key == "mode" or value is None or value is False:
            continue
        
        # Format the flag (e.g., 'sumstats' becomes '--sumstats')
        flag = f"--{key.replace('_', '-')}"
        
        if isinstance(value, bool): # For action="store_true"
            cmd.append(flag)
        else:
            cmd.append(flag)
            cmd.append(str(value))

    subprocess.run(cmd, check=True)