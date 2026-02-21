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



def pleio_direct_runner(args):

    if args.isf is None and not args.create:
        logger.info("No ISF provided → enabling --create automatically")
        args.create = True

    cmd = [
        "micromamba", "run", "-n", "pleio",
        "python", "/opt/pleio/pleio.py"
    ]

    for key, value in vars(args).items():

        if key == "mode" or value is None or value is False:
            continue

        flag = f"--{key}"

        if isinstance(value, bool):
            cmd.append(flag)
        else:
            cmd.extend([flag, str(value)])

    subprocess.run(cmd, check=True)



