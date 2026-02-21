import subprocess
import logging
import os
import shutil
from pathlib import Path
from typing import List, Optional

logger = logging.getLogger(__name__)


def run_cmd_with_log(cmd, log_file, step_name):
    """
    Run subprocess command with tee-style logging (console + file).
    """
    logger.info(f"Running {step_name}")
    print(f"[LOG] Writing {step_name} output → {log_file}", flush=True)

    with open(log_file, "w") as log:
        log.write(f"{step_name} COMMAND:\n")
        log.write(" ".join(map(str, cmd)) + "\n\n")   # safe join

        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1
        )

        if process.stdout is not None:                 # safety fix
            for line in process.stdout:
                print(line, end="")
                log.write(line)

        process.wait()

        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode, cmd)


def validate_paths(paths: List[Optional[str]], description: str):
    """Checks if files exist before starting MTAG."""
    for p in paths:
        if p:
            individual_paths = p.split(",") if "," in p else [p]
            for path_str in individual_paths:
                path_obj = Path(path_str.strip())
                if not path_obj.exists():
                    logger.error(f"Validation Failed: {description} not found at {path_obj}")
                    raise FileNotFoundError(f"Missing {description}: {path_obj}")
    logger.info(f"Validation Passed: {description}")


def asset_direct_runner(args):
    """
    Takes parsed args Namespace and calls ASSET CLI.
    """
    out_path = Path(args.output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    # -------------------------------------------------
    # STEP 1: ASSET
    # -------------------------------------------------

    from importlib.resources import files
    asset_cli = str(files("postgwas_pleio.meta_analysis.asset.R") / "asset_cli.R")
    
    asset_log = out_path / f"{args.run_name}_asset.log"

    cmd_asset = [
        "micromamba", "run", "-n", "postgwas",
        "Rscript", str(asset_cli),   # fix
        "--input_manifest", args.input_manifest,
        "--fastasset_input", args.fastasset_input,
        "--output_dir", str(out_path),
        "--hm3", str(args.hm3),
        "--ld_ref", str(args.ld_ref_panel),
        "--info_filter", str(args.info_filter),
        "--maf_filter", str(args.maf_filter),
        "--chunk_size", str(args.chunk_size),
        "--scr_pthr", str(args.scr_pthr),
        "--meth_pval", str(args.meth_pval),
        "--ncores", str(args.ncores),
        "--run_name", str(args.run_name),
    ]

    print(f"[STEP 1] Starting ASSET...", flush=True)

    run_cmd_with_log(cmd_asset, asset_log, "ASSET")

    # -------------------------------------------------
    # STEP 2: Return metadata
    # -------------------------------------------------
    fastasset_1sided = Path(out_path, "5_fastasset_output", f"{args.run_name}_fastasset_1sided.tsv.gz")
    fastasset_2sided = Path(out_path, "5_fastasset_output", f"{args.run_name}_fastasset_2sided.tsv.gz")

    asset_1sided = Path(out_path, "6_htraits_output", f"{args.run_name}_htraits_1sided.tsv.gz")
    asset_2sided = Path(out_path, "6_htraits_output", f"{args.run_name}_htraits_2sided.tsv.gz")
    asset_meta = Path(out_path, "6_htraits_output", f"{args.run_name}_htraits_meta.tsv.gz")

    return {
        "sample_ids": str(args.run_name),
        "fastasset_1sided": str(fastasset_1sided),   # fixed
        "fastasset_2sided": str(fastasset_2sided),
        "asset_1sided": str(asset_1sided),
        "asset_2sided": str(asset_2sided),
        "asset_meta": str(asset_meta),
        "log_file": str(asset_log),                  # fixed
    }