import subprocess
import logging
import sys  # Added for flushing
from postgwas_pleio.formatter.asset_formatter import asset_formatter
from pathlib import Path
import os


logger = logging.getLogger(__name__)


def run_cmd_with_log(cmd, log_file, step_name):
    """
    Run subprocess command with tee-style logging (console + file).
    """
    logger.info(f"Running {step_name}")
    print(f"[LOG] Writing {step_name} output → {log_file}", flush=True)

    with open(log_file, "w") as log:
        log.write(f"{step_name} COMMAND:\n")
        log.write(" ".join(map(str, cmd)) + "\n\n")

        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1
        )

        if process.stdout is not None:
            for line in process.stdout:
                print(line, end="")     # console
                log.write(line)         # file

        process.wait()

        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode, cmd)


def asset_pipeline_runner(args):
    """
    Executes ASSET and renames the generic trait_N outputs to original sample IDs.
    """
    out_path = Path(args.out)
    out_path.mkdir(parents=True, exist_ok=True)

    # -------------------------------------------------
    # STEP 1: Formatter
    # -------------------------------------------------
    print(f"\n[STEP 1] Formatting VCFs into ASSET format...", flush=True)

    formatter_results = asset_formatter(
        inputfile=args.inputfile,
        output_folder=str(args.out),
        run_name=str(args.run_name)
    )

    log_file = formatter_results["log_file"]

    # -------------------------------------------------
    # STEP 2: ASSET
    # -------------------------------------------------
    from importlib.resources import files
    asset_cli = str(files("postgwas_pleio.meta_analysis.asset.R") / "asset_cli.R")
    
    asset_log = out_path / f"{args.run_name}_asset.log"

    cmd_asset = [
        "micromamba", "run", "-n", "postgwas",
        "Rscript", str(asset_cli),
        "--input_manifest", formatter_results["fastasset_sample_manifest"],
        "--fastasset_input", formatter_results["fastasset_input"],
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

    print(f"[STEP 2] Starting ASSET...", flush=True)

    run_cmd_with_log(cmd_asset, asset_log, "ASSET")

    # -------------------------------------------------
    # STEP 4: Return metadata
    # -------------------------------------------------
    fastasset_1sided = Path(out_path, "5_fastasset_output", f"{args.run_name}_fastasset_1sided.tsv.gz")
    fastasset_2sided = Path(out_path, "5_fastasset_output", f"{args.run_name}_fastasset_2sided.tsv.gz")

    asset_1sided = Path(out_path, "6_htraits_output", f"{args.run_name}_htraits_1sided.tsv.gz")
    asset_2sided = Path(out_path, "6_htraits_output", f"{args.run_name}_htraits_2sided.tsv.gz")
    asset_meta = Path(out_path, "6_htraits_output", f"{args.run_name}_htraits_meta.tsv.gz")

    return {
        "sample_ids": str(args.run_name),
        "fastasset_1sided": str(fastasset_1sided),  # removed trailing space bug
        "fastasset_2sided": str(fastasset_2sided),
        "asset_1sided": str(asset_1sided),
        "asset_2sided": str(asset_2sided),
        "asset_meta": str(asset_meta),
        "log_file": log_file,
    }