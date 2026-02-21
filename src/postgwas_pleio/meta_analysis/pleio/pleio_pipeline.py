import subprocess
import logging
import sys, os
from postgwas_pleio.formatter.pleio_formatter import pleio_formatter
import pandas as pd
from pathlib import Path

logger = logging.getLogger(__name__)


def run_cmd_with_log(cmd, log_file, step_name):
    """
    Run subprocess command with tee-style logging (console + file).
    """
    logger.info(f"Running {step_name}")
    print(f"[LOG] Writing {step_name} output → {log_file}", flush=True)

    with open(log_file, "w") as log:
        log.write(f"{step_name} COMMAND:\n")
        log.write(" ".join(cmd) + "\n\n")

        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1
        )

        for line in process.stdout:
            print(line, end="")     # console
            log.write(line)         # file

        process.wait()

        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode, cmd)


def pleio_pipeline_runner(args):
    """
    Executes Pleio and renames the generic trait_N outputs to original sample IDs.
    """

    input_df = pd.read_csv(args.inputfile, sep=None, engine='python')
    sample_ids = input_df['sample_id'].tolist()

    out_path = Path(args.out)
    out_path.mkdir(parents=True, exist_ok=True)

    # -------------------------------------------------
    # STEP 1: Formatter
    # -------------------------------------------------
    print(f"\n[STEP 1] Formatting VCFs into Pleio format...", flush=True)

    formatter_results = pleio_formatter(
        inputfile=args.inputfile,
        output_folder=str(args.out),
        run_name=str(args.run_name)
    )

    PLEIO_ROOT = "/opt/pleio"
    log_file = formatter_results["log_file"]

    # -------------------------------------------------
    # STEP 2: LDSC preprocess
    # -------------------------------------------------
    ldsc_out_dir = out_path / "pleio_ldsc"
    ldsc_log = out_path / f"{args.run_name}_ldsc_preprocess.log"

    cmd_pleio_ldsc = [
        "micromamba", "run", "-n", "mtag",
        "python", os.path.join(PLEIO_ROOT, "ldsc_preprocess.py"),
        "--input", formatter_results["pleio_sample_manifest"],
        "--ref-ld-chr", str(args.ld_ref_panel),
        "--w-ld-chr", str(args.ld_ref_panel),
        "--out", str(ldsc_out_dir),
    ]

    print(f"[STEP 2] Starting Pleio ldsc_preprocess...", flush=True)

    run_cmd_with_log(cmd_pleio_ldsc, ldsc_log, "LDSC PREPROCESS")

    print("\n[SUCCESS] Pleio ldsc_preprocess finished.", flush=True)

    # -------------------------------------------------
    # STEP 3: Pleio meta
    # -------------------------------------------------
    pleio_log = out_path / f"{args.run_name}_pleio_meta.log"
    cmd_pleio_meta = [
        "micromamba", "run", "-n", "pleio",
        "python", os.path.join(PLEIO_ROOT, "pleio.py"),
        "--metain", str(ldsc_out_dir / "metain.txt.gz"),
        "--sg", str(ldsc_out_dir / "sg.txt.gz"),
        "--ce", str(ldsc_out_dir / "ce.txt.gz"),
        "--nis", str(args.nis),
        "--parallel",
        "--create",
        "--blup",
    ]
    # add optional flag
    if args.flattening_p_values:
        cmd_pleio_meta.append("--flattening_p_values")
    # always last
    cmd_pleio_meta.extend([
        "--out", str(out_path / f"{args.run_name}_pleio_results")
    ])

    print(f"[STEP 3] Starting Pleio meta-analysis...", flush=True)
    run_cmd_with_log(cmd_pleio_meta, pleio_log, "PLEIO META")
    print("\n[SUCCESS] Pleio meta-analysis finished.", flush=True)

    # -------------------------------------------------
    # STEP 4: Return metadata
    # -------------------------------------------------
    return {
        "sample_ids": sample_ids,
        "master_pleio_file": formatter_results["pleio_file"],
        "pleio_results_dir": str(out_path),
        "log_file": log_file,
        "ldsc_log": str(ldsc_log),
        "pleio_log": str(pleio_log)
    }
