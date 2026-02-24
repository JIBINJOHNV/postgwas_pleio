import os
import subprocess
import datetime
import logging
from pathlib import Path
from importlib.resources import files

from postgwas_pleio.utils.input_manifest_val import validate_manifest
from postgwas_pleio.formatter.merge_vcf_to_tsv import run_merge_vcf_to_tsv_pipeline
from postgwas_pleio.formatter.master_tsv_to_tool_format import run_mastertsv_to_toolformat_pipeline

logger = logging.getLogger(__name__)


def run_cmd_with_log(cmd, log_file, step_name):
    """
    Run subprocess command with tee-style logging (console + file).
    Path-safe: converts cmd items and log_file to str.
    """
    cmd = [str(x) for x in cmd]
    log_file = str(log_file)
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
            bufsize=1,
        )
        if process.stdout is not None:
            for line in process.stdout:
                print(line, end="")
                log.write(line)
        process.wait()
        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode, cmd)


def placo_pipeline_runner(args):
    """
    Validate manifest → merge VCFs → create PLACO input → run PLACO (R)
    """

    # ===========================
    # STEP 1: Validate manifest
    # ===========================
    out_path = Path(args.out)
    out_path.mkdir(parents=True, exist_ok=True)

    print("\n[STEP 1] Validating manifest...", flush=True)
    val_input = validate_manifest(args.inputfile, str(out_path))

    if len(val_input["manifest_df"]) > 2:
        raise ValueError("More than 2 samples found in the manifest.")
    
    # ===========================
    # STEP 2: Merge VCF → master TSV
    # ===========================
    print("\n[STEP 2] Merging VCFs → master TSV...", flush=True)

    preprocess_dir = out_path / f"0_{args.run_name}_pre_process"
    preprocess_dir.mkdir(parents=True, exist_ok=True)

    master_tsv = preprocess_dir / f"{args.run_name}_master_sumstats.tsv"
    preprocess_log = preprocess_dir / f"{args.run_name}_master_sumstats.log"

    master_vcf_tsv, success = run_merge_vcf_to_tsv_pipeline(
        input_df=val_input["manifest_df"],
        out_tsv=str(master_tsv),
        fixed_fields=["ID", "CHROM", "POS", "REF", "ALT"],
        format_fields=["AF", "ES", "SE", "NEF", "LP", "SI", "EZ"],
        log_file=str(preprocess_log),
        n_threads=int(args.cores),
        qc_n_rows=5_000_000,
        bcftools_exec="bcftools",
        fail_fast_on_merge=True,
    )

    if not success:
        raise RuntimeError("Master TSV creation failed.")

    # ===========================
    # STEP 3: Create PLACO input
    # ===========================
    print("\n[STEP 3] Creating PLACO input TSV...", flush=True)

    placo_input_dir = out_path / f"1_{args.run_name}_placo_input"
    placo_input_dir.mkdir(parents=True, exist_ok=True)

    placo_inputs = run_mastertsv_to_toolformat_pipeline(
        mode={"joined_samples"},
        tsv_file=master_vcf_tsv,
        input_manifest_df=val_input["manifest_df"],
        joined_out_tsv=str(placo_input_dir / f"{args.run_name}_placo_ready.tsv"),
        joined_remove_duplicate_ids=True,
        joined_create_unique_id=True,
        joined_sample_retained_cols=["ID", "LP", "EZ"],
        joined_header_map=None,
        joined_suffix_replace={":EZ": "_Z", ":LP": "_P"},
        joined_convert_to_raw_p_value=True,
        log_file=str(preprocess_log),
    )

    placo_input_file = placo_inputs["joined"]

    # ===========================
    # STEP 4: Run PLACO (R)
    # ===========================
    print("\n[STEP 4] Running PLACO (R engine)...", flush=True)

    placo_cli = str("/app/src/postgwas_pleio/meta_analysis/placo/R/placo_cli.R")
    placo_script = str("/opt/placo/PLACO_v0.2.0.R")

    placo_output = out_path / f"{args.run_name}_placo_results.csv"
    placo_log = out_path / f"{args.run_name}_placo.log"

    cmd_placo = [
        "micromamba", "run", "-n", "postgwas",
        "Rscript", placo_cli,
        "--input", str(placo_input_file),
        "--placo_script", placo_script,
        "--output", str(placo_output),
        "--method", str(args.method),
        "--id_col", "ID",
        "--z_suffix", "_Z",
        "--p_suffix", "_P",
        "--pthreshold", str(args.pthreshold),
    ]

    run_cmd_with_log(cmd_placo, placo_log, "PLACO")

    print("\n[SUCCESS] PLACO finished.", flush=True)

    return {
        "placo_input": str(placo_input_file),
        "placo_output": str(placo_output),
        "placo_log": str(placo_log),
    }