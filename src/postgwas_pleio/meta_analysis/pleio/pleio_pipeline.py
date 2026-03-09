import subprocess
import logging
import os
from pathlib import Path

from postgwas_pleio.utils.input_manifest_val import validate_manifest
from postgwas_pleio.formatter.merge_vcf_to_tsv import run_merge_vcf_to_tsv_pipeline
from postgwas_pleio.formatter.master_tsv_to_tool_format import run_mastertsv_to_toolformat_pipeline

logger = logging.getLogger(__name__)


def run_cmd_with_log(cmd, log_file, step_name):
    """
    Run subprocess command with tee-style logging (console + file).
    Accepts Path items in cmd/log_file; internally converts to str.
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

        for line in process.stdout:
            print(line, end="")     # console
            log.write(line)         # file

        process.wait()

        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode, cmd)


def pleio_pipeline_runner(args):
    """
    Validate manifest → merge VCFs → create PLEIO inputs → LDSC preprocess → Pleio meta.
    """

    # ---------------------------
    # STEP 1: Validate manifest
    # ---------------------------
    out_path = Path(args.out)
    out_path.mkdir(parents=True, exist_ok=True)

    print("\n[STEP 1] Validating manifest...", flush=True)
    val_input = validate_manifest(args.inputfile, str(out_path))

    # ---------------------------
    # STEP 2: Merge VCFs → master TSV
    # ---------------------------
    print("\n[STEP 2] Merging individual VCFs into a master TSV...", flush=True)

    pre_process_dir = out_path / f"0_{args.run_name}_pre_process"
    pre_process_dir.mkdir(parents=True, exist_ok=True)

    out_tsv = pre_process_dir / f"{args.run_name}_master_sumstats.tsv"
    preprocess_log_file = pre_process_dir / f"{args.run_name}_master_sumstats.log"

    master_vcf_tsv, success_2 = run_merge_vcf_to_tsv_pipeline(
        input_df=val_input["manifest_df"],
        out_tsv=str(out_tsv),
        fixed_fields=["ID", "CHROM", "POS", "REF", "ALT"],
        format_fields=["AF", "ES", "SE", "NEF", "LP", "SI", "EZ"],
        log_file=str(preprocess_log_file),
        n_threads=int(args.cores),
        qc_n_rows=5_000_000,
        bcftools_exec="bcftools",
        fail_fast_on_merge=True,
    )

    if not success_2:
        raise RuntimeError("Master TSV creation failed. Cannot proceed to PLEIO.")

    # ---------------------------
    # STEP 3: Create PLEIO inputs
    # ---------------------------
    print("\n[STEP 3] Creating PLEIO specific file format from master TSV...", flush=True)

    pleio_input_dir = out_path / f"1_{args.run_name}_pleio_input"
    pleio_input_dir.mkdir(parents=True, exist_ok=True)

    pleio_inputs = run_mastertsv_to_toolformat_pipeline(
        mode={"joined_samples", "sample_specific"},
        tsv_file=master_vcf_tsv,
        input_manifest_df=val_input["manifest_df"],
        joined_out_tsv=str(pleio_input_dir / f"{args.run_name}_pleio_ready.tsv"),
        sample_specific_directory=str(pleio_input_dir),
        joined_remove_duplicate_ids=True,
        joined_create_unique_id=False,
        joined_sample_retained_cols=["ID", "CHROM", "POS", "REF", "ALT", "AF", "NEF", "ES", "SE", "LP", "SI", "EZ"],
        joined_header_map=None,
        joined_suffix_replace=None,
        joined_convert_to_raw_p_value=True,
        sample_remove_duplicate_ids=True,
        sample_create_unique_id=False,
        sample_specific_retained_cols=["ID", "CHROM", "POS", "REF", "ALT", "AF", "NEF", "ES", "SE", "LP", "SI", "EZ"],
        sample_header_map={"REF": "A2", "ALT": "A1", "ID": "SNP"},
        sample_suffix_replace={"ES": "BETA", "EZ": "Z", "AF": "FRQ", "NEF": "N", "SI": "INFO"},
        sample_convert_to_raw_p_value=True,
        log_file=str(preprocess_log_file),
    )

    # Create pleio input file (NAME + FILE)
    pleio_input_df = pleio_inputs["sample_specific"].rename({"matrix_tsv": "FILE", "sample_id": "NAME"})
    pleio_input_file = pleio_input_dir / f"{args.run_name}_pleio_input_file.tsv"
    pleio_input_df.write_csv(pleio_input_file, separator="\t")

    sample_ids = pleio_input_df["NAME"].to_list()

    PLEIO_ROOT = "/opt/pleio"

    # ---------------------------
    # STEP 4: LDSC preprocess
    # ---------------------------
    ldsc_out_dir = out_path / "2_pleio_ldsc"
    ldsc_out_dir.mkdir(parents=True, exist_ok=True)

    ldsc_log = out_path / f"{args.run_name}_ldsc_preprocess.log"

    cmd_pleio_ldsc = [
        "micromamba", "run", "-n", "mtag",
        "python", os.path.join(PLEIO_ROOT, "ldsc_preprocess.py"),
        "--input", str(pleio_input_file),
        "--ref-ld-chr", str(args.ld_ref_panel),
        "--w-ld-chr", str(args.ld_ref_panel),
        "--out", str(ldsc_out_dir),
    ]

    print("[STEP 4] Starting Pleio ldsc_preprocess...", flush=True)
    run_cmd_with_log(cmd_pleio_ldsc, ldsc_log, "LDSC PREPROCESS")
    print("\n[SUCCESS] Pleio ldsc_preprocess finished.", flush=True)

    # ---------------------------
    # STEP 5: Pleio meta
    # ---------------------------
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

    if getattr(args, "flattening_p_values", False):
        cmd_pleio_meta.append("--flattening_p_values")

    cmd_pleio_meta.extend(["--out", str(out_path / f"3_{args.run_name}_pleio_results")])

    print("[STEP 5] Starting Pleio meta-analysis...", flush=True)
    run_cmd_with_log(cmd_pleio_meta, pleio_log, "PLEIO META")
    print("\n[SUCCESS] Pleio meta-analysis finished.", flush=True)

    return {
        "sample_ids": sample_ids,
        "pleio_input_file": str(pleio_input_file),
        "ldsc_out_dir": str(ldsc_out_dir),
        "pleio_results_dir": str(out_path / f"{args.run_name}_pleio_results"),
        "preprocess_log_file": str(preprocess_log_file),
        "ldsc_log": str(ldsc_log),
        "pleio_log": str(pleio_log),
    }