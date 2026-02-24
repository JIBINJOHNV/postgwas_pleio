import os
import subprocess
import datetime
import logging
from pathlib import Path

from postgwas_pleio.utils.input_manifest_val import validate_manifest
from postgwas_pleio.formatter.merge_vcf_to_tsv import run_merge_vcf_to_tsv_pipeline
from postgwas_pleio.formatter.master_tsv_to_tool_format import run_mastertsv_to_toolformat_pipeline

logger = logging.getLogger(__name__)


def mtag_pipeline_runner(args):
    """
    Validate manifest → merge VCFs → create MTAG inputs → run MTAG → rename outputs to sample IDs.
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
        raise RuntimeError("Master TSV creation failed. Cannot proceed to MTAG.")

    # ---------------------------
    # STEP 3: Create MTAG inputs
    # ---------------------------
    print("\n[STEP 3] Creating MTAG specific file format from master TSV...", flush=True)

    mtag_input_dir = out_path / f"1_{args.run_name}_mtag_input"
    mtag_input_dir.mkdir(parents=True, exist_ok=True)

    mtag_inputs = run_mastertsv_to_toolformat_pipeline(
        mode={"joined_samples", "sample_specific"},
        tsv_file=master_vcf_tsv,
        input_manifest_df=val_input["manifest_df"],
        joined_out_tsv=str(mtag_input_dir / f"{args.run_name}_mtag_ready.tsv"),
        sample_specific_directory=str(mtag_input_dir),
        joined_remove_duplicate_ids=True,
        joined_create_unique_id=False,
        joined_sample_retained_cols=["ID", "CHROM", "POS", "REF", "ALT", "AF", "NEF", "LP", "SI", "EZ"],
        joined_header_map=None,
        joined_suffix_replace=None,
        joined_convert_to_raw_p_value=True,
        sample_remove_duplicate_ids=True,
        sample_create_unique_id=False,
        sample_specific_retained_cols=["ID", "CHROM", "POS", "REF", "ALT", "AF", "NEF", "LP", "SI", "EZ"],
        sample_header_map=None,
        sample_suffix_replace=None,
        sample_convert_to_raw_p_value=True,
        log_file=str(preprocess_log_file),
    )

    # Use the *actual* order of sample-specific outputs (guarantees trait_1 ↔ sample_id match)
    sample_ids = mtag_inputs["sample_specific"]["sample_id"].to_list()

    # Ensure sumstats list are strings
    sumstats_list = [str(p) for p in mtag_inputs["sample_specific"]["matrix_tsv"].to_list()]
    sumstats_inputs = ",".join(sumstats_list)

    # ---------------------------
    # STEP 4: Run MTAG
    # ---------------------------
    print("\n[STEP 4] Running MTAG statistical engine...", flush=True)

    MTAG_ROOT = "/opt/mtag"
    mtag_output_dir = out_path / f"2_{args.run_name}_mtag_results"
    mtag_output_dir.mkdir(parents=True, exist_ok=True)

    output_prefix = str(mtag_output_dir / f"{args.run_name}_mtag_results")
    mtag_log_file = str(mtag_output_dir / f"{args.run_name}_mtag_results.log")

    cmd_mtag = [
        "micromamba", "run", "-n", "mtag",
        "python", f"{MTAG_ROOT}/mtag.py",
        "--sumstats", sumstats_inputs,
        "--out", output_prefix,
        "--snp_name", "ID",
        "--z_name", "EZ",
        "--beta_name", "ES",
        "--se_name", "SE",
        "--n_name", "NEF",
        "--eaf_name", "AF",
        "--chr_name", "CHROM",
        "--bpos_name", "POS",
        "--a1_name", "ALT",
        "--a2_name", "REF",
        "--p_name", "LP",
        "--cores", str(int(args.cores)),
        "--chunksize", str(int(args.chunksize)),
        "--stream_stdout",
        "--force",
    ]

    # Optional arguments
    if getattr(args, "ld_ref_panel", None):
        cmd_mtag.extend(["--ld_ref_panel", str(args.ld_ref_panel)])
    if getattr(args, "no_overlap", False):
        cmd_mtag.append("--no_overlap")
    if getattr(args, "perfect_gencov", False):
        cmd_mtag.append("--perfect_gencov")
    if getattr(args, "equal_h2", False):
        cmd_mtag.append("--equal_h2")
    if getattr(args, "std_betas", False):
        cmd_mtag.append("--std_betas")
    if getattr(args, "fdr", False):
        cmd_mtag.append("--fdr")

    print(f"DEBUG COMMAND: {' '.join(cmd_mtag)}\n", flush=True)

    try:
        subprocess.run(cmd_mtag, check=True)
        print("\n[SUCCESS] MTAG statistical engine finished. Proceeding to renaming...", flush=True)

        # ---------------------------
        # STEP 5: Rename MTAG outputs
        # ---------------------------
        print("[STEP 5] Mapping generic results to sample names...", flush=True)

        with open(mtag_log_file, "a") as f_log:
            f_log.write("\n" + "=" * 130 + "\n")
            f_log.write(f"MTAG OUTPUT RENAMING - {datetime.datetime.now()}\n")
            f_log.write("=" * 130 + "\n")

            for i, sample_id in enumerate(sample_ids):
                trait_idx = i + 1
                generic_path = f"{output_prefix}_trait_{trait_idx}.txt"

                final_filename = f"{args.run_name}_mtag_results_trait_{trait_idx}_{sample_id}.txt"
                final_full_path = mtag_output_dir / final_filename

                if os.path.exists(generic_path):
                    os.rename(generic_path, str(final_full_path))
                    f_log.write(f"trait_{trait_idx} → {sample_id} → {final_full_path}\n")
                    print(f"   - Mapped trait_{trait_idx} to {sample_id}", flush=True)
                else:
                    f_log.write(f"trait_{trait_idx} → {sample_id} → NOT FOUND ({generic_path})\n")
                    print(f"   - [WARNING] trait_{trait_idx} not found for {sample_id}", flush=True)

            f_log.write("=" * 130 + "\n")

        print(f"\n[FINISH] MTAG workflow completed. Audit trail: {mtag_log_file}", flush=True)

    except subprocess.CalledProcessError as e:
        print(f"\n[ERROR] MTAG failed with exit code {e.returncode}", flush=True)
        raise RuntimeError("MTAG statistical analysis failed.") from e

    return {
        "sample_ids": sample_ids,
        "mtag_input_files": sumstats_list,
        "mtag_results_dir": str(mtag_output_dir),
        "log_file": mtag_log_file,
    }