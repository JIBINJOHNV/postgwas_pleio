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


def asset_pipeline_runner(args):
    """
    Validate manifest → merge VCFs → create ASSET inputs → run ASSET (R) → return result paths.
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
        raise RuntimeError("Master TSV creation failed. Cannot proceed to ASSET.")

    # ---------------------------
    # STEP 3: Create ASSET inputs
    # ---------------------------
    print("\n[STEP 3] Creating ASSET specific inputs from master TSV...", flush=True)

    asset_input_dir = out_path / f"1_{args.run_name}_asset_input"
    asset_input_dir.mkdir(parents=True, exist_ok=True)

    # We need:
    # - joined file for fastASSET screening (joined)
    # - sample-specific files for "input_manifest" (NAME + FILE + optional fields)
    asset_inputs = run_mastertsv_to_toolformat_pipeline(
        mode={"joined_samples", "sample_specific"},
        tsv_file=master_vcf_tsv,
        input_manifest_df=val_input["manifest_df"],
        joined_out_tsv=str(asset_input_dir / f"{args.run_name}_asset_joined.tsv"),
        sample_specific_directory=str(asset_input_dir),

        joined_remove_duplicate_ids=True,
        joined_create_unique_id=True,
        joined_sample_retained_cols=["ID", "NEF", "ES", "SE"],
        joined_header_map=None,
        joined_suffix_replace={":ES": ".Beta", ":SE": ".SE", ":NEF": ".N"},
        joined_convert_to_raw_p_value=False,

        sample_remove_duplicate_ids=True,
        sample_create_unique_id=False,
        sample_specific_retained_cols=["ID", "CHROM", "POS", "REF", "ALT", "AF", "NEF", "ES", "SE", "LP", "SI", "EZ"],
        sample_header_map={"ID": "SNP", "CHROM": "CHR", "POS": "POS", "REF": "A2", "ALT": "A1"},
        sample_suffix_replace={"ES": "effect", "EZ": "Z", "AF": "MAF", "NEF": "N", "SI": "INFO", "LP": "P"},
        sample_convert_to_raw_p_value=True,

        log_file=str(preprocess_log_file),
    )

    # Build ASSET input manifest (NAME + FILE) from sample_specific table
    asset_input_df = asset_inputs["sample_specific"].rename({"matrix_tsv": "FILE", "sample_id": "NAME"})

    # IMPORTANT: Polars does not support pandas-style assignment.
    # If your R script requires PPREV/SPREV columns, they must already exist in asset_input_df,
    # or you should add them using Polars expressions. Here we only fill if columns exist.
    if "PPREV" in asset_input_df.columns:
        asset_input_df = asset_input_df.with_columns(
            asset_input_df["PPREV"].fill_null("nan").alias("PPREV")
        )
    if "SPREV" in asset_input_df.columns:
        asset_input_df = asset_input_df.with_columns(
            asset_input_df["SPREV"].fill_null("nan").alias("SPREV")
        )

    asset_input_file = asset_input_dir / f"{args.run_name}_asset_input_manifest.tsv"
    asset_input_df.write_csv(asset_input_file, separator="\t")

    sample_ids = asset_input_df["NAME"].to_list()

    # joined fastASSET input file path (depends on your toolformat pipeline return schema)
    joined_fastasset_input = asset_inputs.get("joined")
    # Some implementations return a Path; some return a string.
    joined_fastasset_input = str(joined_fastasset_input) if joined_fastasset_input is not None else None
    if not joined_fastasset_input:
        raise RuntimeError("ASSET joined input was not produced (asset_inputs['joined'] is empty).")

    # ---------------------------
    # STEP 4: Run ASSET (R)
    # ---------------------------
    print("\n[STEP 4] Running ASSET pipeline (R)...", flush=True)

    asset_cli = str(files("postgwas_pleio.meta_analysis.asset.R") / "asset_cli.R")
    asset_log = out_path / f"{args.run_name}_asset.log"

    cmd_asset = [
        "micromamba", "run", "-n", "postgwas",
        "Rscript", asset_cli,
        "--input_manifest", str(asset_input_file),
        "--fastasset_input", joined_fastasset_input,
        "--output_dir", str(out_path),
        "--hm3", str(args.hm3),
        "--ld_ref", str(args.ld_ref_panel),
        "--info_filter", str(args.info_filter),
        "--maf_filter", str(args.maf_filter),
        "--chunk_size", str(args.chunk_size),
        "--scr_pthr", str(args.scr_pthr),
        "--meth_pval", str(args.meth_pval),
        "--ncores", str(args.cores),
        "--run_name", str(args.run_name),
    ]

    run_cmd_with_log(cmd_asset, asset_log, "ASSET")
    print("\n[SUCCESS] ASSET finished.", flush=True)

    # ---------------------------
    # STEP 5: Return metadata
    # ---------------------------
    fastasset_1sided = out_path / "5_fastasset_output" / f"{args.run_name}_fastasset_1sided.tsv.gz"
    fastasset_2sided = out_path / "5_fastasset_output" / f"{args.run_name}_fastasset_2sided.tsv.gz"

    asset_1sided = out_path / "6_htraits_output" / f"{args.run_name}_htraits_1sided.tsv.gz"
    asset_2sided = out_path / "6_htraits_output" / f"{args.run_name}_htraits_2sided.tsv.gz"
    asset_meta = out_path / "6_htraits_output" / f"{args.run_name}_htraits_meta.tsv.gz"

    return {
        "sample_ids": sample_ids,
        "asset_input_manifest": str(asset_input_file),
        "fastasset_input": str(joined_fastasset_input),
        "fastasset_1sided": str(fastasset_1sided),
        "fastasset_2sided": str(fastasset_2sided),
        "asset_1sided": str(asset_1sided),
        "asset_2sided": str(asset_2sided),
        "asset_meta": str(asset_meta),
        "preprocess_log_file": str(preprocess_log_file),
        "asset_log": str(asset_log),
    }