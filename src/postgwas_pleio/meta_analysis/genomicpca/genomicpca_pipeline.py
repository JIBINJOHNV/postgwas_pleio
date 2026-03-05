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


def genomicpca_pipeline_runner(args):
    """
    Validate manifest → merge VCFs → create GenomicPCA inputs → run GenomicPCA (R) → return result paths.
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
        raise RuntimeError("Master TSV creation failed. Cannot proceed to GenomicPCA.")
    
    # ---------------------------
    # STEP 3: Create GenomicPCA inputs
    # ---------------------------
    print("\n[STEP 3] Creating GenomicPCA specific inputs from master TSV...", flush=True)
    asset_input_dir = out_path / f"1_{args.run_name}_asset_input"
    ldsc_input_dir = out_path / f"2_{args.run_name}_ldsc_input"
    asset_input_dir.mkdir(parents=True, exist_ok=True)
    ldsc_input_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate sample-specific files used by the GenomicPCA pipeline
    asset_inputs = run_mastertsv_to_toolformat_pipeline(
        mode={"joined_samples", "sample_specific", "sample_specific_ldsc"},
        tsv_file=master_vcf_tsv,
        input_manifest_df=val_input["manifest_df"],
        joined_out_tsv=str(asset_input_dir / f"{args.run_name}_asset_joined.tsv"),
        sample_specific_directory=str(asset_input_dir),
        sample_specific_ldsc_directory=str(ldsc_input_dir),
        
        joined_remove_duplicate_ids=True,
        joined_create_unique_id=True,
        joined_sample_retained_cols=["ID", "NEF", "ES", "SE"],
        joined_header_map=None,
        joined_suffix_replace={":ES": ".Beta", ":SE": ".SE", ":NEF": ".N"},
        joined_convert_to_raw_p_value=False,
        
        sample_remove_duplicate_ids=True,
        sample_create_unique_id=True,
        sample_specific_retained_cols=["ID", "CHROM", "POS", "ALT", "REF", "AF", "NEF", "EZ", "LP"],
        sample_header_map={"ID": "SNPID", "CHROM": "CHR", "POS": "BP", "ALT": "EA", "REF": "OA"},
        sample_suffix_replace={"AF": "EAF", "NEF": "N", "EZ": "Z", "LP": "P"},
        sample_convert_to_raw_p_value=True,
        
        sample_ldsc_remove_duplicate_ids=False,
        sample_ldsc_create_unique_id=False,
        sample_ldsc_specific_retained_cols=["ID", "REF", "ALT", "AF", "NEF", "ES", "LP", "SI"],
        sample_ldsc_header_map={"ID": "SNPID", "REF": "A2", "ALT": "A1"},
        sample_ldsc_suffix_replace={"NEF": "N", "AF": "MAF", "ES": "effect", "LP": "P", "SI": "INFO"},
        sample__ldscconvert_to_raw_p_value=True,
        log_file=str(preprocess_log_file),
    )
    
    # Build GenomicPCA input manifest
    asset_input_df = asset_inputs["sample_specific"].rename(
        {"matrix_tsv": "FILE", "sample_id": "NAME"}
    )
    
    asset_ldsc_input_df = asset_inputs["sample_specific_ldsc"].rename(
        {"ldsc_matrix_tsv": "ldsc_FILE", "sample_id": "NAME"}
    ).select(["sumstat_vcf", "NAME", "ldsc_FILE"])
    
    asset_input_df = asset_input_df.join(
        asset_ldsc_input_df,
        on=["sumstat_vcf", "NAME"],
        how="left"
    )
    
    # Ensure PPREV/SPREV fields exist if required
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
    
    # ---------------------------
    # STEP 4: Run GenomicPCA (R)
    # ---------------------------
    print("\n[STEP 4] Running GenomicPCA pipeline (R)...", flush=True)
    r_script_path = files("postgwas_pleio.meta_analysis.genomicpca") / "R" / "genomicpca_cli.R"

    if not r_script_path.is_file():
        raise FileNotFoundError(
            f"[ERROR] GenomicPCA R script not found:\n{r_script_path}\n"
            "Check package installation or MANIFEST configuration."
        )
    genomicpca_cli = str(r_script_path)
    genomicpca_log = out_path / f"{args.run_name}_genomicpca.log"
    
    cmd_genomicpca = [
        "micromamba", "run", "-n", "postgwas",
        "Rscript", genomicpca_cli,
        "--input_manifest", str(asset_input_file),
        "--output_dir", str(out_path),
        "--hm3", str(args.hm3),
        "--ld_ref", str(args.ld_ref),
        "--info_filter", str(args.info_filter),
        "--maf_filter", str(args.maf_filter),
        "--run_name", str(args.run_name),
    ]
    
    if hasattr(args, "approach"):
        cmd_genomicpca.extend(["--approach", str(args.approach)])
    
    run_cmd_with_log(cmd_genomicpca, genomicpca_log, "GenomicPCA")
    print("\n[SUCCESS] GenomicPCA pipeline finished.", flush=True)
    
    # ---------------------------
    # STEP 5: Return metadata
    # ---------------------------
    gpca_corr_dir = out_path / "4_GenomicPCA_correlation"
    gpca_cov_dir = out_path / "5_GenomicPCA_alternative_covariates"
    
    gpca_corr_results = gpca_corr_dir / f"{args.run_name}_GenomicPCA_correlation_meta.tsv"
    gpca_cov_results = gpca_cov_dir / f"{args.run_name}_GenomicPCA_covariates_meta.tsv"
    return {
        "sample_ids": sample_ids,
        "genomicpca_input_manifest": str(asset_input_file),
        "genomicpca_correlation_results": str(gpca_corr_results),
        "genomicpca_covariance_results": str(gpca_cov_results),
        "preprocess_log_file": str(preprocess_log_file),
        "genomicpca_log": str(genomicpca_log),
    }