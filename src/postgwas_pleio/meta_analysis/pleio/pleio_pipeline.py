import subprocess
import logging
import os
from pathlib import Path
import polars as pl
from typing import List

from postgwas_pleio.utils.input_manifest_val import validate_manifest
from postgwas_pleio.formatter.merge_vcf_to_tsv import run_merge_vcf_to_tsv_pipeline
from postgwas_pleio.formatter.master_tsv_to_tool_format import run_mastertsv_to_toolformat_pipeline

logger = logging.getLogger(__name__)




def run_cmd_with_log(cmd, log_file, step_name):
    cmd = [str(x) for x in cmd]

    log_file = Path(log_file)
    log_file.parent.mkdir(parents=True, exist_ok=True)   # ✅ FIX
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

        if process.stdout is None:
            raise RuntimeError("Failed to capture subprocess output")

        for line in process.stdout:
            print(line, end="")
            log.write(line)

        process.wait()

        if process.returncode != 0:
            log.write(f"\n[ERROR] Return code: {process.returncode}\n")
            raise subprocess.CalledProcessError(process.returncode, cmd)



def _run_harmonisation(args, out_path: Path, manifest_file_path: str, log_file: Path) -> None:
    resource_folder = Path(args.resource_folder).resolve()
    defaults_path = Path(args.defaults_config).resolve()
    manifest_path = Path(manifest_file_path).resolve()
    out_path = Path(out_path).resolve()

    all_paths = [
        str(resource_folder),
        str(defaults_path.parent),
        str(manifest_path.parent),
        str(out_path)
    ]

    common_root = os.path.commonpath(all_paths)

    if common_root == "/":
        print("[WARNING] Common root is '/', using individual mounts instead.")
        mounts = []
        for p in set(all_paths):
            mounts.extend(["-v", f"{p}:{p}"])
    else:
        print(f"[DEBUG] Smart Mount Root: {common_root}")
        mounts = ["-v", f"{common_root}:{common_root}"]

    cmd = [
        "docker", "run", "--rm", "--platform", "linux/amd64",
        "-e", "PYTHONUNBUFFERED=1",
        *mounts,
        "jibinjv/postgwas:1.3",
        "postgwas", "harmonisation",
        "--nthreads", str(args.nthreads),
        "--max-mem", str(args.max_mem),
        "--config", str(manifest_path),
        "--defaults", str(defaults_path),
    ]

    print("[DEBUG] Docker CMD:", " ".join(cmd), flush=True)

    run_cmd_with_log(
        cmd=cmd,
        log_file=log_file,
        step_name="POSTGWAS HARMONISATION"
    )



def generate_pleio_harmonisation_inputs(
    pleio_input_df: pl.DataFrame,
    p_meta_file: str,
    blup_file: str,
    harmonisation_input_dir: str,
    harmonised_output_dir: str,
    pvalue_column: str,
    sample_id_suffix: str,
    args,
) -> pl.DataFrame:
    """
    Generate harmonisation input files using Polars with trait-type aware logic.

    Binary traits:
        Requires N (sample size) and NC (number of cases)
        Computes NCONTROL = N - NC

    Quantitative traits:
        Uses N only

    Returns:
        manifest_df (pl.DataFrame)
    """

    # ------------------------------------------------------------------
    # Setup paths
    # ------------------------------------------------------------------
    harmonisation_input_dir = Path(harmonisation_input_dir)
    harmonised_output_dir = Path(harmonised_output_dir)

    harmonisation_input_dir.mkdir(parents=True, exist_ok=True)
    harmonised_output_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Validate input
    # ------------------------------------------------------------------
    if "TYPE" not in pleio_input_df.columns:
        raise ValueError("pleio_input_df must contain 'TYPE' column")

    # ------------------------------------------------------------------
    # Read core inputs
    # ------------------------------------------------------------------
    p_meta_df = (
        pl.read_csv(p_meta_file, separator="\t")
        .select(["SNP", pvalue_column])
        .with_columns([
            pl.col(pvalue_column).cast(pl.Float64, strict=False)
        ])
        .rename({pvalue_column: "P"})
    )

    blup_df = pl.read_csv(blup_file, separator="\t")

    sample_ids: List[str] = pleio_input_df["NAME"].to_list()
    manifest_rows = []

    # ------------------------------------------------------------------
    # Loop per sample
    # ------------------------------------------------------------------
    for sample_id in sample_ids:
        safe_sample_id = str(sample_id).replace("/", "_").replace(" ", "_")
        # ---------------------------
        # Extract metadata
        # ---------------------------
        row = pleio_input_df.filter(pl.col("NAME") == sample_id)

        input_sample_file = row.select("FILE").item()
        trait_type = row.select("TYPE").item()

        input_sample_df = pl.read_csv(input_sample_file, separator="\t")

        # ---------------------------
        # BLUP columns
        # ---------------------------
        beta_col = sample_id
        se_col = f"{sample_id}_se"

        if beta_col not in blup_df.columns or se_col not in blup_df.columns:
            raise ValueError(f"Missing BLUP columns for {sample_id}")

        sample_blup_df = (
            blup_df
            .select(["SNP", beta_col, se_col])
            .with_columns([
                pl.col(beta_col).cast(pl.Float64, strict=False),
                pl.col(se_col).cast(pl.Float64, strict=False),
            ])
        )

        # ---------------------------
        # Merge BLUP + p_meta
        # ---------------------------
        blup_meta_df = (
            p_meta_df
            .join(sample_blup_df, on="SNP", how="inner")
            .rename({
                beta_col: "Beta",
                se_col: "SE"
            })
        )

        # ---------------------------
        # Trait-specific handling
        # ---------------------------
        base_cols = ['CHROM', 'POS', 'SNP', 'A2', 'A1', 'FRQ', 'INFO']

        if trait_type == "binary":
            required_cols = base_cols + ["N", "NC"]

            missing = [c for c in required_cols if c not in input_sample_df.columns]
            if missing:
                raise ValueError(f"{sample_id} missing columns: {missing}")

            input_sample_df = input_sample_df.select(required_cols)

            # 🔥 CRITICAL: Cast before math
            input_sample_df = input_sample_df.with_columns([
                pl.col("N").cast(pl.Float64, strict=False),
                pl.col("NC").cast(pl.Float64, strict=False),
                pl.col("FRQ").cast(pl.Float64, strict=False),
                pl.col("INFO").cast(pl.Float64, strict=False),
            ])

            input_sample_df = input_sample_df.with_columns(
                pl.when(pl.col("N") >= pl.col("NC"))
                .then(pl.col("N") - pl.col("NC"))
                .otherwise(None)
                .alias("NCONTROL")
            )

            ncase_col = "NC"
            ncontrol_col = "NCONTROL"

        else:
            required_cols = base_cols + ["N"]

            missing = [c for c in required_cols if c not in input_sample_df.columns]
            if missing:
                raise ValueError(f"{sample_id} missing columns: {missing}")

            input_sample_df = input_sample_df.select(required_cols)

            input_sample_df = input_sample_df.with_columns([
                pl.col("N").cast(pl.Float64, strict=False),
                pl.col("FRQ").cast(pl.Float64, strict=False),
                pl.col("INFO").cast(pl.Float64, strict=False),
            ])

            ncase_col = "NA"
            ncontrol_col = "N"

        # ---------------------------
        # Merge with BLUP/meta
        # ---------------------------
        harmonisation_df = input_sample_df.join(blup_meta_df, on="SNP", how="inner")

        # 🔥 FINAL TYPE SAFETY (MOST IMPORTANT FIX)
        numeric_cols = ["Beta", "SE", "P", "FRQ", "INFO", "N"]

        if trait_type == "binary":
            numeric_cols += ["NC", "NCONTROL"]

        harmonisation_df = harmonisation_df.with_columns([
            pl.col(c).cast(pl.Float64, strict=False)
            for c in numeric_cols
            if c in harmonisation_df.columns
        ])

        # ---------------------------
        # Output file
        # ---------------------------
        out_file = harmonisation_input_dir / f"{safe_sample_id}_{sample_id_suffix}_meta_harmonisation_input.tsv" 
        harmonisation_df.write_csv(out_file, separator="\t")

        # ---------------------------
        # Manifest entry
        # ---------------------------
        sample_out_dir = harmonised_output_dir / safe_sample_id
        sample_out_dir.mkdir(parents=True, exist_ok=True)

        manifest_rows.append({
            "sumstat_file": str(out_file),
            "gwas_outputname": f"{safe_sample_id}_{sample_id_suffix}", 
            "chr_col": "CHROM",
            "pos_col": "POS",
            "snp_id_col": "SNP",
            "ea_col": "A1",
            "oa_col": "A2",
            "eaf_col": "FRQ",
            "beta_or_col": "Beta",
            "se_col": "SE",
            "imp_z_col": "NA",
            "pval_col": "P",
            "ncontrol_col": ncontrol_col,
            "ncase_col": ncase_col,
            "ncontrol": "NA",
            "ncase": "NA",
            "imp_info_col": "INFO",
            "infofile": "NA",
            "infocolumn": "NA",
            "eaffile": "NA",
            "eafcolumn": "NA",
            "liftover": "Yes",
            "chr_pos_col": "NA",
            "resourse_folder": str(args.resource_folder),
            "output_folder": str(sample_out_dir),
        })

        print(
            f"[INFO] {safe_sample_id} → TYPE={trait_type}, "
            f"ncase_col={ncase_col}, ncontrol_col={ncontrol_col}",
            flush=True
        )

        print(f"[DEBUG] Dtypes: {harmonisation_df.dtypes}", flush=True)

    # ------------------------------------------------------------------
    # Final manifest
    # ------------------------------------------------------------------
    manifest_df = pl.DataFrame(manifest_rows)

    print("\n[INFO] Final manifest preview:", flush=True)
    print(manifest_df, flush=True)

    return manifest_df


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
        format_fields=["AF", "ES", "SE", "NEF", "LP", "SI", "EZ", "SS", "NC"],
        log_file=str(preprocess_log_file),
        n_threads=int(args.nthreads),
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
        sample_specific_retained_cols=["ID", "CHROM", "POS", "REF", "ALT", "AF", "SS", "ES", "SE", "LP", "SI", "EZ",'NC'],
        sample_header_map={"REF": "A2", "ALT": "A1", "ID": "SNP"},
        sample_suffix_replace={"ES": "BETA", "EZ": "Z", "AF": "FRQ", "SS": "N", "SI": "INFO"},
        sample_calculate_sample_prevalence=True,
        sample_total_sample_size_col="SS",
        sample_number_of_cases_col="NC",
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
    pleio_analysis_dir = out_path / "2_pleio_output"
    pleio_analysis_dir.mkdir(parents=True, exist_ok=True) 

    ldsc_out_dir = pleio_analysis_dir / "pleio_ldsc/"
    pleio_out_dir = pleio_analysis_dir / "pleio_output/"
    pleio_output_file=str(f"{pleio_out_dir}/{args.run_name}_pleio_output")

    ldsc_out_dir.mkdir(parents=True, exist_ok=True)
    pleio_out_dir.mkdir(parents=True, exist_ok=True)

    ldsc_log = ldsc_out_dir / f"{args.run_name}_ldsc_preprocess.log"

    harmonisation_input_dir=out_path / "4_harmonisation_input_files"
    harmonisation_input_dir.mkdir(parents=True, exist_ok=True)  

    harmonised_output_dir=out_path / "5_harmonised_sumstats"
    harmonised_output_dir.mkdir(parents=True, exist_ok=True)

    cmd_pleio_ldsc = [
        "micromamba", "run", "-n", "mtag",
        "python", str(Path(PLEIO_ROOT) / "ldsc_preprocess.py"),
        "--input", str(pleio_input_file),
        "--ref-ld-chr", str(args.ld_ref_panel),
        "--w-ld-chr", str(args.ld_ref_panel),
        "--out", str(ldsc_out_dir),
    ]

    print("[STEP 4] Starting Pleio ldsc_preprocess...", flush=True)
    run_cmd_with_log(cmd_pleio_ldsc, ldsc_log, "LDSC PREPROCESS")
    print("\n[SUCCESS] Pleio ldsc_preprocess finished.", flush=True)

    # ✅ Validate LDSC outputs
    required_files = ["metain.txt.gz", "sg.txt.gz", "ce.txt.gz"]
    for f in required_files:
        fp = ldsc_out_dir / f
        if not fp.exists():
            raise FileNotFoundError(f"LDSC output missing: {fp}")

    # ---------------------------
    # STEP 5: Pleio meta 
    # ---------------------------
    pleio_log = pleio_out_dir / f"{args.run_name}_pleio_meta.log"

    cmd_pleio_meta = [
        "micromamba", "run", "-n", "pleio",
        "python", str(Path(PLEIO_ROOT) / "pleio.py"),
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

    cmd_pleio_meta.extend(["--out", str(pleio_output_file)])

    print("[STEP 5] Starting Pleio meta-analysis...", flush=True) 
    run_cmd_with_log(cmd_pleio_meta, pleio_log, "PLEIO META")
    print("\n[SUCCESS] Pleio meta-analysis finished.", flush=True)

    # ========================= 
    # ADD THIS BEFORE STEP 6 CALL
    # =========================

    harmonisation_log = None
    
    if args.harmonise:
        harmonisation_log = harmonised_output_dir / f"{args.run_name}_pleio_harmonisation.log"

        # ✅ FIX: Generate harmonisation inputs + manifest
        manifest_df_1 = generate_pleio_harmonisation_inputs(
            pleio_input_df=pleio_input_df,
            p_meta_file=pleio_output_file + ".txt.gz",
            blup_file=pleio_output_file + ".blup.gz",  
            harmonisation_input_dir=str(harmonisation_input_dir),
            harmonised_output_dir=str(harmonised_output_dir),
            args=args,
            pvalue_column="pleio_p",
            sample_id_suffix="pleio",
        )

        # ✅ FIX: Generate harmonisation inputs + manifest
        manifest_df_2 = generate_pleio_harmonisation_inputs(
            pleio_input_df=pleio_input_df,
            p_meta_file=pleio_output_file + ".txt.gz",
            blup_file=pleio_output_file + ".blup.gz",  
            harmonisation_input_dir=str(harmonisation_input_dir),
            harmonised_output_dir=str(harmonised_output_dir),
            args=args,
            pvalue_column="LS_p",
            sample_id_suffix="pleio_lsp",
        )

        manifest_df = pl.concat([manifest_df_1, manifest_df_2])
        # ✅ FIX: Save manifest and define path
        manifest_file_path = harmonisation_input_dir / f"{args.run_name}_harmonisation_manifest.tsv"
        manifest_df.write_csv(manifest_file_path, separator="\t")

        try:
            _run_harmonisation(
                args=args,
                out_path=harmonised_output_dir,
                manifest_file_path=str(manifest_file_path),   # ✅ FIX
                log_file=harmonisation_log,
            )
        except Exception as e:
            raise RuntimeError(f"Harmonisation failed. Check log: {harmonisation_log}") from e

    return {
        "sample_ids": sample_ids,
        "pleio_input_file": str(pleio_input_file),
        "ldsc_out_dir": str(ldsc_out_dir),
        "pleio_results_dir": str(pleio_out_dir),
        "preprocess_log_file": str(preprocess_log_file),
        "ldsc_log": str(ldsc_log),
        "pleio_log": str(pleio_log),
        "harmonisation_log": str(harmonisation_log),
    }