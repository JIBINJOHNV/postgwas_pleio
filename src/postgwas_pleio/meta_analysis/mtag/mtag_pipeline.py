import os
import subprocess
import datetime
import logging
from pathlib import Path
import polars as pl

from postgwas_pleio.utils.input_manifest_val import validate_manifest
from postgwas_pleio.formatter.merge_vcf_to_tsv import run_merge_vcf_to_tsv_pipeline
from postgwas_pleio.formatter.master_tsv_to_tool_format import run_mastertsv_to_toolformat_pipeline
from postgwas_pleio.utils.harmonise_input import create_manifest_with_info

logger = logging.getLogger(__name__)


def _append_log(log_file: Path, message: str) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    with log_file.open("a", encoding="utf-8") as f:
        f.write(message.rstrip() + "\n")


def _run_subprocess_with_log(cmd, log_file: Path, step_name: str, cwd: Path | None = None):

    log_file.parent.mkdir(parents=True, exist_ok=True)

    with log_file.open("a", encoding="utf-8") as log:

        log.write("\n" + "=" * 130 + "\n")
        log.write(f"{step_name} - {datetime.datetime.now()}\n")
        log.write("=" * 130 + "\n")
        log.write("COMMAND:\n")
        log.write(" ".join(map(str, cmd)) + "\n\n")
        log.flush()

        try:
            process = subprocess.Popen(
                cmd,
                cwd=str(cwd) if cwd else None,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
            )

            if process.stdout is None:
                raise RuntimeError("Failed to capture subprocess stdout")

            # 🔁 STREAM
            for line in process.stdout:
                print(line, end="")
                log.write(line)

            process.stdout.close()
            returncode = process.wait()

            if returncode != 0:
                log.write(f"\n[ERROR] {step_name} failed with exit code {returncode}\n")
                log.flush()
                raise RuntimeError(f"{step_name} failed with exit code {returncode}")

            log.write(f"\n[SUCCESS] {step_name} completed successfully\n")
            log.flush()

            return subprocess.CompletedProcess(cmd, returncode, None, None)

        except Exception as e:
            log.write(f"\n[EXCEPTION] {str(e)}\n")
            log.flush()
            raise

def _build_empty_manifest_df() -> pl.DataFrame:
    return pl.DataFrame(schema={
        "sumstat_file": pl.Utf8,
        "gwas_outputname": pl.Utf8,
        "chr_col": pl.Utf8,
        "pos_col": pl.Utf8,
        "snp_id_col": pl.Utf8,
        "ea_col": pl.Utf8,
        "oa_col": pl.Utf8,
        "eaf_col": pl.Utf8,
        "beta_or_col": pl.Utf8,
        "se_col": pl.Utf8,
        "imp_z_col": pl.Utf8,
        "pval_col": pl.Utf8,
        "ncontrol_col": pl.Utf8,
        "ncase_col": pl.Utf8,
        "ncontrol": pl.Utf8,
        "ncase": pl.Utf8,
        "imp_info_col": pl.Utf8,
        "infofile": pl.Utf8,
        "infocolumn": pl.Utf8,
        "eaffile": pl.Utf8,
        "eafcolumn": pl.Utf8,
        "liftover": pl.Utf8,
        "chr_pos_col": pl.Utf8,
        "resourse_folder": pl.Utf8,
        "output_folder": pl.Utf8,
    })


def _run_harmonisation(args, out_path: Path, manifest_file_path: str, log_file: Path) -> None:
    """
    Run postgwas harmonisation docker once using the generated manifest.
    """
    print("\n[STEP 6] Running harmonisation (Docker)...", flush=True)

    resource_folder = Path(args.resource_folder).resolve()
    defaults_path = Path(args.defaults_config).resolve()

    cmd = [
    "docker", "run", "--rm", "--platform", "linux/amd64",
    "-e", "PYTHONUNBUFFERED=1",   # 🔥 REQUIRED
    "-v", f"{out_path}:{out_path}",
    "-v", f"{resource_folder}:{resource_folder}",
    "-v", f"{defaults_path.parent}:{defaults_path.parent}",
    "jibinjv/postgwas:1.3",
    "postgwas", "harmonisation",
    "--nthreads", str(args.nthreads),
    "--max-mem", str(args.max_mem),
    "--config", manifest_file_path,
    "--defaults", str(defaults_path),
    ]

    result = _run_subprocess_with_log(
        cmd=cmd,
        log_file=log_file,
        step_name="POSTGWAS HARMONISATION",
    )

    if result.returncode != 0:
        raise RuntimeError("Harmonisation failed")


def mtag_pipeline_runner(args):
    """
    Validate manifest → merge VCFs → create MTAG inputs → run MTAG → rename outputs → harmonise
    """
    out_path = Path(args.out)
    out_path.mkdir(parents=True, exist_ok=True)

    pre_process_dir = out_path / f"0_{args.run_name}_pre_process"
    pre_process_dir.mkdir(parents=True, exist_ok=True)

    mtag_input_dir = out_path / f"1_{args.run_name}_mtag_input"
    mtag_input_dir.mkdir(parents=True, exist_ok=True)

    mtag_output_dir = out_path / f"2_{args.run_name}_mtag_results"
    mtag_output_dir.mkdir(parents=True, exist_ok=True)

    manifest_dir = out_path / "harmonisation_input_files"
    manifest_dir.mkdir(parents=True, exist_ok=True)

    out_tsv = pre_process_dir / f"{args.run_name}_master_sumstats.tsv"
    preprocess_log_file = pre_process_dir / f"{args.run_name}_master_sumstats.log"
    output_prefix = str(mtag_output_dir / f"{args.run_name}_mtag_results")
    mtag_log_file = mtag_output_dir / f"{args.run_name}_mtag_results.log"
    manifest_file_path = str(manifest_dir / f"{args.run_name}_manifest_out.tsv")

    manifest_df = _build_empty_manifest_df()
    sample_ids = []
    sumstats_list = []

    # ------------------------------------------------------------------
    # STEP 1: Validate manifest
    # ------------------------------------------------------------------
    try:
        print("\n[STEP 1] Validating manifest...", flush=True)
        val_input = validate_manifest(args.inputfile, str(out_path))
    except Exception as e:
        raise RuntimeError("Manifest validation failed.") from e

    # ------------------------------------------------------------------
    # STEP 2: Merge VCFs → master TSV
    # ------------------------------------------------------------------
    try:
        print("\n[STEP 2] Merging individual VCFs into a master TSV...", flush=True)

        master_vcf_tsv, success_2 = run_merge_vcf_to_tsv_pipeline(
            input_df=val_input["manifest_df"],
            out_tsv=str(out_tsv),
            fixed_fields=["ID", "CHROM", "POS", "REF", "ALT"],
            format_fields=["AF", "ES", "SE", "NEF", "LP", "SI", "EZ"],
            log_file=str(preprocess_log_file),
            n_threads=int(args.nthreads),
            qc_n_rows=5_000_000,
            bcftools_exec="bcftools",
            fail_fast_on_merge=True,
        )

        if not success_2:
            raise RuntimeError("Master TSV creation failed. Cannot proceed to MTAG.")
    except Exception as e:
        raise RuntimeError("STEP 2 failed: merge VCF to master TSV failed.") from e

    # ------------------------------------------------------------------
    # STEP 3: Create MTAG inputs  
    # ------------------------------------------------------------------
    ## effective sample size is recomended https://github.com/JonJala/mtag/issues/239  ; https://github.com/JonJala/mtag/issues/60 
    try:
        print("\n[STEP 3] Creating MTAG specific file format from master TSV...", flush=True)

        mtag_inputs = run_mastertsv_to_toolformat_pipeline(
            mode=("joined_samples", "sample_specific"),
            tsv_file=master_vcf_tsv,
            input_manifest_df=val_input["manifest_df"],
            joined_out_tsv=str(mtag_input_dir / f"{args.run_name}_mtag_ready.tsv"),
            sample_specific_directory=str(mtag_input_dir),
            joined_remove_duplicate_ids=True,
            joined_create_unique_id=False,
            joined_sample_retained_cols=["ID", "AF", "NEF", "SI"],
            joined_header_map=None,
            joined_suffix_replace=None,
            joined_convert_to_raw_p_value=True,
            sample_remove_duplicate_ids=True,
            sample_create_unique_id=False,
            sample_specific_retained_cols=[
                "ID", "CHROM", "POS", "REF", "ALT", "AF", "NEF", "LP", "SI", "EZ"
            ],
            sample_header_map=None,
            sample_suffix_replace=None,
            sample_convert_to_raw_p_value=True,
            log_file=str(preprocess_log_file),
        )

        sample_ids = mtag_inputs["sample_specific"]["sample_id"].to_list()
        sumstats_list = [
            str(p) for p in mtag_inputs["sample_specific"]["matrix_tsv"].to_list()
        ]
        sumstats_inputs = ",".join(sumstats_list)

    except Exception as e:
        raise RuntimeError("STEP 3 failed: MTAG input creation failed.") from e

    # ------------------------------------------------------------------
    # STEP 4: Run MTAG
    # ------------------------------------------------------------------
    try:
        print("\n[STEP 4] Running MTAG statistical engine...", flush=True)

        MTAG_ROOT = "/opt/mtag"

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
            "--cores", str(int(args.nthreads)),
            "--chunksize", str(int(args.chunksize)),
            "--stream_stdout",
            "--force",
        ]

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

        if getattr(args, "incld_ambig_snps", False):
            cmd_mtag.append("--incld_ambig_snps")

        print(f"DEBUG COMMAND: {' '.join(cmd_mtag)}\n", flush=True)

        _run_subprocess_with_log(
            cmd=cmd_mtag,
            log_file=mtag_log_file,
            step_name="MTAG STATISTICAL ENGINE",
        )

        print("\n[SUCCESS] MTAG statistical engine finished. Proceeding to renaming...", flush=True)

    except Exception as e:
        raise RuntimeError("MTAG statistical analysis failed.") from e

    # ------------------------------------------------------------------
    # STEP 5: Rename outputs + build manifest
    # ------------------------------------------------------------------
    try:
        print("[STEP 5] Preparing input files for harmonisation...", flush=True)

        _append_log(mtag_log_file, "\n" + "=" * 130)
        _append_log(mtag_log_file, f"MTAG OUTPUT RENAMING - {datetime.datetime.now()}")
        _append_log(mtag_log_file, "=" * 130)

        if not (getattr(args, "perfect_gencov", False) and getattr(args, "equal_h2", False)):
            for i, sample_id in enumerate(sample_ids):
                trait_idx = i + 1
                generic_path = Path(f"{output_prefix}_trait_{trait_idx}.txt")
                final_filename = f"{args.run_name}_mtag_results_trait_{trait_idx}_{sample_id}.txt"
                final_full_path = mtag_output_dir / final_filename

                if generic_path.exists():
                    os.rename(generic_path, str(final_full_path))
                    _append_log(mtag_log_file, f"trait_{trait_idx} → {sample_id} → {final_full_path}")
                    print(f"   - Mapped trait_{trait_idx} to {sample_id}", flush=True)
                else:
                    _append_log(mtag_log_file, f"trait_{trait_idx} → {sample_id} → NOT FOUND ({generic_path})")
                    print(f"   - [WARNING] trait_{trait_idx} not found for {sample_id}", flush=True)
                    continue

                _, manifest_df = create_manifest_with_info(
                    sumstat_file=str(final_full_path),
                    gwas_outputname=f"{args.run_name}_MTAG_meta_{sample_id}",

                    chr_col="CHR",
                    pos_col="BP",
                    snp_id_col="SNP",
                    ea_col="A1",
                    oa_col="A2",
                    eaf_col="FRQ",

                    beta_or_col="mtag_beta",
                    se_col="mtag_se",
                    imp_z_col="mtag_z",
                    pval_col="mtag_pval",

                    ncontrol_col="N",
                    ncase_col="NA",
                    ncontrol="NA",
                    ncase="NA",

                    eaffile="NA",
                    eafcolumn="NA",
                    chr_pos_col="NA",
                    liftover="Yes",

                    sumstat_merge_keys=["SNP"],
                    info_merge_keys=["ID"],

                    resource_folder=args.resource_folder,
                    output_folder=args.out,

                    info_files=[str(mtag_input_dir / f"{sample_id}_matrix.tsv")],
                    info_suffix="SI",
                    info_method=args.info_method,
                    final_info_col="info_score", 

                    manifest_df=manifest_df,
                    sep="\t",
                )

        else:
            print(f"[DEBUG] perfect_gencov={args.perfect_gencov}", flush=True)
            print(f"[DEBUG] equal_h2={args.equal_h2}", flush=True)
            final_full_path = mtag_output_dir / f"{args.run_name}_mtag_results_mtag_meta.txt" 
            print(final_full_path)
            print(str(mtag_input_dir / f"{args.run_name}_mtag_ready.tsv"))

            if not final_full_path.exists():
                raise FileNotFoundError(f"Expected MTAG meta output not found: {final_full_path}")

            _, manifest_df = create_manifest_with_info(
                sumstat_file=str(final_full_path),
                gwas_outputname=f"{args.run_name}_MTAG_meta",

                chr_col="CHR",
                pos_col="BP",
                snp_id_col="SNP",
                ea_col="A1",
                oa_col="A2",
                eaf_col="meta_freq",

                beta_or_col="mtag_beta",
                se_col="mtag_se",
                imp_z_col="mtag_z",
                pval_col="mtag_pval",

                ncontrol_col="sample_size",
                ncase_col="NA",
                ncontrol="NA",
                ncase="NA",

                eaffile="NA",
                eafcolumn="NA",
                chr_pos_col="NA",
                liftover="Yes",

                sumstat_merge_keys=["SNP"],
                info_merge_keys=["ID"],
                sample_size_merge_keys=["ID"],

                resource_folder=args.resource_folder,
                output_folder=args.out,

                info_files=[str(mtag_input_dir / f"{args.run_name}_mtag_ready.tsv")],
                info_suffix="SI",
                info_method=args.info_method,
                final_info_col="info_score",

                sample_size_files=[str(mtag_input_dir / f"{args.run_name}_mtag_ready.tsv")],
                sample_size_suffix="NEF",
                sample_size_method=args.n_method,
                final_sample_size_col="sample_size",

                manifest_df=manifest_df,
                sep="\t",
            )

        manifest_df.write_csv(manifest_file_path, separator="\t")
        _append_log(mtag_log_file, f"Manifest written to: {manifest_file_path}")

    except Exception as e:
        raise RuntimeError("STEP 5 failed: output renaming / manifest generation failed.") from e

    # ------------------------------------------------------------------
    # STEP 6: Harmonisation
    # ------------------------------------------------------------------
    try:
        _run_harmonisation(
            args=args,
            out_path=out_path,
            manifest_file_path=manifest_file_path,
            log_file=mtag_log_file,
        )
    except Exception as e:
        raise RuntimeError("Harmonisation failed.") from e

    print(f"\n[FINISH] MTAG workflow completed. Audit trail: {mtag_log_file}", flush=True)

    return {
        "sample_ids": sample_ids,
        "mtag_input_files": sumstats_list,
        "mtag_results_dir": str(mtag_output_dir),
        "manifest_file": manifest_file_path,
        "log_file": str(mtag_log_file),
    }

# import six.moves.BaseHTTPServer 
# from posixpath import sep 
# import os
# import subprocess
# import datetime
# import logging
# from pathlib import Path

# from postgwas_pleio.utils.input_manifest_val import validate_manifest
# from postgwas_pleio.formatter.merge_vcf_to_tsv import run_merge_vcf_to_tsv_pipeline
# from postgwas_pleio.formatter.master_tsv_to_tool_format import run_mastertsv_to_toolformat_pipeline
# from postgwas_pleio.utils.harmonise_input import create_manifest_with_info

# logger = logging.getLogger(__name__)


# def mtag_pipeline_runner(args):
#     """
#     Validate manifest → merge VCFs → create MTAG inputs → run MTAG → rename outputs to sample IDs.
#     """

#     # ---------------------------
#     # STEP 1: Validate manifest
#     # ---------------------------
#     out_path = Path(args.out)
#     out_path.mkdir(parents=True, exist_ok=True)

#     print("\n[STEP 1] Validating manifest...", flush=True)
#     val_input = validate_manifest(args.inputfile, str(out_path))

#     # ---------------------------
#     # STEP 2: Merge VCFs → master TSV
#     # ---------------------------
#     print("\n[STEP 2] Merging individual VCFs into a master TSV...", flush=True)

#     pre_process_dir = out_path / f"0_{args.run_name}_pre_process"
#     pre_process_dir.mkdir(parents=True, exist_ok=True)

#     out_tsv = pre_process_dir / f"{args.run_name}_master_sumstats.tsv"
#     preprocess_log_file = pre_process_dir / f"{args.run_name}_master_sumstats.log"

#     master_vcf_tsv, success_2 = run_merge_vcf_to_tsv_pipeline(
#         input_df=val_input["manifest_df"],
#         out_tsv=str(out_tsv),
#         fixed_fields=["ID", "CHROM", "POS", "REF", "ALT"],
#         format_fields=["AF", "ES", "SE", "NEF", "LP", "SI", "EZ"],
#         log_file=str(preprocess_log_file),
#         n_threads=int(args.nthreads),
#         qc_n_rows=5_000_000,
#         bcftools_exec="bcftools",
#         fail_fast_on_merge=True,
#     )

#     if not success_2:
#         raise RuntimeError("Master TSV creation failed. Cannot proceed to MTAG.")

#     # ---------------------------
#     # STEP 3: Create MTAG inputs
#     # ---------------------------
#     print("\n[STEP 3] Creating MTAG specific file format from master TSV...", flush=True)

#     mtag_input_dir = out_path / f"1_{args.run_name}_mtag_input"
#     mtag_input_dir.mkdir(parents=True, exist_ok=True)

#     mtag_inputs = run_mastertsv_to_toolformat_pipeline(
#         mode={"joined_samples", "sample_specific"},
#         tsv_file=master_vcf_tsv,
#         input_manifest_df=val_input["manifest_df"],
#         joined_out_tsv=str(mtag_input_dir / f"{args.run_name}_mtag_ready.tsv"),
#         sample_specific_directory=str(mtag_input_dir),
#         joined_remove_duplicate_ids=True,
#         joined_create_unique_id=False,
#         joined_sample_retained_cols=["ID", "CHROM", "POS", "REF", "ALT", "AF", "NC","SS", "LP", "SI", "EZ"],
#         joined_header_map=None,
#         joined_suffix_replace=None,
#         joined_convert_to_raw_p_value=True,
#         sample_remove_duplicate_ids=True,
#         sample_create_unique_id=False,
#         sample_specific_retained_cols=["ID", "CHROM", "POS", "REF", "ALT", "AF", "NC","SS", "LP", "SI", "EZ"],
#         sample_header_map=None,
#         sample_suffix_replace=None,
#         sample_convert_to_raw_p_value=True,
#         log_file=str(preprocess_log_file),
#     )

#     # Use the *actual* order of sample-specific outputs (guarantees trait_1 ↔ sample_id match)
#     sample_ids = mtag_inputs["sample_specific"]["sample_id"].to_list()

#     # Ensure sumstats list are strings
#     sumstats_list = [str(p) for p in mtag_inputs["sample_specific"]["matrix_tsv"].to_list()]
#     sumstats_inputs = ",".join(sumstats_list)

#     # ---------------------------
#     # STEP 4: Run MTAG
#     # ---------------------------
#     print("\n[STEP 4] Running MTAG statistical engine...", flush=True)

#     MTAG_ROOT = "/opt/mtag"
#     mtag_output_dir = out_path / f"2_{args.run_name}_mtag_results"
#     mtag_output_dir.mkdir(parents=True, exist_ok=True)

#     output_prefix = str(mtag_output_dir / f"{args.run_name}_mtag_results")
#     mtag_log_file = str(mtag_output_dir / f"{args.run_name}_mtag_results.log")

#     cmd_mtag = [
#         "micromamba", "run", "-n", "mtag",
#         "python", f"{MTAG_ROOT}/mtag.py",
#         "--sumstats", sumstats_inputs,
#         "--out", output_prefix,
#         "--snp_name", "ID",
#         "--z_name", "EZ",
#         "--beta_name", "ES",
#         "--se_name", "SE",
#         "--n_name", "NEF",
#         "--eaf_name", "AF",
#         "--chr_name", "CHROM",
#         "--bpos_name", "POS",
#         "--a1_name", "ALT",
#         "--a2_name", "REF",
#         "--p_name", "LP",
#         "--cores", str(int(args.nthreads)),
#         "--chunksize", str(int(args.chunksize)),
#         "--stream_stdout",
#         "--force",
#     ]

#     # Optional arguments
#     if getattr(args, "ld_ref_panel", None):
#         cmd_mtag.extend(["--ld_ref_panel", str(args.ld_ref_panel)])
#     if getattr(args, "no_overlap", False):
#         cmd_mtag.append("--no_overlap")
#     if getattr(args, "perfect_gencov", False):
#         cmd_mtag.append("--perfect_gencov")
#     if getattr(args, "equal_h2", False):
#         cmd_mtag.append("--equal_h2")
#     if getattr(args, "std_betas", False):
#         cmd_mtag.append("--std_betas")
#     if getattr(args, "fdr", False):
#         cmd_mtag.append("--fdr")
#     if getattr(args, "incld_ambig_snps", False):
#         cmd_mtag.append("--incld_ambig_snps")

#     print(f"DEBUG COMMAND: {' '.join(cmd_mtag)}\n", flush=True)


#     manifest_df = pl.DataFrame(
#         schema={  "sumstat_file": pl.Utf8, "gwas_outputname": pl.Utf8, "chr_col": pl.Utf8, "pos_col": pl.Utf8, "snp_id_col": pl.Utf8, 
#         "ea_col": pl.Utf8, "oa_col": pl.Utf8, "eaf_col": pl.Utf8, "beta_or_col": pl.Utf8, "se_col": pl.Utf8, "imp_z_col": pl.Utf8, 
#         "pval_col": pl.Utf8, "ncontrol_col": pl.Utf8, "ncase_col": pl.Utf8, "ncontrol": pl.Utf8, "ncase": pl.Utf8, "imp_info_col": pl.Utf8,
#         "infofile": pl.Utf8, "infocolumn": pl.Utf8, "eaffile": pl.Utf8, "eafcolumn": pl.Utf8, "liftover": pl.Utf8, "chr_pos_col": pl.Utf8,
#         "resourse_folder": pl.Utf8, "output_folder": pl.Utf8,
#         }
#     )

#     try:
#         subprocess.run(cmd_mtag, check=True)
#         print("\n[SUCCESS] MTAG statistical engine finished. Proceeding to renaming...", flush=True)

#         # ---------------------------
#         # STEP 5: Rename MTAG outputs
#         # ---------------------------
#         print("[STEP 5] Mapping generic results to sample names...", flush=True)
#         mtag_outputs={}
#         with open(mtag_log_file, "a") as f_log:
#             f_log.write("\n" + "=" * 130 + "\n")
#             f_log.write(f"MTAG OUTPUT RENAMING - {datetime.datetime.now()}\n")
#             f_log.write("=" * 130 + "\n")
            
#             if not (getattr(args, "perfect_gencov", False) and getattr(args, "equal_h2", False)):
#                 for i, sample_id in enumerate(sample_ids):
#                     trait_idx = i + 1
#                     generic_path = f"{output_prefix}_trait_{trait_idx}.txt"

#                     final_filename = f"{args.run_name}_mtag_results_trait_{trait_idx}_{sample_id}.txt"
#                     final_full_path = mtag_output_dir / final_filename
#                     mtag_outputs[sample_id]=final_full_path 

#                     # MarkerName	Allele1	Allele2	Freq1	FreqSE	MinFreq	MaxFreq	Weight	Zscore	N	P-value	Direction	HetISq	HetChiSq	HetDf	HetPVal	TotalNEF
#                     try:
#                         # Create manifest row
#                         _, manifest_df = create_manifest_with_info(
#                             sumstat_file=final_full_path,
#                             gwas_outputname=f"{args.run_name}_MTAG_meta",

#                             chr_col="CHR",
#                             pos_col="BP",
#                             snp_id_col="SNP",
#                             ea_col="A1",
#                             oa_col="A2",
#                             eaf_col="FRQ",

#                             beta_or_col="mtag_beta", 
#                             se_col="mtag_se",
#                             imp_z_col="mtag_z", 
#                             pval_col="mtag_pval",

#                             ncontrol_col="N", 
#                             ncase_col="NA",
#                             ncontrol="NA",
#                             ncase="NA",

#                             eaffile="NA",
#                             eafcolumn="NA",
#                             chr_pos_col="NA",
#                             liftover="Yes",

#                             sumstat_merge_keys=["MarkerName"],
#                             info_merge_keys=["ID"],

#                             resource_folder=args.resource_folder,
#                             output_folder=args.out,

#                             info_files=[str(mtag_input_dir / f"{sample_id}_matrix.tsv")],
#                             info_suffix="SI",
#                             info_method=args.info_method,
#                             final_info_col="info_score",

#                             manifest_df=manifest_df,
#                             sep="\t",
#                         )
#                     except Exception as e:
#                         print(f"[ERROR] Harmonisation input manifest file creation failed for {sample_id}: {e}")
#                         sys.exit(1)

#                     manifest_file_path = str(
#                         out_path / f"harmonisation_input_files/{args.run_name}_manifest_out.tsv" )
#                     manifest_df.write_csv(manifest_file_path, separator="\t")


#                     # ------------------------------------------------
#                     # STEP 7: Run harmonisation (The Docker-in-Docker Part)
#                     # ------------------------------------------------
#                     print("\n[STEP 7] Running harmonisation (Docker-in-Docker)...", flush=True)

#                     # Fix 2: Resolve all paths relative to the HOST mapping
#                     resource_folder = Path(args.resource_folder).resolve()
#                     defaults_path = Path(args.defaults_config).resolve()
#                     # metal_out_dir and out_path are already resolved/created above

#                     if not os.path.exists("/var/run/docker.sock"):
#                         print("CRITICAL WARNING: /var/run/docker.sock NOT FOUND.")

#                     cmd = [
#                         "docker", "run", "--rm", "--platform", "linux/amd64",
#                         "-v", f"{out_path}:{out_path}",
#                         "-v", f"{resource_folder}:{resource_folder}",
#                         "-v", f"{defaults_path.parent}:{defaults_path.parent}",
#                         "jibinjv/postgwas:1.3",
#                         "postgwas", "harmonisation",
#                         "--nthreads", str(args.cores),
#                         "--max-mem", str(args.max_mem),
#                         "--config", str(manifest_file_path),
#                         "--defaults", str(defaults_path)
#                     ]

#                     try:
#                         subprocess.run(cmd, check=True)

#                     except subprocess.CalledProcessError:
#                         print("Error: Harmonisation container failed.")
#                         raise

#                     # Rename generic file   
#                     if os.path.exists(generic_path):
#                         os.rename(generic_path, str(final_full_path))
#                         f_log.write(f"trait_{trait_idx} → {sample_id} → {final_full_path}\n")
#                         print(f"   - Mapped trait_{trait_idx} to {sample_id}", flush=True)
#                     else:
#                         f_log.write(f"trait_{trait_idx} → {sample_id} → NOT FOUND ({generic_path})\n")
#                         print(f"   - [WARNING] trait_{trait_idx} not found for {sample_id}", flush=True)
#             else:
#                 pass
                
#             f_log.write("=" * 130 + "\n")

#         print(f"\n[FINISH] MTAG workflow completed. Audit trail: {mtag_log_file}", flush=True)

#     except subprocess.CalledProcessError as e:
#         print(f"\n[ERROR] MTAG failed with exit code {e.returncode}", flush=True)
#         raise RuntimeError("MTAG statistical analysis failed.") from e

#     return {
#         "sample_ids": sample_ids,
#         "mtag_input_files": sumstats_list,
#         "mtag_results_dir": str(mtag_output_dir),
#         "log_file": mtag_log_file,
#     }