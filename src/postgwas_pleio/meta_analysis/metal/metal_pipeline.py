import os
import subprocess
import logging
from pathlib import Path
import shlex
import polars as pl

# Move these to the top to avoid the "UnboundLocalError" you saw earlier
from postgwas_pleio.utils.input_manifest_val import validate_manifest
from postgwas_pleio.formatter.merge_vcf_to_tsv import run_merge_vcf_to_tsv_pipeline
from postgwas_pleio.formatter.master_tsv_to_tool_format import run_mastertsv_to_toolformat_pipeline
from postgwas_pleio.meta_analysis.metal.script_generator import generate_metal_script
#from postgwas_pleio.meta_analysis.metal.harmonise_input import metal_sumstat_to_harmonise_input
from postgwas_pleio.utils.harmonise_input import create_manifest_with_info

logger = logging.getLogger(__name__)

def _run_cmd_with_tee_log(cmd, log_path: Path, *, cwd: Path | None = None) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log:
        log.write("COMMAND:\n" + " ".join(map(str, cmd)) + "\n\n")
        log.flush()
        proc = subprocess.Popen(
            list(map(str, cmd)),
            cwd=str(cwd) if cwd else None,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            print(line, end="", flush=True)
            log.write(line)
        rc = proc.wait()
        if rc != 0:
            raise subprocess.CalledProcessError(rc, cmd)

def metal_pipeline_runner(args):
    # Fix 1: Always resolve the main output path immediately
    out_path = Path(args.out).resolve()
    out_path.mkdir(parents=True, exist_ok=True)

    # --- STEPS 1-6 (Logic remains the same, but use absolute paths) ---
    print("\n[STEP 1] Validating manifest...", flush=True)
    val = validate_manifest(args.inputfile, str(out_path))
    manifest_df = val["manifest_df"]

    print("\n[STEP 2] Merging individual VCFs...", flush=True)
    pre_process_dir = out_path / f"0_{args.run_name}_pre_process"
    pre_process_dir.mkdir(parents=True, exist_ok=True)
    master_tsv = pre_process_dir / f"{args.run_name}_master_sumstats.tsv"
    preprocess_log_file = pre_process_dir / f"{args.run_name}_master_sumstats.log"

    master_vcf_tsv, ok = run_merge_vcf_to_tsv_pipeline(
        input_df=manifest_df,
        out_tsv=str(master_tsv),
        fixed_fields=["ID", "CHROM", "POS", "REF", "ALT"],
        format_fields=["AF", "ES", "SE", "NEF", "LP", "SI", "EZ"],
        log_file=str(preprocess_log_file),
        n_threads=int(args.nthreads),
        qc_n_rows=5_000_000,
        bcftools_exec="bcftools",
        fail_fast_on_merge=True,
        require_variant_in_all_samples=False
    )
    if not ok: raise RuntimeError("Master TSV creation failed.")

    print("\n[STEP 3] Creating METAL inputs...", flush=True)
    metal_input_dir = out_path / f"1_{args.run_name}_metal_input"
    metal_input_dir.mkdir(parents=True, exist_ok=True)
    metal_inputs = run_mastertsv_to_toolformat_pipeline(
        mode=["sample_specific", "joined_samples"],
        tsv_file=master_vcf_tsv,
        input_manifest_df=manifest_df,

        joined_out_tsv=str(metal_input_dir / f"{args.run_name}_joined_samples.tsv"),

        sample_specific_directory=str(metal_input_dir),

        joined_remove_duplicate_ids=True, 
        joined_create_unique_id=True,
        joined_sample_retained_cols=["ID","CHROM","POS", "SI"],
        joined_header_map=None,
        joined_suffix_replace=None, 
        joined_convert_to_raw_p_value=True,

        sample_remove_duplicate_ids=True,
        sample_create_unique_id=True,
        sample_specific_retained_cols=[
            "ID", "REF", "ALT", "ES", "SE", "NEF", "AF", "LP", "EZ"
        ],
        sample_header_map=None,
        sample_suffix_replace=None,
        sample_convert_to_raw_p_value=True,

        log_file=str(preprocess_log_file),
    )

    
    sample_ids = metal_inputs["sample_specific"]["sample_id"].to_list()
    sumstats_list = [str(p) for p in metal_inputs["sample_specific"]["matrix_tsv"].to_list()]

    # Save sample order for reproducibility
    sample_order_path = metal_input_dir / f"{args.run_name}_sample_order.tsv"

    with sample_order_path.open("w", encoding="utf-8") as f:
        f.write("trait_index\tsample_id\tsumstats_tsv\n")
        for i, (sid, p) in enumerate(zip(sample_ids, sumstats_list), start=1):
            f.write(f"{i}\t{sid}\t{p}\n")


    print("\n[STEP 4] Generating METAL script...", flush=True)
    metal_out_dir = out_path / f"2_{args.run_name}_metal_output"
    metal_out_dir.mkdir(parents=True, exist_ok=True)
    output_prefix = metal_out_dir / f"{args.run_name}_metal_out"
    script_path = metal_out_dir / f"{args.run_name}_metal_script.txt"
    metal_log = metal_out_dir / f"{args.run_name}_metal_run.log"

    generate_metal_script(
        input_files=sumstats_list,
        output_prefix=str(output_prefix),
        script_path=str(script_path),
        scheme=args.scheme,
        genomic_control=args.genomic_control,
        overlap_correction=args.overlap_correction,
        track_freq=args.track_freq,
        heterogeneity=args.heterogeneity,
    )

    print("\n[STEP 5] Running METAL...", flush=True)
    _run_cmd_with_tee_log(["metal", str(script_path)], metal_log, cwd=metal_out_dir)

    metal_result_file = metal_out_dir / f"{args.run_name}_metal_out_1.tbl"
    if not metal_result_file.exists(): raise RuntimeError("METAL output not generated")


    # ------------------------------------------------
    # STEP 6: Prepare harmonisation input
    # ------------------------------------------------
    print("\n[STEP 6] Preparing harmonisation input...", flush=True)

    manifest_df = pl.DataFrame(
        schema={  "sumstat_file": pl.Utf8, "gwas_outputname": pl.Utf8, "chr_col": pl.Utf8, "pos_col": pl.Utf8, "snp_id_col": pl.Utf8, 
        "ea_col": pl.Utf8, "oa_col": pl.Utf8, "eaf_col": pl.Utf8, "beta_or_col": pl.Utf8, "se_col": pl.Utf8, "imp_z_col": pl.Utf8, 
        "pval_col": pl.Utf8, "ncontrol_col": pl.Utf8, "ncase_col": pl.Utf8, "ncontrol": pl.Utf8, "ncase": pl.Utf8, "imp_info_col": pl.Utf8,
        "infofile": pl.Utf8, "infocolumn": pl.Utf8, "eaffile": pl.Utf8, "eafcolumn": pl.Utf8, "liftover": pl.Utf8, "chr_pos_col": pl.Utf8,
        "resourse_folder": pl.Utf8, "output_folder": pl.Utf8,
        }
    )
    
    try:    
        if args.scheme.upper() == "SAMPLESIZE":
            if args.sample_size_approach.lower() == "weight":
                ncontrol_col = "Weight"

            elif args.sample_size_approach.lower() == "totalnef":
                ncontrol_col = "TotalNEF"

            elif args.sample_size_approach.lower() == "sampleoverlapcorrected":
                ncontrol_col = "N"

            else:
                raise ValueError(
                    f"Unsupported sample_size_approach: {args.sample_size_approach}. "
                    "Expected one of: weight, totalnef, sampleoverlapcorrected"
                )

            _, manifest_df = create_manifest_with_info(
                sumstat_file=metal_result_file,
                gwas_outputname=f"{args.run_name}_METAL_meta",

                chr_col="CHROM",
                pos_col="POS",
                snp_id_col="MarkerName",
                ea_col="Allele1",
                oa_col="Allele2",
                eaf_col="Freq1",

                beta_or_col="NA", 
                se_col="NA",
                imp_z_col="Zscore", 
                pval_col="P-value",

                ncontrol_col=ncontrol_col, 
                ncase_col="NA",
                ncontrol="NA",
                ncase="NA",

                eaffile="NA",
                eafcolumn="NA",
                chr_pos_col="NA",
                liftover="Yes",

                sumstat_merge_keys=["MarkerName"],
                info_merge_keys=["ID"],

                resource_folder=args.resource_folder,
                output_folder=args.out,

                info_files=[str(metal_input_dir / f"{args.run_name}_joined_samples.tsv")],
                info_suffix="SI",
                extra_info_file_columns=["CHROM","POS"],
                info_method=args.info_method,
                final_info_col="info_score",

                manifest_df=manifest_df,
                sep="\t",
            )

        elif args.scheme.upper() == "STDERR":
            _, manifest_df = create_manifest_with_info(
                sumstat_file=metal_result_file,
                gwas_outputname=f"{args.run_name}_METAL_meta",

                chr_col="CHROM",
                pos_col="POS",
                snp_id_col="MarkerName",
                ea_col="Allele1",
                oa_col="Allele2",
                eaf_col="Freq1",

                beta_or_col="Effect",
                se_col="StdErr",
                imp_z_col="NA",
                pval_col="P-value",

                ncontrol_col='TotalNEF',
                ncase_col="NA",
                ncontrol="NA",
                ncase="NA",

                eaffile="NA",
                eafcolumn="NA",
                chr_pos_col="NA",
                liftover="Yes",

                sumstat_merge_keys=["MarkerName"],
                info_merge_keys=["ID"],

                resource_folder=args.resource_folder,
                output_folder=args.out,

                info_files=[str(metal_input_dir / f"{args.run_name}_joined_samples.tsv")],
                info_suffix="SI",
                extra_info_file_columns=["CHROM","POS"],
                info_method=args.info_method,
                final_info_col="info_score",

                manifest_df=manifest_df,
                sep="\t",
            )

        else:
            raise ValueError(f"Unsupported METAL scheme: {args.scheme}")

        manifest_file_path = str(
            out_path / f"harmonisation_input_files/{args.run_name}_manifest_out.tsv"
        )

        manifest_df.write_csv(manifest_file_path, separator="\t")

    except Exception as e:
        raise RuntimeError(
            f"Failed during manifest creation for METAL scheme '{args.scheme}': {e}"
        )

    # ------------------------------------------------
    # STEP 7: Run harmonisation (The Docker-in-Docker Part)
    # ------------------------------------------------
    print("\n[STEP 7] Running harmonisation (Docker-in-Docker)...", flush=True)

    # Fix 2: Resolve all paths relative to the HOST mapping
    resource_folder = Path(args.resource_folder).resolve()
    defaults_path = Path(args.defaults_config).resolve()
    # metal_out_dir and out_path are already resolved/created above

    if not os.path.exists("/var/run/docker.sock"):
        print("CRITICAL WARNING: /var/run/docker.sock NOT FOUND.")

    cmd = [
        "docker", "run", "--rm", "--platform", "linux/amd64",
        "-v", f"{out_path}:{out_path}",
        "-v", f"{resource_folder}:{resource_folder}",
        "-v", f"{defaults_path.parent}:{defaults_path.parent}",
        "jibinjv/postgwas:1.3",
        "postgwas", "harmonisation",
        "--nthreads", str(args.nthreads),
        "--max-mem", str(args.max_mem),
        "--config", str(manifest_file_path),
        "--defaults", str(defaults_path)
    ]

    # -------------------------------------------------
    # Print command exactly as it would appear in shell
    # -------------------------------------------------
    print("\nRunning command:")
    print(" ".join(shlex.quote(c) for c in cmd))
    print()

    try:
        subprocess.run(cmd, check=True)

    except subprocess.CalledProcessError:
        print("Error: Harmonisation container failed.")
        raise

    except FileNotFoundError:
        print("Error: 'docker' command missing. Did you rebuild the image with docker-cli?")
        raise

    print("\n✅ METAL pipeline completed successfully.", flush=True)