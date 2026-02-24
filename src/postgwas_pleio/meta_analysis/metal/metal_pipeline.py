
import subprocess
import logging
from pathlib import Path

from postgwas_pleio.utils.input_manifest_val import validate_manifest
from postgwas_pleio.formatter.merge_vcf_to_tsv import run_merge_vcf_to_tsv_pipeline
from postgwas_pleio.formatter.master_tsv_to_tool_format import run_mastertsv_to_toolformat_pipeline
from postgwas_pleio.meta_analysis.metal.script_generator import generate_metal_script

logger = logging.getLogger(__name__)


def _run_cmd_with_tee_log(cmd, log_path: Path, *, cwd: Path | None = None) -> None:
    """
    Run a command and tee stdout/stderr to a log file (and console).
    Minimal + robust. Raises on non-zero exit.
    """
    log_path.parent.mkdir(parents=True, exist_ok=True)

    with log_path.open("w", encoding="utf-8") as log:
        log.write("COMMAND:\n")
        log.write(" ".join(map(str, cmd)) + "\n\n")
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
    """
    Validate manifest → merge VCFs → create METAL per-study inputs → generate METAL script → run METAL.

    Notes:
      - METAL produces ONE meta-analysis output (not one output per study).
      - We therefore write sample_order.tsv to preserve input order / provenance.
    """
    out_path = Path(args.out)
    out_path.mkdir(parents=True, exist_ok=True)

    # ---------------------------
    # STEP 1: Validate manifest
    # ---------------------------
    print("\n[STEP 1] Validating manifest...", flush=True)
    val = validate_manifest(args.inputfile, str(out_path))
    manifest_df = val["manifest_df"]

    # ---------------------------
    # STEP 2: Merge VCFs → master TSV
    # ---------------------------
    print("\n[STEP 2] Merging individual VCFs into a master TSV...", flush=True)

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
        n_threads=int(args.cores),
        qc_n_rows=5_000_000,
        bcftools_exec="bcftools",
        fail_fast_on_merge=True,
    )
    if not ok:
        raise RuntimeError("Master TSV creation failed. Cannot proceed to METAL.")

    # ---------------------------
    # STEP 3: Create METAL per-study inputs (sample_specific)
    # ---------------------------
    print("\n[STEP 3] Creating METAL per-study inputs from master TSV...", flush=True)

    metal_input_dir = out_path / f"1_{args.run_name}_metal_input"
    metal_input_dir.mkdir(parents=True, exist_ok=True)

    # Deterministic: use LIST, not set
    metal_inputs = run_mastertsv_to_toolformat_pipeline(
        mode=["sample_specific"],  # ✅ METAL needs one file per study; joined_samples not required
        tsv_file=master_vcf_tsv,
        input_manifest_df=manifest_df,
        # joined_* not used because mode excludes joined_samples
        joined_out_tsv=str(metal_input_dir / f"{args.run_name}_joined_unused.tsv"),
        sample_specific_directory=str(metal_input_dir),
        joined_remove_duplicate_ids=True,
        joined_create_unique_id=False,
        joined_sample_retained_cols=["ID","CHROM", "POS", "REF", "ALT", "ES", "SE", "NEF", "AF", "LP", "EZ", "SI"],
        joined_header_map=None,
        joined_suffix_replace=None,
        joined_convert_to_raw_p_value=True,
        sample_remove_duplicate_ids=True,
        sample_create_unique_id=True,
        # Keep only what METAL needs (you can keep extras, but minimal is best)
        sample_specific_retained_cols=["ID", "REF", "ALT", "ES", "SE", "NEF", "AF", "LP", "EZ"], 
        sample_header_map=None,
        sample_suffix_replace=None,
        sample_convert_to_raw_p_value=True,
        log_file=str(preprocess_log_file),
    )

    # Use the actual order produced by formatter
    sample_ids = metal_inputs["sample_specific"]["sample_id"].to_list()
    sumstats_list = [str(p) for p in metal_inputs["sample_specific"]["matrix_tsv"].to_list()]

    # Write provenance / order sidecar (important for reproducibility)
    sample_order_path = metal_input_dir / f"{args.run_name}_sample_order.tsv"
    with sample_order_path.open("w", encoding="utf-8") as f:
        f.write("trait_index\tsample_id\tsumstats_tsv\n")
        for i, (sid, p) in enumerate(zip(sample_ids, sumstats_list), start=1):
            f.write(f"{i}\t{sid}\t{p}\n")

    # ---------------------------
    # STEP 4: Generate METAL script
    # ---------------------------
    print("\n[STEP 4] Generating METAL script...", flush=True)

    metal_out_dir = out_path / f"2_{args.run_name}_metal_output"
    metal_out_dir.mkdir(parents=True, exist_ok=True)

    output_prefix = metal_out_dir / f"{args.run_name}_metal_out"
    script_path = metal_out_dir / f"{args.run_name}_metal_script.txt"
    metal_log = metal_out_dir / f"{args.run_name}_metal_run.log"

    script_path = generate_metal_script(
        input_files=sumstats_list,
        output_prefix=str(output_prefix),
        script_path=str(script_path),
        scheme=args.scheme,
        genomic_control=args.genomic_control,
        overlap_correction=args.overlap_correction,
        zcutoff=args.zcutoff,
        verbose=args.verbose,
        column_counting=args.column_counting,
        track_freq=args.track_freq,
        heterogeneity=args.heterogeneity,
    )

    # ---------------------------
    # STEP 5: Run METAL
    # ---------------------------
    print("\n[STEP 5] Running METAL...", flush=True)

    # METAL expects: metal <scriptfile>
    # Use cwd=metal_out_dir so METAL writes outputs there reliably.
    _run_cmd_with_tee_log(["metal", str(script_path)], metal_log, cwd=metal_out_dir)

    print("\n✅ METAL pipeline completed.", flush=True)
    print(f" - METAL script : {script_path}")
    print(f" - METAL log    : {metal_log}")
    print(f" - Sample order : {sample_order_path}")
    print(f" - Output prefix: {output_prefix}", flush=True)