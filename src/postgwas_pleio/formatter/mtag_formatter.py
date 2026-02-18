
import subprocess
import textwrap
import os
import re
import logging
import datetime
from pathlib import Path
from typing import Dict, List
import polars as pl

logger = logging.getLogger(__name__)

import polars as pl
import os
import logging
from typing import List

logger = logging.getLogger(__name__)


def split_wide_gwas_to_samples(input_path: str, output_dir: str) -> List[str]:
    print(f"--- Processing Master Wide File: {input_path} ---", flush=True)
    # Load with comment_prefix=None so #[1]ID isn't ignored
    # --- ADD THIS CHECK ---
    if not os.path.exists(input_path):
        print(f"!!! ERROR: File not found inside container at {input_path} !!!", flush=True)
        return []
    file_size = os.path.getsize(input_path) / (1024 * 1024)
    print(f"--- File found! Size: {file_size:.2f} MB. Loading into Polars... ---", flush=True)
    # -----------------------
    df = pl.read_csv(
        input_path, 
        separator="\t", 
        null_values=["", "."], 
        truncate_ragged_lines=True,
        infer_schema_length=10000,
        comment_prefix=None  
    )
    # Filter for columns where the null count is NOT equal to the total number of rows
    df = df.select([
        pl.col(name) for name in df.columns 
        if df[name].null_count() < df.height
    ])  
    # STEP 1: Split by ']' and retain only the second half
    # This transforms "[2]CHROM" -> "CHROM" and "[6]i42:AF" -> "i42:AF"
    cleaned_headers = []
    for col in df.columns:
        if "]" in col:
            cleaned_headers.append(col.split("]")[-1].strip())
        else:
            cleaned_headers.append(col.strip())
    df.columns = cleaned_headers
    # STEP 2: Rename the very first column to "ID"
    # This handles whatever is left of the "#[1]ID" string
    if df.columns[0] != "ID":
        cols = df.columns
        cols[0] = "ID"
        df.columns = cols
    print(f"Standardized Headers: {df.columns[:10]}...", flush=True)
    # Define fixed anchors
    fixed_cols = ["ID", "CHROM", "POS", "REF", "ALT"]
    # STEP 3: Identify sample names (i42, i65, etc.) from columns with ':'
    sample_cols = [col for col in df.columns if ":" in col]
    unique_samples = sorted(list(set(col.split(":")[0] for col in sample_cols)))
    print(f"Detected Samples for splitting: {unique_samples}", flush=True)
    os.makedirs(output_dir, exist_ok=True)
    generated_files = []
    # STEP 4: Process each sample
    for sample in unique_samples:
        # Select fixed columns AND columns where the name contains the sample name
        this_sample_cols = [c for c in sample_cols if f"{sample}:" in c]
        sample_df = df.select([
            pl.col(fixed_cols),
            pl.col(this_sample_cols)
        ])
        # STEP 5: Split by ':' to finalize internal names
        # "i42:AF" -> "AF", "i42:ES" -> "ES"
        new_cols = sample_df.columns
        new_cols = [
            col.split(":")[-1] if col in this_sample_cols else col
            for col in new_cols
        ]
        sample_df.columns = new_cols
        # Standard conversion for MTAG
        if "LP" in sample_df.columns:
            sample_df = sample_df.with_columns([
                (10.0 ** (-pl.col("LP"))).alias("P")
            ])
        # Final write
        output_file = os.path.join(output_dir, f"{sample}_mtag_ready.tsv")
        sample_df.write_csv(output_file, separator="\t")
        generated_files.append(output_file)
        print(f"Written: {sample} ({len(sample_df)} SNPs)", flush=True)
    return generated_files

def mtag_formatter(sumstat_vcfs: List[str], output_folder: str, run_name: str,n_threads=3) -> Dict[str, str]:
    """
    Automated GWAS VCF-to-TSV formatter for MTAG.
    Logs environmental metadata and the specific shell command used.
    """
    output_dir = Path(output_folder)
    output_dir.mkdir(parents=True, exist_ok=True)
    master_tsv = output_dir / f"{run_name}_master_wide.tsv"
    log_file = output_dir / f"{run_name}_pipeline.log"
    split_dir = output_dir / "mtag_inputs"
    # Safely quote file paths for the shell
    vcf_input_str = " ".join([f"'{str(v)}'" for v in sumstat_vcfs])
    # Construct the command
    cmd = textwrap.dedent(f"""
        set -e
        bcftools merge --threads {n_threads} -m none {vcf_input_str} | \\
        bcftools view --threads {n_threads} -e 'FMT/ES == "."' | \\
        bcftools view --threads {n_threads} --min-alleles 2 --max-alleles 2 | \\
        bcftools query --print-header \\
        -f '%ID\\t%CHROM\\t%POS\\t%REF\\t%ALT\\t[%AF\\t]\\t[%ES\\t]\\t[%SE\\t]\\t[%NEF\\t]\\t[%LP\\t]\\t[%EZ\\t]\\n' | \\
        sed 's/\\t$//' > "{master_tsv}" 
    """).strip()
    try:
        with open(log_file, "w") as f_log:
            # --- START LOGGING METADATA ---
            f_log.write(f"PostGWAS-Pleio Pipeline Log\n")
            f_log.write(f"Run Name: {run_name}\n")
            f_log.write(f"Timestamp: {datetime.datetime.now().isoformat()}\n")
            f_log.write("-" * 50 + "\n")
            # Log bcftools version for debugging
            try:
                ver = subprocess.check_output(["bcftools", "--version"], text=True).splitlines()[0]
                f_log.write(f"Software Version: {ver}\n")
            except:
                f_log.write("Software Version: bcftools version unknown\n")
            # LOG THE COMMAND
            f_log.write(f"\nEXECUTING COMMAND:\n{cmd}\n")
            f_log.write("-" * 50 + "\n\n")
            f_log.flush()
            # --- EXECUTE ---
            logger.info(f"Running bcftools merge. Logging to: {log_file}")
            subprocess.run(
                cmd, 
                shell=True, 
                check=True, 
                executable="/bin/bash", 
                stdout=f_log, 
                stderr=f_log
            )
            f_log.write("\n[SUCCESS] bcftools merge and query completed.\n")
    except subprocess.CalledProcessError as e:
        logger.error(f"VCF Formatting failed. Check log: {log_file}")
        with open(log_file, "a") as f_log:
            f_log.write(f"\n[ERROR] Command failed with exit code {e.returncode}\n")
        raise RuntimeError(f"Pipeline failed at bcftools step. See {log_file}")
    # Step 2: Split and Convert -log10(p) -> p
    sample_files = split_wide_gwas_to_samples(str(master_tsv), str(split_dir))
    return {
        "master_tsv": str(master_tsv),
        "individual_files": sample_files,
        "mtag_path_list": ",".join(sample_files),
        "log_file": str(log_file)
    }