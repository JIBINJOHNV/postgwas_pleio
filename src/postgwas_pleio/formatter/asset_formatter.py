import os
import subprocess
import textwrap
import datetime
import logging
from pathlib import Path
from typing import Dict, List, Tuple
import polars as pl
import polars.selectors as cs
import numpy as np
import pandas as pd  # Fixed: Added missing import

logger = logging.getLogger(__name__)

def split_wide_gwas_to_samples(input_path: str, output_dir: str, run_name: str, log_file: str, unique_samples: List[str]) -> Tuple[str, Dict[str, str]]:
    """
    Cleans the master wide file, tracks variant counts, and logs progress.
    Returns: (fastasset_master_file_path, dict_of_individual_trait_paths)
    """
    def log_message(message: str):
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        formatted_msg = f"[{timestamp}] {message}"
        print(formatted_msg, flush=True)
        with open(log_file, "a") as f:
            f.write(formatted_msg + "\n")
    log_message(f"--- STEP: Loading Master Wide File: {input_path} ---")
    if not os.path.exists(input_path):
        log_message(f"CRITICAL ERROR: File {input_path} not found.")
        raise FileNotFoundError(f"Master TSV not found at {input_path}")
    df = pl.read_csv(
        input_path, 
        separator="\t", 
        null_values=["", "."], 
        truncate_ragged_lines=True,
        infer_schema_length=10000,
        comment_prefix=None  
    )
    initial_count = df.height
    log_message(f"Initial variants loaded from VCF merge: {initial_count:,}")
    # Filter for columns where the null count is NOT equal to the total number of rows
    df = df.select([
        pl.col(name) for name in df.columns 
        if df[name].null_count() < df.height
    ])  
    # STEP 1: Split by ']' and retain only the second half
    df.columns=[x.split("]")[1] for x in df.columns]
    log_message("Standardized headers and resolved bcftools column naming.")
    # Note: Based on previous printout, A2 is REF and A1 is ALT
    df = df.with_columns(
        pl.when(pl.col("ID").is_null())
        .then(
            pl.format("{}_{}_{}_{}", 
                pl.col("CHROM"), 
                pl.col("POS"), 
                pl.col("REF"), 
                pl.col("ALT")
            )
        )
        .otherwise(pl.col("ID"))
        .alias("ID")
    )
    log_message("Generated ID column from CHROM_POS_REF_ALT if missing ID COLUMN")
    # STEP 2: FastASSET-specific cleaning
    log_message("Filtering for biallelic SNPs and deduplicating...")
    df = df.with_columns(avg_AF = pl.mean_horizontal(cs.ends_with(":AF")))
    df_final = (
        df.filter(
            (pl.col("ID").count().over("ID") == 1) | 
            ((pl.col("REF").str.len_chars() == 1) & (pl.col("ALT").str.len_chars() == 1))
        )
        .sort("avg_AF", descending=True)
        .unique(subset=["ID"], keep="first")
    )
    final_count = df_final.height
    removed_count = initial_count - final_count
    log_message(f"Variants removed: {removed_count:,} | Retained: {final_count:,}")
    # STEP 3: Create Master FastASSET Input (Wide format)
    os.makedirs(output_dir, exist_ok=True)
    df_fastasset = (
        df.drop("ID")
          .with_columns(
              pl.format("{}_{}_{}_{}", 
                  pl.col("CHROM"), 
                  pl.col("POS"), 
                  pl.col("REF"), 
                  pl.col("ALT")
              ).alias("ID")
          )
        ) 
    df_fastasset = df_fastasset.select([
        pl.col("ID"),
        cs.ends_with(":ES").name.map(lambda x: x.replace(":ES", ".Beta")),
        cs.ends_with(":SE").name.map(lambda x: x.replace(":SE", ".SE")),
        cs.ends_with(":NEF").name.map(lambda x: x.replace(":NEF", ".N"))
    ])
    fastasset_input_file = os.path.abspath(os.path.join(output_dir, f"1_fastasset_inputs/{run_name}_fastasset_input.tsv"))
    df_fastasset.write_csv(fastasset_input_file, separator="\t")
    log_message(f"SUCCESS: Generated FastASSET Master Input: {fastasset_input_file}")
    # STEP 4: Individual file splitting
    log_message("Splitting into individual trait-level TSVs...")
    fixed_cols = ["ID", "CHROM", "POS", "REF", "ALT"]
    header_map = { "REF": "A2", "ALT": "A1", "ID": "SNP", "ES": "effect", "EZ": "Z", "AF": "MAF", "NEF": "N","SI":"INFO","raw_P":"P"}
    sample_cols = [col for col in df_final.columns if ":" in col]
    unique_samples_detected = sorted(list(set(col.split(":")[0] for col in sample_cols)))
    sample_files = {}
    for sample in unique_samples_detected:
        this_sample_cols = [c for c in sample_cols if c.startswith(f"{sample}:")]
        sample_df = df_final.select(fixed_cols + this_sample_cols)
        sample_df.columns = [c.split(":")[-1] if ":" in c else c for c in sample_df.columns]
        if "LP" in sample_df.columns:
            sample_df = sample_df.with_columns([(10.0 ** (-pl.col("LP"))).alias("raw_P")]).drop("LP")
        
        existing_map = {old: new for old, new in header_map.items() if old in sample_df.columns}
        #sample_df = sample_df.rename(existing_map)
        columns_expr = []
        for col in sample_df.columns:
            if col in existing_map:
                columns_expr.append(pl.col(col).alias(existing_map[col]))
            else:
                columns_expr.append(pl.col(col))
        sample_df = sample_df.select(columns_expr)
        output_file = os.path.abspath(os.path.join(output_dir, f"2_fastasset_ldsc_inputs/{sample}_fastasset_ldsc.tsv"))
        sample_df.write_csv(output_file, separator="\t")
        sample_files[sample] = output_file
        log_message(f"   - Saved trait file for: {sample} ({sample_df.height:,} variants)")
    log_message("--- STEP: split_wide_gwas_to_samples completed successfully ---") 
    return fastasset_input_file, sample_files



def asset_formatter(inputfile: str, output_folder: str, run_name: str, n_threads=3) -> Dict[str, str]:
    output_dir = Path(output_folder)
    preprocess_dir = output_dir / "0_preprocess"
    fastasset_input_dir = output_dir / "1_fastasset_inputs"
    fastasset_ldsc_input_dir = output_dir / "2_fastasset_ldsc_inputs"
    preprocess_dir.mkdir(parents=True, exist_ok=True) 
    fastasset_input_dir.mkdir(parents=True, exist_ok=True) 
    fastasset_ldsc_input_dir.mkdir(parents=True, exist_ok=True)
    
    master_tsv = preprocess_dir / f"{run_name}_master_wide.tsv"
    log_file = output_dir / f"{run_name}_pipeline.log"
    input_df = pd.read_csv(inputfile, sep=None, engine='python')
    vcf_files = input_df['sumstat_vcf'].tolist() 
    input_sample_ids = input_df['sample_id'].tolist()
    
    for vcf_file in vcf_files:
        if not os.path.exists(vcf_file):
            raise FileNotFoundError(f"VCF file not found: {vcf_file}")
    
    vcf_input_str = " ".join([f"'{str(v)}'" for v in vcf_files])
    cmd = textwrap.dedent(f"""
        set -e
        bcftools merge --threads {n_threads} -m none {vcf_input_str} | \\
        bcftools view --threads {n_threads} -e 'FMT/ES == "."' | \\
        bcftools view --threads {n_threads} --min-alleles 2 --max-alleles 2 | \\
        bcftools query --print-header \\
        -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\t[%AF\t][%ES\t][%SE\t][%NEF\t][%LP\t][%SI\t][%EZ\t]\n' | \\
        sed 's/\\t$//' > "{master_tsv}" 
    """).strip()
    try:
        with open(log_file, "w") as f_log:
            f_log.write(f"FastASSET Pipeline Initialization: {run_name}\n")
            f_log.write("-" * 50 + "\n")
            f_log.flush()
            subprocess.run(cmd, shell=True, check=True, executable="/bin/bash", stdout=f_log, stderr=f_log)
            # Fixed: Passed unique_samples argument
            fastasset_master_file, sample_files = split_wide_gwas_to_samples(
                str(master_tsv), 
                str(output_folder), 
                run_name, 
                str(log_file),
                unique_samples=input_sample_ids
            )
            # --- MODIFIED SECTION START ---
            sample_files_df = pd.DataFrame(sample_files.items(), columns=["NAME", "FILE"])
            
            # 1. Perform merge
            sample_files_df2 = pd.merge(input_df, sample_files_df, left_on="sample_id", right_on="NAME")
            
            # Use sets to determine mismatches
            input_ids = set(input_sample_ids)
            vcf_ids = set(sample_files_df['NAME'])
            matched_ids = set(sample_files_df2['sample_id'])
            missing_from_vcf = list(input_ids - matched_ids)
            extra_in_vcf = list(vcf_ids - matched_ids)
            # --- CONSOLIDATED MISMATCH PRINTING ---
            if missing_from_vcf or extra_in_vcf:
                print("\n" + "!" * 60)
                print("WARNING: SAMPLE MISMATCH DETECTED")
                if missing_from_vcf:
                    print(f"\n[!] MISSING FROM VCF: {len(missing_from_vcf)} samples requested but NOT found:")
                    for s in missing_from_vcf:
                        print(f"    - {s}")
                if extra_in_vcf:
                    print(f"\n[i] EXTRA IN VCF: {len(extra_in_vcf)} samples found but NOT requested:")
                    for s in extra_in_vcf:
                        print(f"    - {s}")
                # Now decide whether to exit
                if missing_from_vcf:
                    print("\n" + "!" * 60)
                    print("CRITICAL ERROR: REQUIRED SAMPLES MISSING. EXITING PROGRAM.")
                    print("Please align your input CSV with the VCF sample headers.")
                    print("!" * 60 + "\n")
                    sys.exit(1)
                else:
                    print("\nProceeding with matched samples only...")
                    print("!" * 60 + "\n")
            # --- COLUMN SELECTION AND PREVALENCE LOGIC ---
            cols_to_keep = [c for c in sample_files_df2.columns if c in input_df.columns or c == 'FILE']
            sample_files_df2 = sample_files_df2[cols_to_keep]
            # SPREV/PPREV Handling
            sample_files_df2['SPREV'] = np.where(sample_files_df2['TYPE'] == 'binary', 0.5, sample_files_df2.get('SPREV', np.nan))
            # FILL MISSING WITH "nan" STRING
            sample_files_df2['PPREV'] = sample_files_df2['PPREV'].fillna("nan")
            sample_files_df2['SPREV'] = sample_files_df2['SPREV'].fillna("nan")
            mask_missing_pprev = (sample_files_df2['TYPE'] == 'binary') & (sample_files_df2['PPREV'] == "nan")
            missing_pop_prevalence = sample_files_df2[mask_missing_pprev]['sample_id'].tolist()
            if missing_pop_prevalence:
                print("\n" + "#" * 60)
                print("DATA QUALITY WARNING: MISSING PREVALENCE")
                print(f"The following binary traits have MISSING PPREV: {missing_pop_prevalence}")
                print("Note: Both Population and Sample Prevalence will be passed as 'nan'.")
                print("This is NOT IDEAL for liability scale transformation.")
                print("CONTINUING PIPELINE REGARDLESS...")
                print("#" * 60 + "\n")
            # --- MODIFIED SECTION END ---           
            print(f"Missing Pop Prevalence in binary traits: {missing_pop_prevalence}")
            manifest_path = os.path.join(fastasset_input_dir, f"{run_name}_fastasset_ldsc_input.txt")
            sample_files_df2.rename(columns={"sample_id":"NAME"}, inplace=True)
            sample_files_df2.to_csv(manifest_path, index=False,sep="\t")
            return {
                "master_tsv": str(master_tsv),
                "missing_pop_prevalence": missing_pop_prevalence,
                "fastasset_input": fastasset_master_file,
                "log_file": str(log_file),
                "fastasset_sample_manifest": manifest_path
            }
    except Exception as e:
        logger.error(f"FastASSET formatting failed: {e}")
        raise

