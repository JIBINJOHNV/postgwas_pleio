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
        -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\t[%AF\t][%NEF\t][%LP\t][%SI\t][%EZ\t]\n' | \\
        sed 's/\\t$//' > "{master_tsv}" 
    """).strip()
    try:
        with open(log_file, "w") as f_log:
            f_log.write(f"FastASSET Pipeline Initialization: {run_name}\n")
            f_log.write("-" * 50 + "\n")
            f_log.flush()
            subprocess.run(cmd, shell=True, check=True, executable="/bin/bash", stdout=f_log, stderr=f_log)
            # Fixed: Passed unique_samples argument
            master_tsv_head = pd.read_csv(master_tsv, sep="\t", nrows=10)
            vcf_ids=[x for x in master_tsv_head.columns if x.endswith('_LP')]
            vcf_ids=[x.replace('_LP', '') for x in vcf_ids]

            # Use sets to determine mismatches 
            input_ids = set(input_sample_ids)
            missing_from_vcf = list(input_ids - vcf_ids)
            extra_in_vcf = list(vcf_ids - input_ids)
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
            # SPREV/PPREV Handling
            # 1️⃣ Non-binary → NA
            sample_files_df.loc[sample_files_df['TYPE'] != 'binary', ['SPREV','PPREV']] = np.nan
            # 2️⃣ Binary + missing SPREV → default 0.5
            mask_sprev = sample_files_df['TYPE'].eq('binary') & sample_files_df['SPREV'].isna()
            sample_files_df.loc[mask_sprev, 'SPREV'] = 0.5
            # 3️⃣ Binary + missing PPREV → collect IDs
            mask_missing_pprev = sample_files_df['TYPE'].eq('binary') & sample_files_df['PPREV'].isna()
            missing_pop_prevalence = sample_files_df.loc[mask_missing_pprev, 'sample_id'].tolist()
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

