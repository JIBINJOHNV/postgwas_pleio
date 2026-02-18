import subprocess
import logging
import sys # Added for flushing
from postgwas_pleio.formatter.mtag_formatter import mtag_formatter

logger = logging.getLogger(__name__)


import os
import subprocess
import datetime
from pathlib import Path

def mtag_pipeline_runner(args):
    """
    Executes MTAG and renames the generic trait_N outputs to original sample IDs
    while logging full path metadata for a complete audit trail.
    """
    # 1. Run the formatter
    print(f"\n[STEP 1] Formatting VCFs into MTAG format...", flush=True)
    formatter_results = mtag_formatter(
        sumstat_vcfs=args.vcfs, 
        output_folder=args.out, 
        run_name=args.run_name
    )

    # Extract original sample names from formatter paths
    # Example: '/path/to/i42_mtag_ready.tsv' -> 'i42'
    original_paths = formatter_results["individual_files"]
    sample_ids = [os.path.basename(f).replace("_mtag_ready.tsv", "") for f in original_paths]
    
    sumstats_list = formatter_results["mtag_path_list"]
    MTAG_ROOT = "/opt/mtag"
    output_prefix = f"{args.out}/{args.run_name}_mtag_results"
    log_file = formatter_results["log_file"]
    
    # 2. Build the Base Command
    cmd_mtag = [
        "micromamba", "run", "-n", "mtag",
        "python", f"{MTAG_ROOT}/mtag.py",
        "--sumstats", sumstats_list,
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
        "--p_name", "P",
        "--cores", str(args.cores),
        "--chunksize", str(int(args.chunksize)),
        "--stream_stdout",
        "--force",
    ]

    # Add optional arguments
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

    # 3. Execution with Real-Time Feedback
    print(f"[STEP 2] Starting MTAG statistical engine...", flush=True)
    print(f"DEBUG COMMAND: {' '.join(cmd_mtag)}\n", flush=True)
    
    try:
        subprocess.run(cmd_mtag, check=True)
        print("\n[SUCCESS] MTAG statistical engine finished. Proceeding to renaming...", flush=True)
        
        # 4. Final Renaming and Full-Path Logging
        print(f"[STEP 3] Mapping generic results to sample names and logging paths...", flush=True)
        
        with open(log_file, "a") as f_log:
            f_log.write("\n" + "="*130 + "\n")
            f_log.write(f"MTAG OUTPUT RENAMING & PATH MAPPING - {datetime.datetime.now()}\n")
            f_log.write("="*130 + "\n")
            
            # Table Header for the log
            header = f"{'GENERIC ID':<12} | {'ORIGINAL SAMPLE':<20} | {'MTAG OUTPUT FILENAME':<40} | {'FINAL FILENAME'}\n"
            f_log.write(header)
            f_log.write("-" * 130 + "\n")

            for i, sample_id in enumerate(sample_ids):
                trait_idx = i + 1
                generic_path = f"{output_prefix}_trait_{trait_idx}.txt"
                
                # Construct final filename and absolute path
                final_filename = f"{args.run_name}_mtag_results_trait_{trait_idx}_{sample_id}.txt"
                final_full_path = os.path.abspath(os.path.join(args.out, final_filename))
                
                if os.path.exists(generic_path):
                    os.rename(generic_path, final_full_path)
                    
                    # Log the full audit trail
                    log_entry = (
                        f"trait_{trait_idx:<7} | "
                        f"{sample_id:<20} | "
                        f"{generic_path:<40} | "
                        f"{final_full_path}\n"
                    )
                    f_log.write(log_entry)
                    print(f"   - Mapped trait_{trait_idx} to {sample_id}", flush=True)
                else:
                    error_msg = f"trait_{trait_idx:<7} | {sample_id:<20} | ERROR: {generic_path} NOT FOUND\n"
                    f_log.write(error_msg)
                    print(f"   - [!] Warning: trait_{trait_idx} not found for {sample_id}", flush=True)

            f_log.write("="*130 + "\n")
        
        print(f"\n[FINISH] MTAG workflow completed. Audit trail saved to: {log_file}", flush=True)

    except subprocess.CalledProcessError as e:
        print(f"\n[ERROR] MTAG failed with exit code {e.returncode}", flush=True)
        raise RuntimeError("MTAG statistical analysis failed.")

    return {
        "sample_ids": sample_ids,
        "formatted_files": formatter_results["individual_files"],
        "mtag_results_dir": args.out,
        "log_file": log_file
    }

    
# def mtag_pipeline_runner(args):
#     # 1. Run the formatter
#     print(f"\n[STEP 1] Formatting VCFs into MTAG format...", flush=True)
#     formatter_results = mtag_formatter(
#         sumstat_vcfs=args.vcfs, 
#         output_folder=args.out, 
#         run_name=args.run_name
#     )

#     sumstats_list = formatter_results["mtag_path_list"]
#     MTAG_ROOT = "/opt/mtag"
    
#     # 2. Build the Base Command
#     cmd_mtag = [
#         "micromamba", "run", "-n", "mtag",
#         "python", f"{MTAG_ROOT}/mtag.py",
#         "--sumstats", sumstats_list,
#         "--out", f"{args.out}/{args.run_name}_mtag_results",
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
#         "--p_name", "P",
#         "--cores", str(args.cores),
#         "--chunksize", str(int(args.chunksize)),
#         "--stream_stdout",
#         "--force",
#     ]

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
#     # 3. Execution with Real-Time Feedback
#     print(f"[STEP 2] Starting MTAG statistical engine...", flush=True)
#     print(f"DEBUG COMMAND: {' '.join(cmd_mtag)}\n", flush=True)
    
#     try:
#         # We use a direct call. MTAG's --stream_stdout will now show up because of ENTRYPOINT logic
#         subprocess.run(cmd_mtag, check=True)
#         print("\n[SUCCESS] MTAG workflow completed successfully.", flush=True)
#     except subprocess.CalledProcessError as e:
#         print(f"\n[ERROR] MTAG failed with exit code {e.returncode}", flush=True)
#         raise RuntimeError("MTAG statistical analysis failed.")

#     return {
#         "formatted_files": formatter_results["individual_files"],
#         "mtag_results_prefix": f"{args.out}/{args.run_name}_mtag_results"
#     }