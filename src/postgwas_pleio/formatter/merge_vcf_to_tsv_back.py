from postgwas_pleio.utils.logging_utils import setup_logger, log_action, log_detected, log_warning, log_error
import subprocess,sys
import textwrap
import polars as pl
from pathlib import Path
import pandas as pd



def normalize_input_df(input_df):

    # pandas → polars
    try:
        import pandas as pd
        if isinstance(input_df, pd.DataFrame):
            return pl.from_pandas(input_df)
    except ImportError:
        pass
    # polars → keep
    if isinstance(input_df, pl.DataFrame):
        return input_df
    raise TypeError("input_df must be pandas or polars DataFrame")

    
# =========================================================
# 1️⃣ Variant count per VCF
# =========================================================
def count_vcf_variants(vcf, logger):

    log_action(logger, f"Counting variants → {vcf}")

    try:
        out = subprocess.check_output(
            ["bcftools", "index", "-n", vcf],
            text=True
        ).strip()

        n = int(out)
        log_detected(logger, f"Variant count (index) → {n:,}")
        return n

    except Exception:

        log_warning(logger, f"Index missing → fallback scan → {vcf}")

        cmd = f"bcftools view -H {vcf} | wc -l"
        out = subprocess.check_output(cmd, shell=True, text=True)

        n = int(out.strip())
        log_detected(logger, f"Variant count (scan) → {n:,}")
        return n


# =========================================================
# 2️⃣ Build bcftools query format
# =========================================================
def build_query_format(fixed_fields, format_fields):

    fixed_map = {
        "ID": "%ID",
        "CHROM": "%CHROM",
        "POS": "%POS",
        "REF": "%REF",
        "ALT": "%ALT"
    }

    fixed_part = "\t".join(fixed_map[f] for f in fixed_fields)
    format_part = "".join([f"[%{f}\\t]" for f in format_fields])

    return f"{fixed_part}\\t{format_part}\\n"


# =========================================================
# 3️⃣ Merge + query
# =========================================================
def merge_vcf_to_tsv(vcf_list, out_tsv, fixed_fields, format_fields, logger, n_threads=4):

    log_action(logger, "Starting VCF merge")

    query_fmt = build_query_format(fixed_fields, format_fields)
    vcf_str = " ".join(vcf_list)

    cmd = textwrap.dedent(f"""
        set -e
        bcftools merge --threads {n_threads} -m none {vcf_str} | \\
        bcftools view --threads {n_threads} -e 'FMT/ES == "."' | \\
        bcftools view --threads {n_threads} --min-alleles 2 --max-alleles 2 | \\
        bcftools query --print-header -f '{query_fmt}' | \\
        sed 's/\\t$//' > "{out_tsv}"
    """)

    log_action(logger, "Executing bcftools merge")

    try:
        subprocess.run(cmd, shell=True, executable="/bin/bash", check=True)
        log_detected(logger, "VCF merge completed successfully")

    except subprocess.CalledProcessError:
        log_error(logger, "VCF merge failed")
        raise


# =========================================================
# 4️⃣ Count TSV variants
# =========================================================
def count_tsv_rows(tsv, logger):

    log_action(logger, f"Counting TSV rows → {tsv}")

    out = subprocess.check_output(f"wc -l {tsv}", shell=True, text=True)
    n = int(out.split()[0]) - 1

    log_detected(logger, f"Total variants in TSV → {n:,}")
    return n


# =========================================================
# 5️⃣ Partial QC scan
# =========================================================


def partial_qc_scan(tsv, logger, n_rows=1_000_000):

    log_action(logger, f"Partial QC scan → first {n_rows:,} rows")

    try:
        df = pl.read_csv(
            tsv,
            separator="\t",
            null_values=["", "."],
            truncate_ragged_lines=True,
            infer_schema_length=10000
        ).head(n_rows)

    except Exception as e:

        log_error(logger, f"TSV read failed → {tsv}")
        log_error(logger, str(e))

        # return safe values → caller will mark QC failure
        return [], ["QC_READ_FAILED"]

    # header cleanup
    try:
        df.columns = [x.split("]")[-1] for x in df.columns]
    except Exception:
        log_warning(logger, "Header cleanup failed")

    # empty column detection
    try:
        empty_cols = [c for c in df.columns if df[c].null_count() == df.height]
    except Exception as e:
        log_error(logger, "Empty column detection failed")
        log_error(logger, str(e))
        return df.columns, ["QC_EMPTY_CHECK_FAILED"]

    log_detected(logger, f"Empty columns detected → {len(empty_cols)}")

    return df.columns, empty_cols


# =========================================================
# 6️⃣ Sample ID extraction
# =========================================================
def extract_samples(columns, fixed_fields):

    samples = set()

    for c in columns:
        if c not in fixed_fields and ":" in c:
            samples.add(c.split(":")[0])
    return samples


# =========================================================
# 7️⃣ Full QC wrapper (PRODUCTION)
# =========================================================
def run_merge_vcf_to_tsv_pipeline(
    input_df,
    out_tsv,
    fixed_fields=None,
    format_fields=None,
    log_file=None
):

    # ⭐ normalize input type
    input_df = normalize_input_df(input_df)
    
    if fixed_fields is None:
        fixed_fields = ["ID", "CHROM", "POS", "REF", "ALT"]

    if format_fields is None:
        format_fields = ["AF", "ES", "SE", "NEF", "LP", "SI", "EZ"]

    if log_file is None:
        log_file = Path(out_tsv).with_suffix(".log")

    logger = setup_logger(log_file)

    qc_failures = []
    cols = []

    # ---------- STEP 1 ----------
    log_action(logger, "STEP 1 → Variant count per sample")

    for row in input_df.iter_rows(named=True):
        try:
            n = count_vcf_variants(row["sumstat_vcf"], logger)
            log_detected(logger, f"{row['sample_id']} → {n:,}")
        except Exception as e:
            qc_failures.append(f"Variant count failed → {row['sample_id']}")
            log_error(logger, str(e))

    # ---------- STEP 2 ----------
    try:
        merge_vcf_to_tsv(
            input_df["sumstat_vcf"].to_list(),
            out_tsv,
            fixed_fields,
            format_fields,
            logger
        )
    except Exception as e:
        qc_failures.append("Merge failed")
        log_error(logger, str(e))

    # ---------- STEP 3 ----------
    try:
        count_tsv_rows(out_tsv, logger)
    except Exception:
        qc_failures.append("TSV count failed")

    # ---------- STEP 4 ----------
    try:
        cols, empty_cols = partial_qc_scan(out_tsv, logger)

        for c in empty_cols:
            qc_failures.append(f"Empty column → {c}")

    except Exception:
        qc_failures.append("Partial QC failed")

    # ---------- STEP 5 ----------
    try:
        tsv_samples = extract_samples(cols, fixed_fields)
        input_samples = set(input_df["sample_id"])

        missing = input_samples - tsv_samples
        extra = tsv_samples - input_samples

        if missing:
            qc_failures.append(f"Missing samples → {missing}")

        if extra:
            qc_failures.append(f"Unexpected samples → {extra}")

    except Exception:
        qc_failures.append("Sample extraction failed")

    # ---------- FINAL ----------
    log_action(logger, "FINAL QC SUMMARY")

    if qc_failures:
        for f in qc_failures:
            log_error(logger, f)

        logger.error("Merging vcf and converting to TSV FAILED ")
        logger.info("=" * 80)
        sys.exit(1)

    else:
        log_action(logger, "PIPELINE COMPLETED SUCCESSFULLY")
        return Path(out_tsv), True


# input_df = pd.read_csv("/Users/JJOHN41/Downloads/kadoorie_biobank/vcf_files/pleio_input.txt", sep="\t")
# run_merge_vcf_to_tsv_pipeline(
#     input_df=input_df,
#     out_tsv="/Users/JJOHN41/Downloads/kadoorie_biobank/vcf_files/kadoorie_biobank_test.tsv",
#     fixed_fields=None,
#     format_fields=None,
#     log_file=None
# )