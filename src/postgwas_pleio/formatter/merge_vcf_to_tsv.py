from postgwas_pleio.utils.logging_utils import (
    setup_logger,
    log_action,
    log_detected,
    log_warning,
    log_error,
)
import subprocess
import sys
import textwrap
import shlex
import polars as pl
from pathlib import Path

# Optional dependency: pandas (only needed if caller passes pandas DataFrame)
try:
    import pandas as pd  # type: ignore
except Exception:  # pragma: no cover
    pd = None


# =========================================================
# 0️⃣ Normalize input DataFrame type
# =========================================================
def normalize_input_df(input_df):
    """
    Accept pandas or polars DataFrame; return polars DataFrame.
    """
    if pd is not None and isinstance(input_df, pd.DataFrame):
        return pl.from_pandas(input_df)

    if isinstance(input_df, pl.DataFrame):
        return input_df

    raise TypeError("input_df must be a pandas.DataFrame or polars.DataFrame")


# =========================================================
# 1️⃣ Variant count per VCF
# =========================================================
def count_vcf_variants(vcf, logger, bcftools_exec="bcftools"):
    """
    Count variants in a bgzipped/indexed VCF using:
      1) bcftools index -n
      2) fallback: bcftools view -H | wc -l
    """
    vcf = str(vcf)
    log_action(logger, f"Counting variants → {vcf}")

    try:
        out = subprocess.check_output(
            [bcftools_exec, "index", "-n", vcf],
            text=True,
        ).strip()

        n = int(out)
        log_detected(logger, f"Variant count (index) → {n:,}")
        return n

    except Exception:
        log_warning(logger, f"Index missing → fallback scan → {vcf}")

        # shell-safe quoting
        vcf_q = shlex.quote(vcf)
        bcftools_q = shlex.quote(str(bcftools_exec))
        cmd = f"{bcftools_q} view -H {vcf_q} | wc -l"

        out = subprocess.check_output(cmd, shell=True, text=True)
        n = int(out.strip())
        log_detected(logger, f"Variant count (scan) → {n:,}")
        return n


# =========================================================
# 2️⃣ Build bcftools query format
# =========================================================
def build_query_format(fixed_fields, format_fields):
    """
    Build format string for bcftools query -f ...
    """
    fixed_map = {
        "ID": "%ID",
        "CHROM": "%CHROM",
        "POS": "%POS",
        "REF": "%REF",
        "ALT": "%ALT",
    }

    fixed_part = "\t".join(fixed_map[f] for f in fixed_fields)
    format_part = "".join([f"[%{f}\\t]" for f in format_fields])

    return f"{fixed_part}\\t{format_part}\\n"


# =========================================================
# 3️⃣ Merge + query
# =========================================================
def merge_vcf_to_tsv(
    vcf_list,
    out_tsv,
    fixed_fields,
    format_fields,
    logger,
    n_threads=4,
    bcftools_exec="bcftools",
):
    """
    Merge VCFs with bcftools merge, apply filters, and query to TSV.
    Keeps your hard-coded filter on ES:
        bcftools view -e 'FMT/ES == "."'
    """
    log_action(logger, "Starting VCF merge")

    query_fmt = build_query_format(fixed_fields, format_fields)

    # shell-safe quoting for VCF paths and bcftools executable
    vcf_str = " ".join(shlex.quote(str(v)) for v in vcf_list)
    out_tsv = str(out_tsv)
    out_tsv_q = shlex.quote(out_tsv)
    bcftools_q = shlex.quote(str(bcftools_exec))

    cmd = textwrap.dedent(
        f"""
        set -euo pipefail
        {bcftools_q} merge --threads {n_threads} -m none {vcf_str} | \\
        {bcftools_q} view --threads {n_threads} -e 'FMT/ES == "."' | \\
        {bcftools_q} view --threads {n_threads} --min-alleles 2 --max-alleles 2 | \\
        {bcftools_q} query --print-header -f '{query_fmt}' | \\
        sed 's/\\t$//' > {out_tsv_q}
        """
    )

    log_action(logger, "Executing bcftools merge")

    try:
        subprocess.run(cmd, shell=True, executable="/bin/bash", check=True)
        log_detected(logger, "VCF merge completed successfully")
    except subprocess.CalledProcessError as e:
        log_error(logger, f"VCF merge failed (returncode={e.returncode})")
        raise


# =========================================================
# 4️⃣ Count TSV variants
# =========================================================
def count_tsv_rows(tsv, logger):
    """
    Count variants (rows) in TSV excluding header.
    """
    tsv = str(tsv)
    log_action(logger, f"Counting TSV rows → {tsv}")

    tsv_q = shlex.quote(tsv)
    out = subprocess.check_output(f"wc -l {tsv_q}", shell=True, text=True)
    n = int(out.split()[0]) - 1  # minus header

    if n < 0:
        n = 0

    log_detected(logger, f"Total variants in TSV → {n:,}")
    return n


# =========================================================
# 5️⃣ Partial QC scan
# =========================================================
def partial_qc_scan(tsv, logger, n_rows=1_000_000):
    """
    Read first n_rows and detect empty columns.
    Returns: (columns, empty_cols)
    """
    tsv = str(tsv)
    log_action(logger, f"Partial QC scan → first {n_rows:,} rows")

    try:
        df = pl.read_csv(
            tsv,
            separator="\t",
            null_values=["", "."],
            truncate_ragged_lines=True,
            infer_schema_length=10000,
        ).head(n_rows)

    except Exception as e:
        log_error(logger, f"TSV read failed → {tsv}")
        log_error(logger, str(e))
        return [], ["QC_READ_FAILED"]

    # header cleanup (you said this is correct)
    try:
        # keep your behavior; rename safely
        old_cols = df.columns
        new_cols = [x.split("]")[-1] for x in old_cols]
        if new_cols != old_cols:
            df = df.rename(dict(zip(old_cols, new_cols)))
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
    """
    Extract sample IDs from column names like sample:AF, sample:ES, etc.
    (You confirmed this logic is correct for your header style.)
    """
    samples = set()
    for c in columns:
        if c not in fixed_fields and ":" in c:
            samples.add(c.split(":")[0])
    return samples


# =========================================================
# 7️⃣ Full wrapper (PRODUCTION; sys.exit(1) retained)
# =========================================================
def run_merge_vcf_to_tsv_pipeline(
    input_df,
    out_tsv,
    fixed_fields=None,
    format_fields=None,
    log_file=None,
    *,
    n_threads=4,
    qc_n_rows=1_000_000,
    bcftools_exec="bcftools",
    fail_fast_on_merge=True,
):
    """
    End-to-end:
      1) Count variants per input VCF
      2) Merge VCFs and query TSV
      3) Count TSV rows
      4) Partial QC scan (empty cols)
      5) Sample check

    Keeps sys.exit(1) behavior on failure (as you prefer).
    """
    # normalize input type (pandas → polars)
    input_df = normalize_input_df(input_df)

    # defaults
    if fixed_fields is None:
        fixed_fields = ["ID", "CHROM", "POS", "REF", "ALT"]

    if format_fields is None:
        format_fields = ["AF", "ES", "SE", "NEF", "LP", "SI", "EZ"]

    out_tsv = Path(out_tsv)

    # log file default
    if log_file is None:
        log_file = out_tsv.with_suffix(".log")
    else:
        log_file = Path(log_file)

    # setup logger early
    logger = setup_logger(log_file)

    # parent responsibility: ensure output dir exists
    try:
        out_tsv.parent.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        log_error(logger, f"Failed to create output directory: {out_tsv.parent}")
        log_error(logger, str(e))
        sys.exit(1)

    # cheap schema guard (prevents confusing crashes)
    required_cols = {"sumstat_vcf", "sample_id"}
    missing = required_cols - set(input_df.columns)
    if missing:
        log_error(logger, f"Missing required columns in input_df: {sorted(missing)}")
        sys.exit(1)

    # validate fixed_fields against build_query_format map (prevents KeyError)
    allowed_fixed = {"ID", "CHROM", "POS", "REF", "ALT"}
    bad_fixed = set(fixed_fields) - allowed_fixed
    if bad_fixed:
        log_error(logger, f"Invalid fixed_fields: {sorted(bad_fixed)}")
        sys.exit(1)

    qc_failures = []
    cols = []

    # ---------- STEP 1 ----------
    log_action(logger, "STEP 1 → Variant count per sample")

    for row in input_df.iter_rows(named=True):
        try:
            n = count_vcf_variants(row["sumstat_vcf"], logger, bcftools_exec=bcftools_exec)
            log_detected(logger, f"{row['sample_id']} → {n:,}")
        except Exception as e:
            qc_failures.append(f"Variant count failed → {row.get('sample_id', 'UNKNOWN')}")
            log_error(logger, str(e))

    # ---------- STEP 2 ----------
    log_action(logger, "STEP 2 → Merge VCFs and convert to TSV")

    try:
        merge_vcf_to_tsv(
            input_df["sumstat_vcf"].to_list(),
            out_tsv,
            fixed_fields,
            format_fields,
            logger,
            n_threads=n_threads,
            bcftools_exec=bcftools_exec,
        )
    except Exception as e:
        qc_failures.append("Merge failed")
        log_error(logger, str(e))

        if fail_fast_on_merge:
            log_action(logger, "FINAL QC SUMMARY")
            for f in qc_failures:
                log_error(logger, f)
            logger.error("Merging vcf and converting to TSV FAILED")
            logger.info("=" * 80)
            sys.exit(1)

    # ---------- STEP 3 ----------
    log_action(logger, "STEP 3 → Count TSV variants")

    try:
        count_tsv_rows(out_tsv, logger)
    except Exception as e:
        qc_failures.append("TSV count failed")
        log_error(logger, str(e))

    # ---------- STEP 4 ----------
    log_action(logger, "STEP 4 → Partial QC scan")

    try:
        cols, empty_cols = partial_qc_scan(out_tsv, logger, n_rows=qc_n_rows)
        for c in empty_cols:
            qc_failures.append(f"Empty column → {c}")
    except Exception as e:
        qc_failures.append("Partial QC failed")
        log_error(logger, str(e))

    # ---------- STEP 5 ----------
    log_action(logger, "STEP 5 → Sample presence check")

    try:
        tsv_samples = extract_samples(cols, fixed_fields)
        input_samples = set(input_df["sample_id"].to_list())

        missing_samples = input_samples - tsv_samples
        extra_samples = tsv_samples - input_samples

        if missing_samples:
            qc_failures.append(f"Missing samples → {sorted(missing_samples)}")

        if extra_samples:
            qc_failures.append(f"Unexpected samples → {sorted(extra_samples)}")

    except Exception as e:
        qc_failures.append("Sample extraction failed")
        log_error(logger, str(e))

    # ---------- FINAL ----------
    log_action(logger, "FINAL QC SUMMARY")

    if qc_failures:
        for f in qc_failures:
            log_error(logger, f)

        logger.error("Merging vcf and converting to TSV FAILED")
        logger.info("=" * 80)
        sys.exit(1)

    log_action(logger, "PIPELINE COMPLETED SUCCESSFULLY")
    return out_tsv, True


# input_df = pd.read_csv("/Users/JJOHN41/Downloads/kadoorie_biobank/vcf_files/pleio_input.txt", sep="\t")
# run_merge_vcf_to_tsv_pipeline(
#     input_df=input_df,
#     out_tsv="/Users/JJOHN41/Downloads/kadoorie_biobank/vcf_files/kadoorie_biobank_test.tsv",
#     fixed_fields=None,
#     format_fields=None,
#     log_file=None
# )