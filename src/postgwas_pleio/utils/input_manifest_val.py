# =========================================================
# QC REPORT CLASS
# =========================================================
import pandas as pd
import numpy as np
import os
import gzip
import sys
from datetime import datetime


# =========================================================
# QC REPORT CLASS (PRODUCTION)
# =========================================================
class QCReport:

    def __init__(self, log_path):
        self.log_path = log_path
        self.lines = []
        self.errors = []
        self.warnings = []

        self._write("\n==== INPUT MANIFEST VALIDATION REPORT ====")
        self._write(f"Timestamp: {datetime.now()}\n")

    def _write(self, msg):
        self.lines.append(msg)

    def step(self, name):
        self._write(f"\n===== STEP: {name} =====")

    def action(self, msg):
        self._write(f"ACTION: {msg}")

    def detect(self, msg):
        self._write(f"DETECTED: {msg}")

    def correct(self, msg):
        self._write(f"CORRECTION: {msg}")

    def summary(self, msg):
        self._write(f"SUMMARY: {msg}")

    def warn(self, msg):
        self.warnings.append(msg)
        self._write(f"WARNING: {msg}")
        print(f"⚠ {msg}")

    def error(self, msg):
        self.errors.append(msg)
        self._write(f"ERROR: {msg}")
        print(f"❌ {msg}")

    def fail(self, msg):
        self.error(msg)
        self.finalize()
        print(f"\n❌ Validation failed — see log\n")
        sys.exit(1)

    def finalize(self):
        self._write("\n===== FINAL SUMMARY =====")
        self._write(f"Total warnings: {len(self.warnings)}")
        self._write(f"Total errors: {len(self.errors)}")

        with open(self.log_path, "w") as f:
            f.write("\n".join(self.lines))

        print(f"\n📝 Validation log saved → {self.log_path}")


# =========================================================
# VALIDATION FUNCTION (FULL PRODUCTION)
# =========================================================
def validate_manifest(manifest_path, output_dir):

    os.makedirs(output_dir, exist_ok=True)
    report = QCReport(os.path.join(output_dir, "input_manifest_validation.log"))

    # -----------------------------------------------------
    # Load manifest
    # -----------------------------------------------------
    report.step("Load input file manifest")
    report.action("Detecting delimiter and loading file")

    delimiters = ["\t", ",", ";", "|"]
    df = None
    used = None

    for d in delimiters:
        try:
            tmp = pd.read_csv(manifest_path, sep=d, nrows=5)
            if tmp.shape[1] > 1:
                df = pd.read_csv(manifest_path, sep=d)
                used = d
                break
        except Exception:
            continue

    if df is None:
        report.detect("Delimiter detection failed")
        report.fail("Could not detect manifest delimiter")

    delim_display = {"\t": "\\t", ",": ",", ";": ";", "|": "|"}[used]
    report.detect(f"Delimiter detected: {delim_display}")
    report.correct("Manifest loaded successfully")
    report.summary(f"{df.shape[0]} rows, {df.shape[1]} columns loaded")

    # -----------------------------------------------------
    # Column normalization
    # -----------------------------------------------------
    report.step("Column normalization")
    report.action("Trimming whitespace from column names")

    old = df.columns.tolist()
    df.columns = [c.strip() for c in df.columns]

    if old != df.columns.tolist():
        report.detect(f"Whitespace found in columns: {old}")
        report.correct("Column whitespace removed")
    else:
        report.detect("No column whitespace detected")
        report.correct("No correction required")

    report.summary(f"{len(df.columns)} columns present")

    # -----------------------------------------------------
    # Required columns
    # -----------------------------------------------------
    report.step("Required column validation")
    report.action("Checking required columns")

    required = ["sumstat_vcf", "TYPE", "SPREV", "PPREV", "sample_id"]
    missing = [c for c in required if c not in df.columns]

    if missing:
        report.detect(f"Missing columns: {missing}")
        report.fail("Required columns missing")

    report.detect("All required columns present")
    report.correct("No correction required")
    report.summary("Required column validation complete")

    # -----------------------------------------------------
    # Value trimming
    # -----------------------------------------------------
    report.step("Value trimming")
    report.action("Trimming whitespace in string columns")

    trim_counts = {}
    for col in df.select_dtypes(include="object"):
        before = df[col].copy()
        df[col] = df[col].astype(str).str.strip()
        trim_counts[col] = int((before != df[col]).sum())

    report.detect(f"Trim counts: {trim_counts}")
    report.correct("Whitespace removed where needed")
    report.summary("Value trimming complete")

    # -----------------------------------------------------
    # Row completeness
    # -----------------------------------------------------
    report.step("Row completeness")
    report.action("Checking missing key fields")

    mask = df[["sumstat_vcf","sample_id","TYPE"]].isna().any(axis=1)

    if mask.any():
        report.detect(f"Incomplete rows: {df.loc[mask,'sample_id'].tolist()}")
        report.fail("Incomplete rows present")

    report.detect("No incomplete rows detected")
    report.correct("No correction required")
    report.summary("Row completeness validation complete")

    # -----------------------------------------------------
    # Duplicate sample_id
    # -----------------------------------------------------
    report.step("Duplicate sample_id check")
    report.action("Checking uniqueness")

    dup = df["sample_id"][df["sample_id"].duplicated()].tolist()

    if dup:
        report.detect(f"Duplicate sample_id: {dup}")
        report.fail("Duplicate sample_id detected")

    report.detect("No duplicate sample_id detected")
    report.correct("No correction required")
    report.summary("sample_id uniqueness confirmed")

    # -----------------------------------------------------
    # Duplicate VCF logic
    # -----------------------------------------------------
    report.step("Duplicate VCF path check")
    report.action("Checking VCF reuse across sample_id")

    vcfd = df.groupby("sumstat_vcf")["sample_id"].nunique()
    bad = vcfd[vcfd>1].index.tolist()

    if bad:
        report.detect(f"VCF reused across sample_id: {bad}")
        report.fail("VCF reuse detected")

    report.detect("No duplicate VCF misuse detected")
    report.correct("No correction required")
    report.summary("VCF path uniqueness confirmed")

    # -----------------------------------------------------
    # TYPE normalization
    # -----------------------------------------------------
    report.step("TYPE normalization")
    report.action("Lowercasing TYPE values")

    before = df["TYPE"].unique().tolist()
    df["TYPE"] = df["TYPE"].str.lower()
    after = df["TYPE"].unique().tolist()

    report.detect(f"TYPE before: {before}")
    report.correct(f"TYPE after: {after}")
    report.summary("TYPE normalization complete")

    # -----------------------------------------------------
    # Allowed TYPE
    # -----------------------------------------------------
    report.step("Allowed TYPE validation")
    report.action("Checking allowed TYPE values")

    allowed = {"binary","quantitative"}
    bad = df.loc[~df["TYPE"].isin(allowed),"sample_id"].tolist()

    if bad:
        report.detect(f"Invalid TYPE rows: {bad}")
        report.fail("Invalid TYPE detected")

    report.detect("All TYPE values valid")
    report.correct("No correction required")
    report.summary("TYPE validation complete")

    # -----------------------------------------------------
    # Prevalence rules
    # -----------------------------------------------------
    report.step("Prevalence validation")
    report.action("Applying prevalence rules")

    non_binary = df["TYPE"]!="binary"
    if non_binary.any():
        report.detect(f"Non-binary rows: {df.loc[non_binary,'sample_id'].tolist()}")
        df.loc[non_binary,["SPREV","PPREV"]] = np.nan
        report.correct("Set SPREV/PPREV to NaN")
    else:
        report.detect("No non-binary rows detected")
        report.correct("No correction required")

    binary = df["TYPE"]=="binary"

    wrong_sprev = df.loc[binary & (df["SPREV"].isna()), "sample_id"].tolist()
    sprev_half = df.loc[binary & (df["SPREV"] == 0.5), "sample_id"].tolist()

    if wrong_sprev:
        report.warn(f"Missing SPREV for binary traits: {wrong_sprev}")
        report.warn(
            "SPREV will be estimated from the GWAS summary statistics as "
            "median(Ncases) / median(total sample size)."
        )
    elif sprev_half:
        report.warn(f"SPREV provided as 0.5 for these binary traits: {sprev_half}")
        report.warn(
            "This assumes equal case and control sample sizes. "
            "If that is not correct, please provide the true SPREV."
        )
    else:
        report.detect("SPREV is already available for all binary traits.")

    miss = df.loc[binary & df["PPREV"].isna(),"sample_id"].tolist()
    if miss:
        report.warn(f"Missing PPREV for binary: {miss}")

    bad = df.loc[binary & ((df["PPREV"]<=0)|(df["PPREV"]>=1)),"sample_id"].tolist()
    if bad:
        report.detect(f"PPREV outside range: {bad}")
        report.fail("PPREV outside (0,1)")

    report.summary("Prevalence validation complete")

    # -----------------------------------------------------
    # File existence
    # -----------------------------------------------------
    report.step("File existence")
    report.action("Checking file existence")

    missing = [p for p in df["sumstat_vcf"] if not os.path.exists(p)]

    if missing:
        report.detect(f"Missing files: {missing}")
        report.fail("Missing VCF files")

    report.detect("All files exist")
    report.correct("No correction required")
    report.summary("File existence validation complete")

    # -----------------------------------------------------
    # Gzip integrity
    # -----------------------------------------------------
    report.step("Gzip integrity")
    report.action("Checking gzip integrity")

    bad=[]
    for f in df["sumstat_vcf"]:
        if f.endswith(".gz"):
            try:
                with gzip.open(f) as g:
                    g.read(1)
            except:
                bad.append(f)

    if bad:
        report.detect(f"Corrupted gzip: {bad}")
        report.fail("Gzip integrity failed")

    report.detect("All gzip files valid")
    report.correct("No correction required")
    report.summary("Gzip integrity validation complete")

    # -----------------------------------------------------
    # Save snapshot
    # -----------------------------------------------------
    report.step("Save validated manifest")
    report.action("Writing validated snapshot")

    out=os.path.join(output_dir,"validated_input_manifest.tsv")
    df.to_csv(out,sep="\t",index=False)

    report.detect(f"Snapshot saved: {out}")
    report.correct("Validated manifest created")
    report.summary("Snapshot saved")

    report.finalize()
    return {
        "manifest_df": df,
        "manifest_path": out
    }


# validate_manifest(manifest_path="/Users/JJOHN41/Downloads/kadoorie_biobank/vcf_files/pleio_input.tsv", 
#                  output_dir="/Users/JJOHN41/Downloads/kadoorie_biobank/vcf_files/")