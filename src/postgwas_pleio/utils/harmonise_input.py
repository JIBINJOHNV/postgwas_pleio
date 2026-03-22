
from __future__ import annotations
import polars as pl
from pathlib import Path
from typing import List, Tuple


def process_external_files(
    files: List[str] | str,
    file_merge_keys: List[str],
    sumstat_merge_keys: List[str],
    column_suffix: str,
    method: str = "mean",
    final_col: str = "merged_value",
    extra_columns: List[str] | None = None,
    sep: str = "\t",
    debug: bool = False,
) -> Tuple[pl.DataFrame, List[str]]:

    def dbg(msg):
        if debug:
            print(f"[DEBUG] {msg}", flush=True)

    if extra_columns is None:
        extra_columns = []
    if isinstance(files, str):
        files = [files]
    if not files:
        raise ValueError("No files provided")

    # ---------------------------
    # Aggregation function mapper
    # ---------------------------
    # FIX: Use a more compatible way to handle median if median_horizontal is missing
    agg_map = {
        "mean": lambda cols: pl.mean_horizontal(cols),
        "min": lambda cols: pl.min_horizontal(cols),
        "max": lambda cols: pl.max_horizontal(cols),
        "sum": lambda cols: pl.sum_horizontal(cols),
        "median": lambda cols: pl.concat_list(cols).list.median(), # Compatible median
    }

    if method not in agg_map:
        raise ValueError(f"Invalid method: {method}")

    dfs = []
    per_file_cols = []

    for i, file in enumerate(files):
        df = pl.read_csv(file, separator=sep, infer_schema_length=10000)
        df.columns = [c.strip() for c in df.columns]

        # Suffix detection (case-insensitive)
        s_low = column_suffix.lower()
        cols = [
            c for c in df.columns
            if (c.lower().endswith(f":{s_low}") or c.lower().endswith(s_low))
            and c not in file_merge_keys
        ]

        if not cols:
            raise ValueError(f"No columns matching suffix '{column_suffix}' in {file}")

        # Process and Deduplicate
        temp_col = f"{final_col}__file{i}"
        df = df.unique(subset=file_merge_keys)
        
        # Apply the aggregation from the map
        df = df.with_columns(agg_map[method]([pl.col(c) for c in cols]).alias(temp_col))

        df = df.select(list(dict.fromkeys(file_merge_keys + [temp_col] + extra_columns)))
        dfs.append(df)
        per_file_cols.append(temp_col)

    # Merge Phase
    merged = dfs[0]
    for i, d in enumerate(dfs[1:], start=1):
        merged = merged.join(d, on=file_merge_keys, how="outer", suffix="_dup")
        for key in file_merge_keys:
            dup = f"{key}_dup"
            if dup in merged.columns:
                merged = merged.with_columns(pl.coalesce([pl.col(key), pl.col(dup)]).alias(key))
                merged = merged.drop(dup)

    # Aggregate across files
    merged = merged.with_columns(
        agg_map[method]([pl.col(c) for c in per_file_cols]).alias(final_col)
    )

    # Rename keys (ID -> SNP)
    rename_map = dict(zip(file_merge_keys, sumstat_merge_keys))
    merged = merged.rename(rename_map)

    final_cols = list(dict.fromkeys(sumstat_merge_keys + [final_col] + extra_columns))
    return merged.select(final_cols), files


    
def create_manifest_with_info(
    sumstat_file: str,
    sumstat_merge_keys: List[str],
    info_merge_keys: List[str],
    gwas_outputname: str,
    resource_folder: Path | str,
    output_folder: Path | str,
    chr_col: str,
    pos_col: str,
    snp_id_col: str,
    ea_col: str,
    oa_col: str,
    eaf_col: str,
    beta_or_col: str,
    se_col: str,
    pval_col: str,
    ncontrol_col: str,
    ncase_col: str,
    ncontrol: str = "NA",
    ncase: str = "NA",
    eaffile: str = "NA",
    eafcolumn: str = "NA",
    infofile: str = "NA",
    infocolumn: str = "NA",
    chr_pos_col: str = "NA",
    imp_z_col: str = "NA",
    liftover: str = "Yes",

    # INFO
    info_files: List[str] | str | None = None,
    info_suffix: str = "_info",
    info_method: str = "mean",
    final_info_col: str = "info_score",
    extra_info_file_columns: List[str] | None = None,

    # SAMPLE SIZE
    sample_size_files: List[str] | str | None = None,
    sample_size_merge_keys: List[str] | None = None,   # ✅ NEW
    sample_size_suffix: str = "_n",
    sample_size_method: str = "mean",
    final_sample_size_col: str = "sample_size",
    extra_sample_size_file_columns: List[str] | None = None,

    # FREQUENCY
    freq_files: List[str] | str | None = None,
    freq_merge_keys: List[str] | None = None,          # ✅ NEW
    freq_suffix: str = "_eaf",
    freq_method: str = "mean",
    final_freq_col: str = "freq",
    extra_freq_file_columns: List[str] | None = None,

    manifest_df: pl.DataFrame | None = None,
    sep: str = "\t",
):
    def is_real_col(col):
        return col not in (None, "", "NA")

    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    df = pl.read_csv(sumstat_file, separator=sep, infer_schema_length=10000)

    # ---------------------------------------------------
    # Uppercase alleles
    # ---------------------------------------------------
    allele_exprs = []
    if is_real_col(ea_col) and ea_col in df.columns:
        allele_exprs.append(pl.col(ea_col).cast(pl.Utf8).str.to_uppercase().alias(ea_col))
    if is_real_col(oa_col) and oa_col in df.columns:
        allele_exprs.append(pl.col(oa_col).cast(pl.Utf8).str.to_uppercase().alias(oa_col))
    if allele_exprs:
        df = df.with_columns(allele_exprs)
    # ---------------------------------------------------
    # INFO MERGE
    # ---------------------------------------------------
    if info_files: 
        print("INFO MERGE")
        merged_info, info_files_list = process_external_files(
            files=info_files,
            file_merge_keys=info_merge_keys,
            sumstat_merge_keys=sumstat_merge_keys,
            column_suffix=info_suffix,
            method=info_method,
            final_col=final_info_col,
            extra_columns=extra_info_file_columns,
            sep=sep,
        )
        df = df.join(merged_info, on=sumstat_merge_keys, how="left")

        # update manifest values safely
        infofile = ",".join(info_files_list)
        infocolumn = final_info_col
    else:
        if final_info_col not in df.columns:
            df = df.with_columns(pl.lit(None).alias(final_info_col))
    # ---------------------------------------------------
    # SAMPLE SIZE MERGE
    # ---------------------------------------------------
    if sample_size_files:
        print("SAMPLE SIZE MERGE")
        ss_keys = sample_size_merge_keys or sumstat_merge_keys

        merged_n, _ = process_external_files(
            files=sample_size_files,
            file_merge_keys=ss_keys,
            sumstat_merge_keys=sumstat_merge_keys,
            column_suffix=sample_size_suffix,
            method=sample_size_method,
            final_col=final_sample_size_col,
            extra_columns=extra_sample_size_file_columns,
            sep=sep,
        )
        df = df.join(merged_n, on=sumstat_merge_keys, how="left")
        ncontrol_col = final_sample_size_col

    # ---------------------------------------------------
    # FREQUENCY MERGE
    # ---------------------------------------------------
    if freq_files:
        f_keys = freq_merge_keys or sumstat_merge_keys

        merged_freq, _ = process_external_files(
            files=freq_files,
            file_merge_keys=f_keys,
            sumstat_merge_keys=sumstat_merge_keys,
            column_suffix=freq_suffix,
            method=freq_method,
            final_col=final_freq_col,
            extra_columns=extra_freq_file_columns,
            sep=sep,
        )
        df = df.join(merged_freq, on=sumstat_merge_keys, how="left") 
        eaf_col = final_freq_col

    # ---------------------------------------------------
    # Validate required columns
    # ---------------------------------------------------
    required_cols = [snp_id_col, ea_col, oa_col, pval_col]
    for col in required_cols:
        if is_real_col(col) and col not in df.columns:
            raise ValueError(f"Missing column: {col}")

    # ---------------------------------------------------
    # Effect validation
    # ---------------------------------------------------
    has_beta = is_real_col(beta_or_col) and beta_or_col in df.columns
    has_se = is_real_col(se_col) and se_col in df.columns
    has_z = is_real_col(imp_z_col) and imp_z_col in df.columns

    if not ((has_beta and has_se) or has_z):
        raise ValueError("Need beta+se OR z")

    # ---------------------------------------------------
    # Sorting
    # ---------------------------------------------------
    if is_real_col(chr_col) and is_real_col(pos_col):
        df = (
            df.with_columns(
                pl.when(pl.col(chr_col).cast(pl.Utf8).str.to_uppercase() == "X").then(23)
                .when(pl.col(chr_col).cast(pl.Utf8).str.to_uppercase() == "Y").then(24)
                .otherwise(pl.col(chr_col).cast(pl.Int64))
                .alias("_chr_sort")
            )
            .sort(["_chr_sort", pos_col])
            .drop("_chr_sort")
        )

    # ---------------------------------------------------
    # Write output
    # ---------------------------------------------------
    harmonisation_dir = output_folder / "harmonisation_input_files"
    harmonisation_dir.mkdir(parents=True, exist_ok=True)

    updated_file = str(harmonisation_dir / f"{gwas_outputname}_with_info.tsv")
    df.write_csv(updated_file, separator="\t")

    # ---------------------------------------------------
    # Manifest (UNCHANGED STRUCTURE)
    # ---------------------------------------------------
    row = pl.DataFrame({
        "sumstat_file": [updated_file],
        "gwas_outputname": [gwas_outputname],
        "chr_col": [chr_col],
        "pos_col": [pos_col],
        "snp_id_col": [snp_id_col],
        "ea_col": [ea_col],
        "oa_col": [oa_col],
        "eaf_col": [eaf_col],
        "beta_or_col": [beta_or_col],
        "se_col": [se_col],
        "imp_z_col": [imp_z_col],
        "pval_col": [pval_col],
        "ncontrol_col": [ncontrol_col],
        "ncase_col": [ncase_col],
        "ncontrol": [ncontrol],
        "ncase": [ncase], 
        "imp_info_col": [infocolumn],
        "infofile": ["NA"],
        "infocolumn": ["NA"],
        "eaffile": [eaffile],
        "eafcolumn": [eafcolumn],
        "liftover": [liftover],
        "chr_pos_col": [chr_pos_col],
        "resourse_folder": [str(resource_folder)],
        "output_folder": [str(output_folder)],
    })

    if manifest_df is None:
        manifest_df = row
    else:
        manifest_df = pl.concat([manifest_df, row], how="vertical_relaxed")

    return df, manifest_df