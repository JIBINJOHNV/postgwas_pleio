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
    debug: bool = False,   # ✅ NEW
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

    if len(file_merge_keys) != len(sumstat_merge_keys):
        raise ValueError("Merge keys mismatch")

    dbg(f"Files: {files}")
    dbg(f"File merge keys: {file_merge_keys}")
    dbg(f"Sumstat merge keys: {sumstat_merge_keys}")
    dbg(f"Suffix: {column_suffix}, Method: {method}, Final col: {final_col}")

    # ---------------------------
    # Aggregation function mapper
    # ---------------------------
    agg_map = {
        "mean": pl.mean_horizontal,
        "median": pl.median_horizontal,
        "min": pl.min_horizontal,
        "max": pl.max_horizontal,
        "sum": pl.sum_horizontal,
    }

    if method not in agg_map:
        raise ValueError(f"Invalid method: {method}")

    dfs = []
    per_file_cols = []

    for i, file in enumerate(files):
        dbg(f"\n--- Processing file {i+1}/{len(files)}: {file} ---")

        df = pl.read_csv(file, separator=sep, infer_schema_length=10000)
        print("infor file")
        print(df.head())
        print(df.columns)

        dbg(f"Columns in file: {df.columns[:20]}{'...' if len(df.columns) > 20 else ''}")
        dbg(f"Row count: {df.height}")

        # ---------------------------
        # Validate merge keys
        # ---------------------------
        for k in file_merge_keys:
            if k not in df.columns:
                raise ValueError(f"Missing key '{k}' in {file}")

        dbg(f"Merge keys validated")

        # ---------------------------
        # Detect suffix columns
        # ---------------------------
        suffix = column_suffix.lower()

        cols = [
            c for c in df.columns
            if c.lower().endswith(f":{suffix}") or c.lower() == suffix or c.lower().endswith(suffix)
        ]

        dbg(f"Detected columns for suffix '{column_suffix}': {cols}")

        if not cols:
            raise ValueError(f"No column ending with '{column_suffix}' in {file}")

        # ---------------------------
        # Validate extra columns
        # ---------------------------
        for col in extra_columns:
            if col not in df.columns:
                raise ValueError(f"Missing extra column '{col}' in {file}")

        if extra_columns:
            dbg(f"Extra columns: {extra_columns}")

        # ---------------------------
        # Select required columns
        # ---------------------------
        select_cols = list(dict.fromkeys(file_merge_keys + cols + extra_columns))
        df = df.select(select_cols)

        dbg(f"Selected columns: {select_cols}")

        # ---------------------------
        # Deduplicate keys
        # ---------------------------
        before = df.height
        df = df.unique(subset=file_merge_keys, keep="first")
        after = df.height

        if before != after:
            dbg(f"Deduplicated rows: {before} → {after}")

        # ---------------------------
        # Aggregate within file
        # ---------------------------
        temp_col = f"{final_col}__file{i}"

        agg_expr = agg_map[method]([pl.col(c) for c in cols])

        df = df.with_columns(agg_expr.alias(temp_col))

        dbg(f"Created temp aggregated column: {temp_col}")

        df = df.select(
            list(dict.fromkeys(file_merge_keys + [temp_col] + extra_columns))
        )

        dfs.append(df)
        per_file_cols.append(temp_col)

    # ---------------------------
    # Merge all files
    # ---------------------------
    dbg("\n--- Merging files ---")

    merged = dfs[0]
    for i, d in enumerate(dfs[1:], start=2):
        dbg(f"Merging file {i}")
        merged = merged.join(d, on=file_merge_keys, how="outer", suffix="_dup")

    dbg(f"Merged columns: {merged.columns}")

    # ---------------------------
    # Aggregate across files
    # ---------------------------
    dbg(f"Aggregating across files using columns: {per_file_cols}")

    exprs = [pl.col(c) for c in per_file_cols]

    merged = merged.with_columns(
        agg_map[method](exprs).alias(final_col)
    )

    # ---------------------------
    # Resolve duplicate extra columns
    # ---------------------------
    for col in extra_columns:
        candidates = [
            c for c in merged.columns
            if c == col or c.startswith(f"{col}_dup")
        ]

        if len(candidates) > 1:
            dbg(f"Resolving duplicates for column '{col}': {candidates}")

            merged = merged.with_columns(
                pl.coalesce([pl.col(c) for c in candidates]).alias(col)
            )
            merged = merged.drop([c for c in candidates if c != col])

    # ---------------------------
    # Drop intermediate columns
    # ---------------------------
    merged = merged.drop(per_file_cols)

    # ---------------------------
    # Final select
    # ---------------------------
    final_cols = list(dict.fromkeys(file_merge_keys + [final_col] + extra_columns))
    merged = merged.select(final_cols)

    dbg(f"Final selected columns: {final_cols}")

    # ---------------------------
    # Rename keys
    # ---------------------------
    rename_map = dict(zip(file_merge_keys, sumstat_merge_keys))
    merged = merged.rename(rename_map)

    dbg(f"Renamed columns: {rename_map}")
    dbg(f"Final output columns: {merged.columns}")
    dbg(f"Final row count: {merged.height}")

    return merged, files





    
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
    print("create_manifest_with_info")
    print(sumstat_file)
    print(output_folder) 
    
    def is_real_col(col):
        return col not in (None, "", "NA")

    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    df = pl.read_csv(sumstat_file, separator=sep, infer_schema_length=10000)
    print(df.head())
    print(df.columns)

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
    print(df.head())
    print(df.columns) 
    # ---------------------------------------------------
    # INFO MERGE
    # ---------------------------------------------------
    print("INFO MERGE")
    if info_files:
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
    print(df.head())
    print(df.columns)
    # ---------------------------------------------------
    # SAMPLE SIZE MERGE
    # ---------------------------------------------------
    print("SAMPLE SIZE MERGE")
    if sample_size_files:
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
        "imp_info_col": [final_info_col],
        "infofile": [infofile],
        "infocolumn": [infocolumn],
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

    

# from __future__ import annotations
# import polars as pl
# from pathlib import Path
# from typing import List


# def create_manifest_with_info(
#     sumstat_file: str,
#     sumstat_merge_keys: List[str],
#     info_merge_keys: List[str],
#     gwas_outputname: str,
#     resource_folder: Path | str,
#     output_folder: Path | str,
#     chr_col: str,
#     pos_col: str,
#     snp_id_col: str,
#     ea_col: str,
#     oa_col: str,
#     eaf_col: str,
#     beta_or_col: str,
#     se_col: str,
#     pval_col: str,
#     ncontrol_col: str,
#     ncase_col: str,
#     ncontrol: str = "NA",
#     ncase: str = "NA",
#     eaffile: str = "NA",
#     eafcolumn: str = "NA",
#     chr_pos_col: str = "NA",
#     imp_z_col: str = "NA",
#     liftover: str = "Yes",
#     info_files: List[str] | str | None = None,
#     info_suffix: str = "_info",
#     info_method: str = "mean",
#     final_info_col: str = "info_score",
#     extra_info_file_columns: List[str] | None = None,
#     manifest_df: pl.DataFrame | None = None,
#     sep: str = "\t",
# ):
#     output_folder = Path(output_folder)
#     output_folder.mkdir(parents=True, exist_ok=True)

#     df = pl.read_csv(sumstat_file, separator=sep, infer_schema_length=10000)

#     if extra_info_file_columns is None:
#         extra_info_file_columns = []

#     if info_files is not None and isinstance(info_files, str):
#         info_files = [info_files]

#     if len(info_merge_keys) != len(sumstat_merge_keys):
#         raise ValueError(
#             f"info_merge_keys and sumstat_merge_keys must have the same length. "
#             f"Got {len(info_merge_keys)} and {len(sumstat_merge_keys)}."
#         )

#     info_cols: List[str] = []
#     infofile_str = "NA"

#     def is_real_col(colname: str) -> bool:
#         return colname not in (None, "", "NA")

#     # ---------------------------------------------------
#     # Uppercase alleles only if columns are real and present
#     # ---------------------------------------------------
#     allele_exprs = []
#     if is_real_col(ea_col):
#         if ea_col not in df.columns:
#             raise ValueError(
#                 f"ea_col='{ea_col}' not found in sumstat file. "
#                 f"Available columns: {list(df.columns)}"
#             )
#         allele_exprs.append(pl.col(ea_col).cast(pl.Utf8).str.to_uppercase().alias(ea_col))

#     if is_real_col(oa_col):
#         if oa_col not in df.columns:
#             raise ValueError(
#                 f"oa_col='{oa_col}' not found in sumstat file. "
#                 f"Available columns: {list(df.columns)}"
#             )
#         allele_exprs.append(pl.col(oa_col).cast(pl.Utf8).str.to_uppercase().alias(oa_col))

#     if allele_exprs:
#         df = df.with_columns(allele_exprs)

#     # ---------------------------------------------------
#     # Process INFO files
#     # ---------------------------------------------------
#     if info_files:
#         info_dfs = []

#         for file in info_files:
#             info_df = pl.read_csv(file, separator=sep, infer_schema_length=10000)

#             for k in info_merge_keys:
#                 if k not in info_df.columns:
#                     raise ValueError(
#                         f"Merge key '{k}' missing in INFO file: {file}. "
#                         f"Available columns: {list(info_df.columns)}"
#                     )

#             detected_cols = [
#                 c for c in info_df.columns
#                 if c.lower().endswith(info_suffix.lower())
#             ]

#             if not detected_cols:
#                 raise ValueError(
#                     f"No column ending with '{info_suffix}' found in INFO file: {file}. "
#                     f"Available columns: {list(info_df.columns)}"
#                 )

#             info_col = detected_cols[0]

#             for col in extra_info_file_columns:
#                 if col not in info_df.columns:
#                     raise ValueError(
#                         f"Extra INFO column '{col}' missing in INFO file: {file}. "
#                         f"Available columns: {list(info_df.columns)}"
#                     )

#             selected_cols = list(dict.fromkeys(info_merge_keys + [info_col] + extra_info_file_columns))
#             info_df = info_df.select(selected_cols)

#             info_cols.append(info_col)
#             info_dfs.append(info_df)

#         merged_info = info_dfs[0]

#         for d in info_dfs[1:]:
#             merged_info = merged_info.join(
#                 d,
#                 on=info_merge_keys,
#                 how="outer",
#                 suffix="_dup"
#             )

#         exprs = [pl.col(c) for c in info_cols]

#         if info_method == "mean":
#             merged_info = merged_info.with_columns(pl.mean_horizontal(exprs).alias(final_info_col))
#         elif info_method == "median":
#             merged_info = merged_info.with_columns(pl.median_horizontal(exprs).alias(final_info_col))
#         elif info_method == "min":
#             merged_info = merged_info.with_columns(pl.min_horizontal(exprs).alias(final_info_col))
#         elif info_method == "max":
#             merged_info = merged_info.with_columns(pl.max_horizontal(exprs).alias(final_info_col))
#         else:
#             raise ValueError("info_method must be one of: mean, median, min, max")

#         # resolve duplicated extra columns more safely
#         for col in extra_info_file_columns:
#             candidates = [c for c in merged_info.columns if c == col or c.startswith(f"{col}_dup")]

#             if len(candidates) > 1:
#                 merged_info = merged_info.with_columns(
#                     pl.coalesce([pl.col(c) for c in candidates]).alias(col)
#                 )
#                 merged_info = merged_info.drop([c for c in candidates if c != col])

#         merged_info = merged_info.select(
#             list(dict.fromkeys(info_merge_keys + [final_info_col] + extra_info_file_columns))
#         )

#         rename_map = dict(zip(info_merge_keys, sumstat_merge_keys))
#         merged_info = merged_info.rename(rename_map)

#         missing_sumstat_keys = [k for k in sumstat_merge_keys if k not in df.columns]
#         if missing_sumstat_keys:
#             raise ValueError(
#                 f"Sumstat merge keys missing in sumstat file: {missing_sumstat_keys}. "
#                 f"Available columns: {list(df.columns)}"
#             )

#         df = df.join(merged_info, on=sumstat_merge_keys, how="left")
#         infofile_str = ",".join(info_files)

#     else:
#         if final_info_col not in df.columns:
#             df = df.with_columns(pl.lit(None).alias(final_info_col))

#     # ---------------------------------------------------
#     # Validate mandatory columns AFTER schema changes
#     # ---------------------------------------------------
#     base_column_map = {
#         "snp_id_col": snp_id_col,
#         "ea_col": ea_col,
#         "oa_col": oa_col,
#         "eaf_col": eaf_col,
#         "pval_col": pval_col,
#         "ncontrol_col": ncontrol_col,
#         "chr_col": chr_col,
#         "pos_col": pos_col,
#     }

#     missing_base = []
#     for param, col in base_column_map.items():
#         if not is_real_col(col):
#             continue
#         if col not in df.columns:
#             missing_base.append(f"{param} -> '{col}'")

#     if missing_base:
#         raise ValueError(
#             "Mandatory columns missing in sumstat after processing:\n"
#             + "\n".join(missing_base)
#             + f"\n\nAvailable columns:\n{list(df.columns)}"
#         )

#     # ---------------------------------------------------
#     # Validate effect column logic
#     # Require: (beta + se) OR z
#     # ---------------------------------------------------
#     has_beta = is_real_col(beta_or_col) and beta_or_col in df.columns
#     has_se = is_real_col(se_col) and se_col in df.columns
#     has_z = is_real_col(imp_z_col) and imp_z_col in df.columns

#     if not ((has_beta and has_se) or has_z):
#         raise ValueError(
#             "GWAS effect columns missing.\n"
#             f"Expected either (beta_or_col='{beta_or_col}' AND se_col='{se_col}') "
#             f"OR imp_z_col='{imp_z_col}'.\n"
#             f"Detected: beta={has_beta}, se={has_se}, z={has_z}\n"
#             f"Available columns: {list(df.columns)}"
#         )

#     # ---------------------------------------------------
#     # Chromosome sorting only when chr and pos are real columns
#     # ---------------------------------------------------
#     if is_real_col(chr_col) and is_real_col(pos_col):
#         if chr_col not in df.columns or pos_col not in df.columns:
#             raise ValueError(
#                 f"Cannot sort by chromosome/position because "
#                 f"chr_col='{chr_col}' or pos_col='{pos_col}' is missing."
#             )

#         df = (
#             df.with_columns(
#                 pl.when(pl.col(chr_col).cast(pl.Utf8).str.to_uppercase() == "X")
#                 .then(23)
#                 .when(pl.col(chr_col).cast(pl.Utf8).str.to_uppercase() == "Y")
#                 .then(24)
#                 .when(pl.col(chr_col).cast(pl.Utf8).str.to_uppercase() == "MT")
#                 .then(25)
#                 .otherwise(pl.col(chr_col).cast(pl.Int64))
#                 .alias("_chr_sort")
#             )
#             .sort(["_chr_sort", pos_col])
#             .drop("_chr_sort")
#         )

#     # ---------------------------------------------------
#     # Write updated sumstat
#     # ---------------------------------------------------
#     harmonisation_dir = output_folder / "harmonisation_input_files"
#     harmonisation_dir.mkdir(parents=True, exist_ok=True)

#     updated_file = str(harmonisation_dir / f"{gwas_outputname}_with_info.tsv")
#     df.write_csv(updated_file, separator="\t")

#     # ---------------------------------------------------
#     # Create manifest row
#     # ---------------------------------------------------
#     row = pl.DataFrame({
#         "sumstat_file": [updated_file],
#         "gwas_outputname": [gwas_outputname],
#         "chr_col": [chr_col],
#         "pos_col": [pos_col],
#         "snp_id_col": [snp_id_col],
#         "ea_col": [ea_col],
#         "oa_col": [oa_col],
#         "eaf_col": [eaf_col],
#         "beta_or_col": [beta_or_col],
#         "se_col": [se_col],
#         "imp_z_col": [imp_z_col],
#         "pval_col": [pval_col],
#         "ncontrol_col": [ncontrol_col],
#         "ncase_col": [ncase_col],
#         "ncontrol": [ncontrol],
#         "ncase": [ncase],
#         "imp_info_col": [final_info_col],
#         "infofile": [infofile_str],
#         "infocolumn": [final_info_col],
#         "eaffile": [eaffile],
#         "eafcolumn": [eafcolumn],
#         "liftover": [liftover],
#         "chr_pos_col": [chr_pos_col],
#         "resourse_folder": [str(resource_folder)],
#         "output_folder": [str(output_folder)],
#     })

#     if manifest_df is None:
#         manifest_df = row
#     else:
#         manifest_df = pl.concat([manifest_df, row], how="vertical_relaxed")

#     return df, manifest_df