from __future__ import annotations
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union
import polars as pl
from postgwas_pleio.utils.logging_utils import setup_logger, log_action, log_detected, log_warning, log_error
import math
import pandas as pd

PathLike = Union[str, Path]



def _convert_lp_to_raw_pvalue(df: pl.DataFrame, logger) -> pl.DataFrame:
    """
    Convert columns ending with ':LP' to raw p-values.
    Overwrites original column.

    Rules:
      - LP > 323 → set to 323 (numeric stability)
      - LP must be finite
      - Reports all adjustments to log + screen
    """
    MAX_LP = 323.0  # float64 safe boundary
    lp_cols = [c for c in df.columns if c.endswith(":LP")]
    if not lp_cols:
        return df
    log_action(logger, f"Converting {len(lp_cols)} LP columns to raw p-values")
    print(f"[INFO] Converting {len(lp_cols)} LP columns to raw p-values", flush=True)
    # -----------------------------
    # Detect very large LP
    # -----------------------------
    large_lp = df.select([
        (pl.col(c) > MAX_LP).sum().alias(c)
        for c in lp_cols
    ])
    total_large_lp = sum(large_lp.row(0))
    if total_large_lp > 0:
        msg = f"{total_large_lp:,} LP values exceed {MAX_LP} and will be capped"
        log_warning(logger, msg)
        print(f"[WARNING] {msg}", flush=True)
    # -----------------------------
    # Detect negative LP
    # -----------------------------
    negative_lp = df.select([
        (pl.col(c) < 0).sum().alias(c)
        for c in lp_cols
    ])
    total_negative_lp = sum(negative_lp.row(0))
    if total_negative_lp > 0:
        msg = f"{total_negative_lp:,} negative LP values detected (will produce p > 1)"
        log_warning(logger, msg)
        print(f"[WARNING] {msg}", flush=True)
    # -----------------------------
    # Apply cap + convert
    # -----------------------------
    df = df.with_columns([
        pl.when(pl.col(c).is_finite())
          .then(10.0 ** (-pl.col(c).clip_max(MAX_LP)))
          .otherwise(None)
          .alias(c)
        for c in lp_cols
    ])
    # -----------------------------
    # Post-check p-values
    # -----------------------------
    over_one = df.select([
        (pl.col(c) > 1.0).sum().alias(c)
        for c in lp_cols
    ])
    total_over_one = sum(over_one.row(0))
    below_zero = df.select([
        (pl.col(c) < 0).sum().alias(c)
        for c in lp_cols
    ])
    total_below_zero = sum(below_zero.row(0))
    if total_over_one > 0:
        msg = f"{total_over_one:,} raw p-values > 1 detected after conversion"
        log_warning(logger, msg)
        print(f"[WARNING] {msg}", flush=True)
    if total_below_zero > 0:
        msg = f"{total_below_zero:,} raw p-values < 0 detected after conversion"
        log_warning(logger, msg)
        print(f"[WARNING] {msg}", flush=True)
    if (
        total_large_lp == 0 and
        total_negative_lp == 0 and
        total_over_one == 0 and
        total_below_zero == 0
    ):
        log_action(logger, "All LP values converted successfully within numeric range")
        print("[INFO] All LP values converted successfully within numeric range", flush=True)
    return df


def _print_main(msg: str) -> None:
    print(msg, flush=True)

def _clean_header_columns(cols: Sequence[str]) -> List[str]:
    return [c.split("]")[-1] for c in cols]

def _drop_fully_null_columns(df: pl.DataFrame) -> pl.DataFrame:
    kept = [name for name in df.columns if df[name].null_count() < df.height]
    return df.select([pl.col(name) for name in kept])

def _ensure_id_column(df: pl.DataFrame, logger) -> pl.DataFrame:
    if "ID" not in df.columns:
        raise ValueError("ID column must exist in input TSV")
    required = ["CHROM", "POS", "REF", "ALT"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        msg = f"Cannot fill missing ID values. Missing columns: {missing}"
        log_error(logger, msg)
        raise ValueError(msg)
    null_count = df.select(pl.col("ID").is_null().sum()).item()
    if null_count > 0:
        log_action(logger, f"Filling {null_count} missing IDs using CHROM_POS_REF_ALT")
        df = df.with_columns(
            pl.when(pl.col("ID").is_null())
            .then(
                pl.format("{}_{}_{}_{}",
                          pl.col("CHROM"),
                          pl.col("POS"),
                          pl.col("REF"),
                          pl.col("ALT"))
            )
            .otherwise(pl.col("ID"))
            .alias("ID")
        )
    return df

def _enforce_unique_id(df: pl.DataFrame, logger) -> pl.DataFrame:
    required = ["CHROM", "POS", "REF", "ALT"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Cannot create unique ID. Missing columns: {missing}")
    log_action(logger, "Recreating ID using CHROM_POS_REF_ALT")
    if "ID" in df.columns:
        df = df.drop("ID")
    df = df.with_columns(
        pl.format("{}_{}_{}_{}",
                  pl.col("CHROM"),
                  pl.col("POS"),
                  pl.col("REF"),
                  pl.col("ALT")
        ).alias("ID")
    )
    return df

def _remove_duplicate_ids(df: pl.DataFrame, logger) -> pl.DataFrame:
    if df.height==0:
        return df
    log_action(logger,"Removing duplicate IDs")
    return df.unique(subset=["ID"],keep="first")

def _apply_header_map(df: pl.DataFrame, header_map: Optional[Dict[str,str]], logger) -> pl.DataFrame:
    if not header_map:
        return df
    existing = {k:v for k,v in header_map.items() if k in df.columns}
    if not existing:
        log_warning(logger,"header_map provided but no columns matched")
        return df
    df = df.rename(existing)
    return df

def _transform_retained_columns(
    retained_cols,
    header_map,
    suffix_replace,
):
    transformed = []
    for col in retained_cols:
        if header_map and col in header_map:
            col = header_map[col]
        if suffix_replace:
            for old_suf, new_suf in suffix_replace.items():
                cleaned_old = old_suf.lstrip(":_")
                cleaned_new = new_suf.lstrip(":_")
                if col == cleaned_old:
                    col = cleaned_new
                    break
        transformed.append(col)
    return transformed

def _replace_suffixes(df: pl.DataFrame, 
        suffix_replace: Optional[Dict[str,str]], logger) -> pl.DataFrame:
    if not suffix_replace:
        return df
    new_cols=[]
    changed=0
    for c in df.columns:
        new_c=c
        for old,new in suffix_replace.items():
            if new_c.endswith(old):
                new_c=new_c[:-len(old)]+new
                changed+=1
                break  # prevent multiple suffix changes
        new_cols.append(new_c)
    if changed:
        log_action(logger,f"Suffix replacement updated {changed} columns")
    df.columns=new_cols
    return df

def _prepare_base_dataframe(tsv_file: Path, logger) -> pl.DataFrame:
    log_action(logger,"Reading TSV")
    df=pl.read_csv(tsv_file,separator="\t",null_values=["","."],truncate_ragged_lines=True,infer_schema_length=10000)
    df.columns=_clean_header_columns(df.columns)
    df=_drop_fully_null_columns(df)
    df=_ensure_id_column(df,logger)
    return df

# def _filter_by_manifest(df: pl.DataFrame,input_manifest_df: pl.DataFrame,logger)->pl.DataFrame:
#     if "sample_id" not in input_manifest_df.columns:
#         raise ValueError("input_manifest_df must contain 'sample_id'")
#     samples=input_manifest_df["sample_id"].drop_nulls().unique().to_list()
#     common_cols=[c for c in df.columns if not any(c.startswith(f"{s}:") for s in samples)]
#     sample_cols=[c for c in df.columns if any(c.startswith(f"{s}:") for s in samples)]
#     return df.select(common_cols+sample_cols)

def _filter_by_manifest(
    df: pl.DataFrame,
    input_manifest_df: pl.DataFrame,
    logger
) -> pl.DataFrame:
    """
    Keep:
      - All base columns (columns without sample prefix)
      - Only sample-prefixed columns matching manifest sample_id

    Sample-prefixed columns are assumed to follow pattern:
        sample_id:METRIC
    """
    if "sample_id" not in input_manifest_df.columns:
        raise ValueError("input_manifest_df must contain 'sample_id'")
    # Collect unique sample IDs
    samples = (
        input_manifest_df["sample_id"]
        .drop_nulls()
        .unique()
        .to_list()
    )
    if not samples:
        log_warning(logger, "Manifest contains no valid sample_id entries")
        return df
    # Identify sample-prefixed columns
    allowed_sample_cols = [
        c for c in df.columns
        if any(c.startswith(f"{s}:") for s in samples)
    ]
    # Base columns = no colon in name
    base_cols = [c for c in df.columns if ":" not in c]
    kept_cols = base_cols + allowed_sample_cols
    removed_cols = len(df.columns) - len(kept_cols)
    log_action(
        logger,
        f"Manifest filtering → keeping {len(allowed_sample_cols)} "
        f"sample-prefixed columns; removed {removed_cols} columns"
    )
    return df.select(kept_cols)

def run_postprocess_joined_samples_from_df(
    df: pl.DataFrame,
    *,
    joined_out_tsv: PathLike,
    joined_remove_duplicate_ids: bool,
    joined_create_unique_id: bool,
    joined_sample_retained_cols: Optional[Sequence[str]],
    joined_header_map: Optional[Dict[str,str]],
    joined_suffix_replace: Optional[Dict[str,str]],
    joined_convert_to_raw_p_value: bool,
    logger,
) -> Path:
    _print_main("Processing JOINED samples")
    log_action(logger, "JOINED MODE → Start")
    df2 = df.clone()
    # 1️⃣ Create unique ID only if requested (does NOT overwrite existing non-null IDs)
    if joined_create_unique_id:
        log_action(logger, "JOINED MODE → Enforcing unique IDs")
        df2 = _enforce_unique_id(df2,logger)
    # 2️⃣ Convert LP to raw p-value if requested
    if joined_convert_to_raw_p_value:
        df2 = _convert_lp_to_raw_pvalue(df2, logger)
    # 3️⃣ Remove duplicate rows if requested
    if joined_remove_duplicate_ids:
        log_action(logger, "JOINED MODE → Removing duplicate IDs")
        df2 = _remove_duplicate_ids(df2, logger)
    # 4️⃣ Apply header renaming
    if joined_header_map:
        df2 = _apply_header_map(df2, joined_header_map, logger)
    # 5️⃣ Apply suffix replacement
    if joined_suffix_replace:
        df2 = _replace_suffixes(df2, joined_suffix_replace, logger)
    # 6️⃣ Apply retention AFTER renaming
    if joined_sample_retained_cols:
        transformed_cols = _transform_retained_columns(
            joined_sample_retained_cols,
            joined_header_map,
            joined_suffix_replace,
        )
        retained_cols = []
        for col in transformed_cols:
            # Normalize (remove leading dot or colon artifacts)
            normalized_col = col.lstrip(".:_")
            # 1️⃣ Direct match
            if col in df2.columns:
                retained_cols.append(col)
                continue
            # 2️⃣ Suffix match (support :, _, and . delimiters)
            matches = [
                c for c in df2.columns
                if (
                    c.endswith(f":{normalized_col}") or
                    c.endswith(f"_{normalized_col}") or
                    c.endswith(f".{normalized_col}")
                )
            ]
            retained_cols.extend(matches)
        if not retained_cols:
            log_warning(logger, "No retained columns matched after renaming")
        df2 = df2.select(retained_cols)
    # 6️⃣ Ensure ID is first column
    id_column = None
    if "ID" in df2.columns:
        id_column = "ID"
    elif joined_header_map and "ID" in joined_header_map:
        mapped_id = joined_header_map["ID"]
        if mapped_id in df2.columns:
            id_column = mapped_id
    if id_column:
        df2 = df2.select([id_column] + [c for c in df2.columns if c != id_column])
    # 7️⃣ Write output
    out_path = Path(joined_out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df2.write_csv(out_path, separator="\t")
    log_detected(logger, f"Joined output written → {out_path}")
    log_detected(logger, f"Final shape → rows={df2.height}, cols={df2.width}")
    _print_main(f"Joined output written: {out_path}")
    return out_path

def run_postprocess_sample_specific_from_df(
    df: pl.DataFrame,
    *,
    sample_specific_directory: PathLike,
    input_manifest_df: pl.DataFrame,
    sample_remove_duplicate_ids: bool,
    sample_create_unique_id: bool,
    sample_specific_retained_cols: Optional[Sequence[str]],
    sample_header_map: Optional[Dict[str, str]],
    sample_suffix_replace: Optional[Dict[str, str]],
    sample_convert_to_raw_p_value: bool,
    logger,
    output_map="matrix_tsv",
) -> Path:
    outdir = Path(sample_specific_directory)
    outdir.mkdir(parents=True, exist_ok=True)
    samples = (
        input_manifest_df["sample_id"]
        .drop_nulls()
        .unique()
        .to_list()
    )
    log_action(logger, "SAMPLE-SPECIFIC MODE → Start")
    sample_output_map = []
    for s in samples:
        cols = ["ID", "CHROM", "POS", "REF", "ALT"] + [
            c for c in df.columns if c.startswith(f"{s}:")
        ]
        if len(cols) == 5:
            log_warning(logger, f"No metric columns found for sample {s}")
            continue
        df_sample = df.select(cols)
        # 1️⃣ Unique ID
        if sample_create_unique_id:
            log_action(logger, f"{s} → Enforcing unique IDs")
            df_sample = _enforce_unique_id(df_sample, logger)
        # 3️⃣ Convert LP
        if sample_convert_to_raw_p_value:
            df_sample = _convert_lp_to_raw_pvalue(df_sample, logger)
        # 2️⃣ Remove duplicates
        if sample_remove_duplicate_ids:
            log_action(logger, f"{s} → Removing duplicate IDs")
            df_sample = _remove_duplicate_ids(df_sample, logger)
        # 4️⃣ Header rename
        if sample_header_map:
            df_sample = _apply_header_map(df_sample, sample_header_map, logger)
        # 5️⃣ Suffix replacement
        if sample_suffix_replace:
            df_sample = _replace_suffixes(df_sample, sample_suffix_replace, logger)
        # 6️⃣ Retention
        # if sample_specific_retained_cols:
        #     transformed_cols = _transform_retained_columns(
        #         sample_specific_retained_cols,
        #         sample_header_map,
        #         sample_suffix_replace,
        #     )
        #     retained_cols = []
        #     seen = set()
        #     for col in transformed_cols:
        #         normalized_col = col.lstrip(".")
        #         # 1️⃣ Direct match
        #         if col in df_sample.columns and col not in seen:
        #             retained_cols.append(col)
        #             seen.add(col)
        #             continue
        #         # 2️⃣ Suffix match (support :, _, and . delimiters)
        #         matches = [
        #             c for c in df_sample.columns
        #             if (
        #                 c.endswith(f":{normalized_col}") or
        #                 c.endswith(f"_{normalized_col}") or
        #                 c.endswith(f".{normalized_col}")
        #             )
        #         ]
        #         for m in matches:
        #             if m not in seen:
        #                 retained_cols.append(m)
        #                 seen.add(m)
        #     if not retained_cols:
        #         log_warning(logger, f"{s} → No retained columns matched after renaming")
        #     else:
        #         df_sample = df_sample.select(retained_cols)
        # 6️⃣ Retention (preserve order from sample_specific_retained_cols)
        if sample_specific_retained_cols:
            transformed_cols = _transform_retained_columns(
                sample_specific_retained_cols,
                sample_header_map,
                sample_suffix_replace,
            )
            retained_cols = []
            seen = set()
            for requested_col in transformed_cols:
                normalized_col = requested_col.lstrip(".")
                # 1️⃣ Direct match
                if requested_col in df_sample.columns and requested_col not in seen:
                    retained_cols.append(requested_col)
                    seen.add(requested_col)
                    continue
                # 2️⃣ Suffix matches
                matches = [
                    c for c in df_sample.columns
                    if (
                        c.endswith(f":{normalized_col}") or
                        c.endswith(f"_{normalized_col}") or
                        c.endswith(f".{normalized_col}")
                    )
                ]
                # IMPORTANT: preserve df column order but attach them
                # under the requested column position
                for m in matches:
                    if m not in seen:
                        retained_cols.append(m)
                        seen.add(m)
            if not retained_cols:
                log_warning(logger, f"{s} → No retained columns matched after renaming")
            else:
                # enforce exact order
                retained_cols = [c for c in retained_cols if c in df_sample.columns]
                df_sample = df_sample.select(retained_cols)
        # 7️⃣ Strip sample prefix ONLY at the very end
        rename_map = {
            c: c.replace(f"{s}:", "")
            for c in df_sample.columns
            if c.startswith(f"{s}:")
        }
        if rename_map:
            df_sample = df_sample.rename(rename_map)
        # 8️⃣ Write file
        out_file = outdir / f"{s}_matrix.tsv"
        df_sample.write_csv(out_file, separator="\t")
        sample_output_map.append(
            {"sample_id": s, output_map : str(out_file)}
        )
        log_detected(logger, f"{s} → Wrote {out_file}")
    output_df = pl.DataFrame(sample_output_map)
    input_manifest_df = input_manifest_df.join(
        output_df,
        on="sample_id",
        how="left",
    )
    return input_manifest_df

def run_mastertsv_to_toolformat_pipeline(*,mode: Sequence[str],
        tsv_file: PathLike,
        input_manifest_df: pl.DataFrame,
        joined_out_tsv: Optional[PathLike]=None,
        sample_specific_directory: Optional[PathLike]=None,
        sample_specific_ldsc_directory: Optional[PathLike]=None,
        joined_remove_duplicate_ids: bool=False,
        joined_create_unique_id: bool=False,
        joined_sample_retained_cols: Optional[Sequence[str]]=None,
        joined_header_map: Optional[Dict[str,str]]=None,
        joined_suffix_replace: Optional[Dict[str,str]]=None,
        joined_convert_to_raw_p_value: bool = False,
        sample_remove_duplicate_ids: bool=False,
        sample_create_unique_id: bool=False,
        sample_specific_retained_cols: Optional[Sequence[str]]=None,
        sample_header_map: Optional[Dict[str,str]]=None,
        sample_suffix_replace: Optional[Dict[str,str]]=None,
        sample_convert_to_raw_p_value: bool = False,
        sample_ldsc_remove_duplicate_ids: bool=False,
        sample_ldsc_create_unique_id: bool=False,
        sample_ldsc_specific_retained_cols: Optional[Sequence[str]]=None,
        sample_ldsc_header_map: Optional[Dict[str,str]]=None,
        sample_ldsc_suffix_replace: Optional[Dict[str,str]]=None,
        sample__ldscconvert_to_raw_p_value: bool = False,
        log_file: Optional[PathLike]=None)->Dict[str,Path]:
    
    if isinstance(input_manifest_df, pd.DataFrame):
        input_manifest_df = pl.from_pandas(input_manifest_df)
    elif not isinstance(input_manifest_df, pl.DataFrame):
        raise TypeError("input_manifest_df must be pandas or polars DataFrame")
    mode_set={m.strip().lower() for m in mode}
    valid={"joined_samples","sample_specific","sample_specific_ldsc"}
    
    if not mode_set or not mode_set.issubset(valid):
        raise ValueError(f"mode must be subset of {valid}")
    tsv_path=Path(tsv_file)
    if log_file is None:
        log_file=tsv_path.with_suffix(".postprocess.log")
    logger=setup_logger(Path(log_file))
    log_action(logger,"POSTPROCESS PIPELINE STARTED")
    df=_prepare_base_dataframe(tsv_path,logger)
    df=_filter_by_manifest(df,input_manifest_df,logger)
    results={}
    
    if "joined_samples" in mode_set:
        if joined_out_tsv is None:
            raise ValueError("joined_out_tsv required")
        results["joined"]=run_postprocess_joined_samples_from_df(
            df,
            joined_out_tsv=joined_out_tsv,
            joined_remove_duplicate_ids=joined_remove_duplicate_ids,
            joined_create_unique_id=joined_create_unique_id,
            joined_sample_retained_cols=joined_sample_retained_cols,
            joined_header_map=joined_header_map,    
            joined_suffix_replace=joined_suffix_replace,
            joined_convert_to_raw_p_value=joined_convert_to_raw_p_value,
            logger=logger,
        )
    if "sample_specific" in mode_set:
        if sample_specific_directory is None:
            raise ValueError("sample_specific_directory required")
        results["sample_specific"]=run_postprocess_sample_specific_from_df(
            df,
            sample_specific_directory=sample_specific_directory,
            input_manifest_df=input_manifest_df,
            sample_remove_duplicate_ids=sample_remove_duplicate_ids,
            sample_create_unique_id=sample_create_unique_id,
            sample_specific_retained_cols=sample_specific_retained_cols,
            sample_header_map=sample_header_map,
            sample_suffix_replace=sample_suffix_replace,
            sample_convert_to_raw_p_value=sample_convert_to_raw_p_value,
            logger=logger,
            output_map="matrix_tsv",
        )
    if "sample_specific_ldsc" in mode_set:
        if sample_specific_ldsc_directory is None:
            raise ValueError("sample_specific_directory required")
        results["sample_specific_ldsc"]=run_postprocess_sample_specific_from_df(
            df,
            sample_specific_directory=sample_specific_ldsc_directory,
            input_manifest_df=input_manifest_df,
            sample_remove_duplicate_ids=sample_ldsc_remove_duplicate_ids,
            sample_create_unique_id=sample_ldsc_create_unique_id,
            sample_specific_retained_cols=sample_ldsc_specific_retained_cols,
            sample_header_map=sample_ldsc_header_map,
            sample_suffix_replace=sample_ldsc_suffix_replace,
            sample_convert_to_raw_p_value=sample__ldscconvert_to_raw_p_value,
            logger=logger,
            output_map="ldsc_matrix_tsv"
        )
    log_action(logger,"POSTPROCESS PIPELINE COMPLETED")
    return results