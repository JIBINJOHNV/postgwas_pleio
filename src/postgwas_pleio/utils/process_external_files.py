
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



# -------------------------------------------------------
# GENERIC MERGE FUNCTION (REQUIRED)
# -------------------------------------------------------

def process_external_files(
    files: List[str] | str,
    file_merge_keys: List[str],
    sumstat_merge_keys: List[str],
    column_suffix: str,
    method: str = "mean",
    final_col: str = "merged_value",
    extra_columns: List[str] | None = None,
    sep: str = "\t",
) -> Tuple[pl.DataFrame, List[str]]:

    if extra_columns is None:
        extra_columns = []

    if isinstance(files, str):
        files = [files]

    if len(file_merge_keys) != len(sumstat_merge_keys):
        raise ValueError("Merge keys mismatch")

    dfs = []
    detected_cols = []

    # -------------------------------------------------------
    # READ + EXPAND ALL MATCHING COLUMNS
    # -------------------------------------------------------
    for file_idx, file in enumerate(files):
        df = pl.read_csv(file, separator=sep, infer_schema_length=10000)

        # Validate merge keys
        for k in file_merge_keys:
            if k not in df.columns:
                raise ValueError(f"Missing key '{k}' in {file}")

        # Detect ALL matching columns
        cols = [c for c in df.columns if c.lower().endswith(column_suffix.lower())]
        if not cols:
            raise ValueError(f"No column ending with '{column_suffix}' in {file}")

        # Validate extra columns
        for col in extra_columns:
            if col not in df.columns:
                raise ValueError(f"Missing extra column '{col}' in {file}")

        # -------------------------------------------------------
        # Expand each detected column as separate dataset
        # -------------------------------------------------------
        for col_idx, target_col in enumerate(cols):

            tmp = df.select(
                list(dict.fromkeys(file_merge_keys + [target_col] + extra_columns))
            )

            # Ensure numeric (safe cast)
            tmp = tmp.with_columns(
                pl.col(target_col).cast(pl.Float64, strict=False)
            )

            # Rename uniquely
            new_col = f"{target_col}__f{file_idx}_c{col_idx}"
            tmp = tmp.rename({target_col: new_col})

            dfs.append(tmp)
            detected_cols.append(new_col)

    # Safety check
    if not dfs:
        raise ValueError("No valid dataframes generated from input files")

    # -------------------------------------------------------
    # MERGE ALL
    # -------------------------------------------------------
    merged = dfs[0]
    for d in dfs[1:]:
        merged = merged.join(d, on=file_merge_keys, how="outer")

    # -------------------------------------------------------
    # AGGREGATION
    # -------------------------------------------------------
    exprs = [pl.col(c) for c in detected_cols]

    if method == "mean":
        agg = pl.mean_horizontal(exprs)
    elif method == "median":
        agg = pl.median_horizontal(exprs)
    elif method == "min":
        agg = pl.min_horizontal(exprs)
    elif method == "max":
        agg = pl.max_horizontal(exprs)
    else:
        raise ValueError(f"Invalid method '{method}'. Use mean/median/min/max")

    merged = merged.with_columns(agg.alias(final_col))

    # -------------------------------------------------------
    # CONTRIBUTION COUNT (VERSION-SAFE)
    # -------------------------------------------------------
    merged = merged.with_columns(
        pl.sum_horizontal([
            pl.when(pl.col(c).is_not_null() & pl.col(c).is_finite())
            .then(1)
            .otherwise(0)
            for c in detected_cols
        ]).alias(f"{final_col}_n")
    )

    # -------------------------------------------------------
    # RESOLVE EXTRA COLUMNS
    # -------------------------------------------------------
    for col in extra_columns:
        candidates = [c for c in merged.columns if c == col or col in c]
        if len(candidates) > 1:
            merged = merged.with_columns(
                pl.coalesce([pl.col(c) for c in candidates]).alias(col)
            )
            merged = merged.drop([c for c in candidates if c != col])

    # -------------------------------------------------------
    # FINAL SELECT
    # -------------------------------------------------------
    merged = merged.select(
        list(dict.fromkeys(
            file_merge_keys + [final_col, f"{final_col}_n"] + extra_columns
        ))
    )

    # Rename keys
    rename_map = dict(zip(file_merge_keys, sumstat_merge_keys))
    merged = merged.rename(rename_map)

    return merged, detected_cols





def process_external_files(
    files: List[str] | str,
    file_merge_keys: List[str],
    sumstat_merge_keys: List[str],
    column_suffix: str,
    method: str = "mean",
    final_col: str = "merged_value",
    extra_columns: List[str] | None = None,
    sep: str = "\t",
) -> Tuple[pl.DataFrame, List[str]]:

    if extra_columns is None:
        extra_columns = []

    if isinstance(files, str):
        files = [files]

    if len(file_merge_keys) != len(sumstat_merge_keys):
        raise ValueError("Merge keys mismatch")

    dfs = []
    detected_cols = []

    for file in files:
        df = pl.read_csv(file, separator=sep, infer_schema_length=10000)

        # Validate merge keys
        for k in file_merge_keys:
            if k not in df.columns:
                raise ValueError(f"Missing key '{k}' in {file}")

        # Detect column
        cols = [c for c in df.columns if c.lower().endswith(column_suffix.lower())]
        if not cols:
            raise ValueError(f"No column ending with '{column_suffix}' in {file}")

        target_col = cols[0]

        # Validate extra columns
        for col in extra_columns:
            if col not in df.columns:
                raise ValueError(f"Missing extra column '{col}' in {file}")

        select_cols = list(dict.fromkeys(
            file_merge_keys + [target_col] + extra_columns
        ))

        df = df.select(select_cols)

        dfs.append(df)
        detected_cols.append(target_col)

    # Merge
    merged = dfs[0]
    for d in dfs[1:]:
        merged = merged.join(d, on=file_merge_keys, how="outer", suffix="_dup")

    exprs = [pl.col(c) for c in detected_cols]

    if method == "mean":
        agg = pl.mean_horizontal(exprs)
    elif method == "median":
        agg = pl.median_horizontal(exprs)
    elif method == "min":
        agg = pl.min_horizontal(exprs)
    elif method == "max":
        agg = pl.max_horizontal(exprs)
    else:
        raise ValueError("Invalid method")

    merged = merged.with_columns(agg.alias(final_col))

    # Resolve duplicates
    for col in extra_columns:
        candidates = [c for c in merged.columns if c == col or c.startswith(f"{col}_dup")]
        if len(candidates) > 1:
            merged = merged.with_columns(
                pl.coalesce([pl.col(c) for c in candidates]).alias(col)
            )
            merged = merged.drop([c for c in candidates if c != col])

    merged = merged.select(
        list(dict.fromkeys(file_merge_keys + [final_col] + extra_columns))
    )

    rename_map = dict(zip(file_merge_keys, sumstat_merge_keys))
    merged = merged.rename(rename_map)

    return merged, files

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

    if len(file_merge_keys) != len(sumstat_merge_keys):
        raise ValueError("Merge keys mismatch")

    ✅ aggregation map (FIXED)
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
    detected_cols = []

    for file_idx, file in enumerate(files):
        dbg(f"\nProcessing file: {file}")

        df = pl.read_csv(file, separator=sep, infer_schema_length=10000)

        Validate keys
        for k in file_merge_keys:
            if k not in df.columns:
                raise ValueError(f"Missing key '{k}' in {file}")

        Detect suffix columns (robust)
        suffix = column_suffix.lower()
        cols = [
            c for c in df.columns
            if c.lower().endswith(f":{suffix}") or c.lower().endswith(suffix)
        ]

        if not cols:
            raise ValueError(f"No column ending with '{column_suffix}' in {file}")

        dbg(f"Detected columns: {cols}")

        for col_idx, target_col in enumerate(cols):
            tmp = df.select(
                list(dict.fromkeys(file_merge_keys + [target_col] + extra_columns))
            )

            tmp = tmp.with_columns(
                pl.col(target_col).cast(pl.Float64, strict=False)
            )

            new_col = f"{target_col}__f{file_idx}_c{col_idx}"
            tmp = tmp.rename({target_col: new_col})

            dfs.append(tmp)
            detected_cols.append(new_col)

    if not dfs:
        raise ValueError("No valid dataframes generated")

    Merge
    merged = dfs[0]
    for d in dfs[1:]:
        merged = merged.join(d, on=file_merge_keys, how="outer")

    dbg(f"Merged columns: {merged.columns}")

    Aggregate
    exprs = [pl.col(c).fill_null(0) for c in detected_cols]

    merged = merged.with_columns(
        agg_map[method](exprs).alias(final_col)
    )

    Contribution count
    merged = merged.with_columns(
        pl.sum_horizontal([
            pl.when(pl.col(c).is_not_null() & pl.col(c).is_finite())
            .then(1)
            .otherwise(0)
            for c in detected_cols
        ]).alias(f"{final_col}_n")
    )

    Resolve extra columns
    for col in extra_columns:
        candidates = [c for c in merged.columns if col in c]
        if len(candidates) > 1:
            merged = merged.with_columns(
                pl.coalesce([pl.col(c) for c in candidates]).alias(col)
            )
            merged = merged.drop([c for c in candidates if c != col])

    Final select
    merged = merged.select(
        list(dict.fromkeys(
            file_merge_keys + [final_col, f"{final_col}_n"] + extra_columns
        ))
    )

    Rename keys
    rename_map = dict(zip(file_merge_keys, sumstat_merge_keys))
    merged = merged.rename(rename_map)

    dbg(f"Final columns: {merged.columns}")

    return merged, detected_cols