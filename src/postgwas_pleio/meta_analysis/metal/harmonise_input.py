

import polars as pl
from pathlib import Path


def metal_sumstat_to_harmonise_input(
    metal_file,
    info_file=None,
    info_method="mean",
    sample_size_approach="totalnef",
    out_tsv=None,
    manifest_out_tsv=None,
    gwas_outputname=None,
    resource_folder=None,
    output_folder=None,
):

    # -----------------------------
    # Allowed options
    # -----------------------------
    allowed_info_methods = ["mean", "median", "min", "max"]
    allowed_sample_sizes = ["weight", "sample_overlap_corrected", "totalnef"]

    info_method = info_method.lower()
    sample_size_approach = sample_size_approach.lower()

    if info_method not in allowed_info_methods:
        raise ValueError(f"info_method must be one of {allowed_info_methods}")

    if sample_size_approach not in allowed_sample_sizes:
        raise ValueError(
            f"sample_size_approach must be one of {allowed_sample_sizes}"
        )

    # -----------------------------
    # Map CLI option → METAL column
    # -----------------------------
    sample_map = {
        "weight": "Weight",
        "sample_overlap_corrected": "N",
        "totalnef": "TotalNEF",
    }

    sample_col = sample_map[sample_size_approach]

    # -----------------------------
    # Load METAL output
    # -----------------------------
    df = pl.read_csv(metal_file, separator="\t")

    if sample_col not in df.columns:
        raise ValueError(f"{sample_col} not found in METAL output")


    z_col = "Zscore" if "Zscore" in df.columns else None
    beta_col = "Effect" if "Effect" in df.columns else None
    se_col = "StdErr" if "StdErr" in df.columns else None

    cols = [
        "MarkerName",
        "Allele1",
        "Allele2",
        "Freq1",
        z_col,
        "P-value",
        sample_col,
        beta_col,
        se_col,
    ] 

    df = df.select([c for c in cols if c in df.columns and c is not None])

    df = df.rename({sample_col: "sample_size"})

    # -----------------------------
    # INFO score aggregation
    # -----------------------------
    if info_file:

        info_df = pl.read_csv(info_file, separator="\t")
        info_cols = info_df.columns[1:]

        if info_method == "mean":
            info_expr = pl.mean_horizontal(pl.col(info_cols))

        elif info_method == "min":
            info_expr = pl.min_horizontal(pl.col(info_cols))

        elif info_method == "max":
            info_expr = pl.max_horizontal(pl.col(info_cols))

        elif info_method == "median":
            info_expr = pl.concat_list(pl.col(info_cols)).list.median()

        info_df = (
            info_df.with_columns(info_expr.alias("info_score"))
            .select([info_df.columns[0], "info_score"])
            .rename({info_df.columns[0]: "MarkerName"})
        )

        df = df.join(info_df, on="MarkerName", how="left")

    # -----------------------------
    # Extract CHR / POS
    # -----------------------------
    parts = pl.col("MarkerName").str.split("_")

    df = df.with_columns([
        parts.list.get(0).alias("CHR"),
        parts.list.get(1).cast(pl.Int64).alias("POS"),
    ])

    # -----------------------------
    # Chromosome sorting
    # -----------------------------
    df = (
        df.with_columns(
            pl.when(pl.col("CHR") == "X")
            .then(23)
            .when(pl.col("CHR") == "Y")
            .then(24)
            .otherwise(pl.col("CHR").cast(pl.Int64))
            .alias("_chr_sort")
        )
        .sort(["_chr_sort", "POS"])
        .drop("_chr_sort")
    )

    # -----------------------------
    # Uppercase alleles
    # -----------------------------
    df = df.with_columns([
        pl.col("Allele1").str.to_uppercase(),
        pl.col("Allele2").str.to_uppercase(),
    ])

    # -----------------------------
    # Save harmonised summary stats
    # -----------------------------
    if out_tsv:
        out_path = Path(out_tsv)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        df.write_csv(out_path, separator="\t")

    # -----------------------------
    # Create manifest
    # -----------------------------
    manifest_df = None

    if out_tsv and gwas_outputname:

        manifest_df = pl.DataFrame({
            "sumstat_file": [str(out_tsv)],
            "gwas_outputname": [gwas_outputname],
            "chr_col": ["CHR"],
            "pos_col": ["POS"],
            "snp_id_col": ["MarkerName"],
            "ea_col": ["Allele1"],
            "oa_col": ["Allele2"],
            "eaf_col": ["Freq1"],
            "beta_or_col": [beta_col],
            "se_col": [se_col],
            "imp_z_col": [z_col],
            "pval_col": ["P-value"],
            "ncontrol_col": ["sample_size"],
            "ncase_col": ["NA"],
            "ncontrol": ["NA"],
            "ncase": ["NA"],
            "imp_info_col": ["info_score"] if "info_score" in df.columns else ["NA"],
            "infofile": ["NA"],
            "infocolumn": ["NA"],
            "eaffile": ["NA"],
            "eafcolumn": ["NA"],
            "liftover": ["Yes"],
            "chr_pos_col": ["NA"],
            "resourse_folder": [resource_folder],
            "output_folder": [output_folder],
        })

        if manifest_out_tsv:
            manifest_path = Path(manifest_out_tsv)
            manifest_path.parent.mkdir(parents=True, exist_ok=True)
            manifest_df.write_csv(manifest_path, separator="\t")

    return df, manifest_df