###############################################################################
# ASSET CLI
# Minimal, stack-safe, production CLI
###############################################################################

options(expressions = 5e5)

###############################################################################
# 0) Libraries (ONLY CLI-level)
###############################################################################
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(glue))

# If this is packaged, consider making script_dir relative to this file.
# For now you hard-coded it:
script_dir <- "/app/src/postgwas_pleio/meta_analysis/asset/R/"

###############################################################################
# 1) Parser
###############################################################################
create_asset_cli <- function() {

  parser <- ArgumentParser(
    prog = "asset_pipeline",
    formatter_class = "argparse.RawTextHelpFormatter",
    description =
"ASSET MULTI-TRAIT PLEIOTROPY PIPELINE
=====================================

Pipeline stages:
  1) GenomicSEM::munge
  2) LDSC
  3) fastASSET
  4) ASSET h.traits

INPUTS:
  --input_manifest
  --fastasset_input
  --hm3
  --ld_ref

OUTPUT:
  3_munge_output/
  4_ldsc_output/
  5_fastasset_output/
  6_htraits_output/
"
  )

  # ---------- REQUIRED ----------
  parser$add_argument("--input_manifest", required = TRUE)
  parser$add_argument("--fastasset_input", required = TRUE)
  parser$add_argument("--output_dir", required = TRUE)
  parser$add_argument("--run_name", required = TRUE)
  parser$add_argument("--hm3", required = TRUE)
  parser$add_argument("--ld_ref", required = TRUE)

  # ---------- OPTIONAL ----------
  parser$add_argument("--info_filter", type = "double", default = 0.9)
  parser$add_argument("--maf_filter",  type = "double", default = 0.01)
  parser$add_argument("--chunk_size",  type = "integer", default = 100000)
  parser$add_argument("--scr_pthr",    type = "double", default = 0.999999)
  parser$add_argument("--meth_pval",   default = "DLM")

  parser$add_argument(
    "--ncores",
    type = "integer",
    default = max(1, floor(parallel::detectCores() * 0.5))
  )

  parser
}

###############################################################################
# 2) Validation
###############################################################################
validate_inputs <- function(args) {

  message("▶ Validating inputs")

  for (f in c(args$input_manifest, args$fastasset_input, args$hm3, args$ld_ref)) {
    # ld_ref can be a directory; handle that:
    if (f == args$ld_ref) {
      if (!dir.exists(f)) stop("Directory not found: ", f)
    } else {
      if (!file.exists(f)) stop("File not found: ", f)
    }
  }

  # ---- manifest header only ----
  man <- data.table::fread(args$input_manifest, nrows = 0)
  req_cols <- c("FILE", "NAME", "SPREV", "PPREV")
  if (!all(req_cols %in% names(man))) {
    stop("Manifest missing columns: ",
         paste(setdiff(req_cols, names(man)), collapse = ", "))
  }

  # ---- fastasset header only ----
  fa <- data.table::fread(args$fastasset_input, nrows = 0)
  if (!"ID" %in% names(fa)) stop("fastASSET input missing ID")
  if (!any(grepl("\\.Beta$", names(fa)))) stop("No .Beta columns")
  if (!any(grepl("\\.SE$",   names(fa)))) stop("No .SE columns")
  if (!any(grepl("\\.N$",    names(fa)))) stop("No .N columns")

  message("✅ Input validation passed")
}

###############################################################################
# 3) Main runner
###############################################################################
run_cli <- function() {

  parser <- create_asset_cli()

  if (length(commandArgs(trailingOnly = TRUE)) == 0) {
    parser$print_help()
    quit(status = 0)
  }

  # ---- Load pipeline layer ONLY here ----
  pipeline_file <- file.path(script_dir, "fast_asset_analysis.R")
  helpers_file  <- file.path(script_dir, "helpers.R")

  if (!file.exists(pipeline_file))
    stop("Missing fast_asset_analysis.R in: ", script_dir)
  if (!file.exists(helpers_file))
    stop("Missing helpers.R in: ", script_dir)

  # local=TRUE keeps symbols scoped and avoids polluting global env
  source(helpers_file)
  source(pipeline_file)

  # ---- Parse args ----
  args <- parser$parse_args()

  # ---- Validate ----
  validate_inputs(args)

  # ---- Convert ----
  params <- args_to_params(args)

  # ---- Run ----
  run_asset_pipeline(params)
}

###############################################################################
# 4) Entry
###############################################################################
if (!interactive()) {
  run_cli()
}