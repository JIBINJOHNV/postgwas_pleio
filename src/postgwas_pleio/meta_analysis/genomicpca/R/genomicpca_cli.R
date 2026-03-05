###############################################################################
# GENOMIC PCA CLI
# Minimal, stack-safe, production CLI
###############################################################################

options(expressions = 5e5)

###############################################################################
# 0) Libraries (ONLY CLI-level)
###############################################################################
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(glue))

script_dir <- "/app/src/postgwas_pleio/meta_analysis/genomicpca/R/"

###############################################################################
# 1) Parser
###############################################################################
create_genomicpca_cli <- function() {

  parser <- ArgumentParser(
    prog = "genomicpca_pipeline",
    formatter_class = "argparse.RawTextHelpFormatter",
    description =
"GENOMIC PCA MULTI-TRAIT GWAMA PIPELINE
=======================================

Pipeline stages:
  1) Munge (GenomicSEM)
  2) LDSC
  3) GenomicPCA
  4) Multivariate GWAMA

INPUTS:
  --input_manifest
  --hm3
  --ld_ref

OUTPUT:
  3_ldsc_analysis/
  4_GenomicPCA_correlation/
  5_GenomicPCA_alternative_covariates/
"
  )

  # ---------- REQUIRED ----------
  parser$add_argument("--input_manifest", required = TRUE)
  parser$add_argument("--output_dir", required = TRUE)
  parser$add_argument("--run_name", required = TRUE)
  parser$add_argument("--hm3", required = TRUE)
  parser$add_argument("--ld_ref", required = TRUE)

  # ---------- OPTIONAL ----------
  parser$add_argument("--info_filter", type = "double", default = 0.7)
  parser$add_argument("--maf_filter",  type = "double", default = 0.3)

  parser$add_argument(
    "--approach",
    type = "character",
    default = "both",
    choices = c("correlation","covariance","both"),
    help = "GenomicPCA approach"
  )

  parser
}

###############################################################################
# 2) Validation
###############################################################################
validate_inputs <- function(args) {

  message("▶ Validating inputs")

  # ---- check files ----
  for (f in c(args$input_manifest, args$hm3)) {
    if (!file.exists(f))
      stop("File not found: ", f)
  }

  if (!dir.exists(args$ld_ref))
    stop("Directory not found: ", args$ld_ref)

  # ---- manifest header only ----
  man <- data.table::fread(args$input_manifest, nrows = 0)

  req_cols <- c("FILE","NAME","SPREV","PPREV")

  if (!all(req_cols %in% names(man))) {
    stop(
      "Manifest missing columns: ",
      paste(setdiff(req_cols, names(man)), collapse = ", ")
    )
  }

}

###############################################################################
# 3) Main runner
###############################################################################
run_cli <- function() {

  parser <- create_genomicpca_cli()

  if (length(commandArgs(trailingOnly = TRUE)) == 0) {
    parser$print_help()
    quit(status = 0)
  }

  # ---- Load pipeline layer ONLY here ----
  pipeline_file <- file.path(script_dir, "genomicpca_analysis.R")
  helpers_file  <- file.path(script_dir, "helpers.R")
  gwama_file    <- file.path(script_dir, "N_weighted_GWAMA.function.1_2_6.R")

  if (!file.exists(pipeline_file))
    stop("Missing genomicpca_analysis.R in: ", script_dir)

  if (!file.exists(helpers_file))
    stop("Missing helpers.R in: ", script_dir)

  if (!file.exists(gwama_file))
    stop("Missing my_GWAMA_26032020.R in: ", script_dir)

  source(helpers_file)
  source(gwama_file)
  source(pipeline_file)

  # ---- Parse args ----
  args <- parser$parse_args()

  # ---- Validate ----
  validate_inputs(args)

  # ---- Convert to params ----
  params <- args_to_params(args)

  # ---- Run pipeline ----
  run_genomicpca_pipeline(params)

}

###############################################################################
# 4) Entry
###############################################################################
if (!interactive()) {
  run_cli()
}