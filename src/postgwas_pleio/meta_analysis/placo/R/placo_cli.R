###############################################################################
# PLACO CLI
# Minimal, stack-safe, production CLI
###############################################################################

options(expressions = 5e5)

###############################################################################
# 0) Libraries (CLI-level ONLY)
###############################################################################
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(glue))

# Adjust this to your module location
script_dir <- "/app/src/postgwas_pleio/meta_analysis/placo/R/"

###############################################################################
# 1) Parser
###############################################################################
create_placo_cli <- function() {

  parser <- ArgumentParser(
    prog = "placo_pipeline",
    formatter_class = "argparse.RawTextHelpFormatter",
    description =
"PLACO / PLACO+ PLEIOTROPY PIPELINE
===================================

Variant-level pleiotropy test for two traits.

INPUT:
  Harmonised master TSV containing:
    ID
    TRAIT1_Z
    TRAIT2_Z
    TRAIT1_P
    TRAIT2_P

OUTPUT:
  CSV with PLACO statistics and p-values
"
  )

  # ---------- REQUIRED ----------
  parser$add_argument("--input", required = TRUE)
  parser$add_argument("--placo_script", required = TRUE)
  parser$add_argument("--output", required = TRUE)

  # ---------- OPTIONAL ----------
  parser$add_argument("--method", default = "placo.plus",
                      choices = c("placo","placo.plus"))

  parser$add_argument("--id_col", default = "ID")
  parser$add_argument("--z_suffix", default = "_Z")
  parser$add_argument("--p_suffix", default = "_P")

  parser$add_argument("--pthreshold", type = "double", default = 1e-4)

  parser
}

###############################################################################
# 2) Validation
###############################################################################
validate_inputs <- function(args) {

  message("▶ Validating inputs")

  if(!file.exists(args$input))
    stop("Input file not found: ", args$input)

  if(!file.exists(args$placo_script))
    stop("PLACO script not found: ", args$placo_script)

  # ---- header check ----
  hdr <- data.table::fread(args$input, nrows = 0)

  if(!args$id_col %in% names(hdr))
    stop("ID column missing")

  if(!any(grepl(paste0(args$z_suffix,"$"), names(hdr))))
    stop("No *_Z columns")

  if(!any(grepl(paste0(args$p_suffix,"$"), names(hdr))))
    stop("No *_P columns")

  message("✅ Input validation passed")
}

###############################################################################
# 3) Main runner
###############################################################################
run_cli <- function() {

  parser <- create_placo_cli()

  if (length(commandArgs(trailingOnly = TRUE)) == 0) {
    parser$print_help()
    quit(status = 0)
  }

  # ---- Load pipeline layer ONLY here ----
  pipeline_file <- file.path(script_dir, "placo_pipeline.R")
  helpers_file  <- file.path(script_dir, "helpers.R")

  if (!file.exists(pipeline_file))
    stop("Missing placo_pipeline.R in: ", script_dir)

  if (!file.exists(helpers_file))
    stop("Missing helpers.R in: ", script_dir)

  source(helpers_file, local = TRUE)
  source(pipeline_file, local = TRUE)

  # ---- Parse args ----
  args <- parser$parse_args()

  # ---- Validate ----
  validate_inputs(args)

  # ---- Convert ----
  params <- args_to_params(args)

  # ---- Run ----

  run_placo_pipeline(
    input = params$input,
    placo_script = params$placo_script,
    output_csv = params$output,
    method = params$method,
    id_col = params$id_col,
    z_suffix = params$z_suffix,
    p_suffix = params$p_suffix,
    p.threshold = params$p.threshold
  )
}



###############################################################################
# 4) Entry
###############################################################################
if (!interactive()) {
  run_cli()
}