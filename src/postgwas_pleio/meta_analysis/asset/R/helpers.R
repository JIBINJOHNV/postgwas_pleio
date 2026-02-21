###############################################################################
# CLI args → pipeline params
###############################################################################
args_to_params <- function(args) {

  suppressPackageStartupMessages(library(glue))

  message("▶ Converting CLI args → pipeline params")

  # ---- Normalize ----
  input_manifest  <- normalizePath(args$input_manifest, mustWork = TRUE)
  fastasset_input <- normalizePath(args$fastasset_input, mustWork = TRUE)
  hm3             <- normalizePath(args$hm3, mustWork = TRUE)
  ld_ref          <- normalizePath(args$ld_ref, mustWork = TRUE)

  if (!dir.exists(ld_ref))
    stop("ld_ref must be a directory: ", ld_ref)

  output_dir <- normalizePath(args$output_dir, mustWork = FALSE)

  if (!dir.exists(output_dir))
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # ---- Derived dirs ----
  munge_dir    <- file.path(output_dir, "3_munge_output")
  ldsc_dir     <- file.path(output_dir, "4_ldsc_output")
  fastasset_dir<- file.path(output_dir, "5_fastasset_output")
  htraits_dir  <- file.path(output_dir, "6_htraits_output")

  dir.create(munge_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(ldsc_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fastasset_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(htraits_dir, recursive = TRUE, showWarnings = FALSE)

  params <- list(
    input_manifest   = input_manifest,
    fastasset_input  = fastasset_input,
    output_dir       = output_dir,
    hm3              = hm3,
    ld_ref           = ld_ref,
    run_name         = args$run_name,   # ⭐ FIX

    info_filter      = args$info_filter,
    maf_filter       = args$maf_filter,
    chunk_size       = args$chunk_size,
    scr_pthr         = args$scr_pthr,
    meth_pval        = args$meth_pval,
    ncores           = args$ncores,

    munge_dir        = munge_dir,
    ldsc_dir         = ldsc_dir,
    fastasset_dir    = fastasset_dir,
    htraits_dir      = htraits_dir,

    run_time         = Sys.time(),
    working_dir      = getwd(),
    host             = Sys.info()[["nodename"]]
  )

  message("✅ Parameter object created")
  params
}

###############################################################################
# Environment capture
###############################################################################
capture_environment <- function(output_dir) {

  message("▶ Capturing environment metadata")

  si <- sessionInfo()

  base_info <- data.frame(
    key = c("R.version","Platform","Running"),
    value = c(
      R.version.string,
      si$platform,
      paste(Sys.info()[["sysname"]], Sys.info()[["release"]])
    ),
    stringsAsFactors = FALSE
  )

  loaded_df <- if (length(si$otherPkgs) > 0) {
    do.call(rbind, lapply(names(si$otherPkgs), function(p) {
      data.frame(package=p, version=si$otherPkgs[[p]]$Version)
    }))
  } else NULL

  ns_df <- if (length(si$loadedOnly) > 0) {
    do.call(rbind, lapply(names(si$loadedOnly), function(p) {
      data.frame(package=p, version=si$loadedOnly[[p]]$Version)
    }))
  } else NULL

  pkg_df <- rbind(loaded_df, ns_df)

  fwrite(base_info, file.path(output_dir,"environment_base.tsv"), sep="\t")

  if (!is.null(pkg_df))
    fwrite(pkg_df, file.path(output_dir,"environment_packages.tsv"), sep="\t")

  writeLines(capture.output(sessionInfo()),
             file.path(output_dir,"sessionInfo.txt"))

  message("✅ Environment captured")
}

###############################################################################
# External tools
###############################################################################
capture_external_tools <- function(output_dir) {

  tools <- c("bcftools","plink","python","micromamba","gcloud")

  get_ver <- function(tool) {
    tryCatch(
      system2(tool,"--version",stdout=TRUE,stderr=TRUE)[1],
      error=function(e) NA
    )
  }

  df <- data.frame(tool=tools, version=sapply(tools,get_ver))
  df <- df[!is.na(df$version), ]

  fwrite(df,file.path(output_dir,"external_tools.tsv"))
}

###############################################################################
# safe_stage
###############################################################################
safe_stage <- function(stage_name, fn, report_env = NULL) {

  message("\n▶ Running stage: ", stage_name)

  start_time <- Sys.time()
  status <- "SUCCESS"

  result <- tryCatch(
    fn(),  # ✅ call the function (closure)
    error = function(e) {
      status <<- "FAILED"
      message("❌ ", stage_name, " failed")
      message("Reason: ", e$message)
      stop(e)
    }
  )

  end_time <- Sys.time()

  if (!is.null(report_env)) {
    report_env[[stage_name]] <- data.frame(
      stage = stage_name,
      status = status,
      start = start_time,
      end = end_time,
      runtime_sec = as.numeric(end_time - start_time, units = "secs"),
      stringsAsFactors = FALSE
    )
  }

  message("✅ ", stage_name, " completed in ",
          round(as.numeric(end_time - start_time, units = "secs"), 1), " sec")

  result
}

###############################################################################
# Report helpers
###############################################################################
init_pipeline_report <- function() {
  new.env(parent = emptyenv())
}

finalize_pipeline_report <- function(report_env) {

  if (length(ls(report_env)) == 0) return(NULL)

  report <- do.call(rbind, mget(ls(report_env), envir=report_env))
  rownames(report) <- NULL
  report
}