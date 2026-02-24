###############################################################################
# CLI args → pipeline params
###############################################################################
args_to_params <- function(args) {

  message("▶ Converting CLI args → pipeline params")

  input <- normalizePath(args$input, mustWork = TRUE)
  placo_script <- normalizePath(args$placo_script, mustWork = TRUE)

  output <- normalizePath(args$output, mustWork = FALSE)

  params <- list(
    input = input,
    placo_script = placo_script,
    output_csv = output,
    method = args$method,
    id_col = args$id_col,
    z_suffix = args$z_suffix,
    p_suffix = args$p_suffix,
    p.threshold = args$pthreshold,
    run_time = Sys.time(),
    host = Sys.info()[["nodename"]],
    working_dir = getwd()
  )

  message("✅ Parameter object created")
  params
}

###############################################################################
# Pipeline wrapper
# ###############################################################################
# run_placo_pipeline <- function(params) {

#   message("▶ Running PLACO pipeline")

#   placo_pipeline(
#     input = params$input,
#     placo_script = params$placo_script,
#     output_csv = params$output_csv,
#     method = params$method,
#     id_col = params$id_col,
#     z_suffix = params$z_suffix,
#     p_suffix = params$p_suffix,
#     p.threshold = params$p.threshold
#   )

#   message("✅ PLACO pipeline completed")
# }


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