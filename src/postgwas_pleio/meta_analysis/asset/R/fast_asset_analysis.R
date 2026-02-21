############################################################
# 0) Libraries
############################################################
library(data.table)
library(glue)
library(parallel)
library(GenomicSEM)
library(ASSET)

options(expressions = 5e5)




  # ⭐ Load pipeline implementation FIRST
  # analysis_file_1 <- "/app/src/postgwas_pleio/meta_analysis/asset/R/fast_asset_analysis.R"
  # analysis_file_2 <- "/app/src/postgwas_pleio/meta_analysis/asset/R/helpers.R"

  # source(analysis_file_1)
  # source(analysis_file_2)


  


  run_munge_stage <- function(input_manifest,
                            hm3,
                            output_dir,
                            info_filter,
                            maf_filter,
                            run_name = "munge") {

  # ---- Read manifest ----
  df <- data.table::fread(input_manifest)

  # ---- Expected output files ----
  munged_expected <- file.path(output_dir, paste0(df$NAME, ".sumstats.gz"))

  if (all(file.exists(munged_expected))) {
    message("Munge exists → skipping")
    df$mung_file <- munged_expected
    return(df)
  }

  # ---- Ensure output directory ----
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # ---- Run munge with explicit paths ----
  munge(
    files       = df$FILE,
    hm3         = hm3,
    trait.names = file.path(output_dir, df$NAME),  # ⭐ KEY FIX
    log.name    = file.path(output_dir, run_name), # ⭐ LOG FIX
    info.filter = info_filter,
    maf.filter  = maf_filter,
    parallel    = TRUE
  )

  # ---- Record outputs ----
  df$mung_file <- munged_expected

  return(df)
}




run_ldsc_stage <- function(df, ld_ref, output_dir, run_name) {

  # ---- Ensure output dir ----
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  rdata <- file.path(output_dir, paste0(run_name, "_LDSCoutput.RData"))

  # ---- Resume logic ----
  if (file.exists(rdata)) {
    message("LDSC exists → loading")

    env <- new.env()
    load(rdata, envir = env)

    return(env$LDSCoutput)
  }

  # ---- Safety check ----
  if (!all(file.exists(df$mung_file))) {
    stop("Some munged files are missing:\n",
         paste(df$mung_file[!file.exists(df$mung_file)], collapse = "\n"))
  }

  # ---- Run LDSC ----
  LDSCoutput <- ldsc(
    traits          = df$mung_file,
    sample.prev     = df$SPREV,
    population.prev = df$PPREV,
    ld              = ld_ref,
    wld             = ld_ref,
    trait.names     = df$NAME
  )

  # ---- Fix dimension names (GenomicSEM quirk) ----
  dimnames(LDSCoutput$I)[[1]] <- dimnames(LDSCoutput$S)[[2]]
  dimnames(LDSCoutput$I)[[2]] <- dimnames(LDSCoutput$S)[[2]]

  # ---- Save ----
  save(LDSCoutput, file = rdata)

  return(LDSCoutput)
}





run_fastasset_stage <- function(asset_input_df,
                                LDSC_cor,
                                block,
                                output_dir,
                                chunk_size,
                                scr_pthr,
                                ncores,
                                run_name) {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # ---- Detect columns once ----
  beta_cols <- grep("\\.Beta$", names(asset_input_df), value = TRUE)
  se_cols   <- grep("\\.SE$",   names(asset_input_df), value = TRUE)
  n_cols    <- grep("\\.N$",    names(asset_input_df), value = TRUE)

  if (length(beta_cols) == 0 || length(se_cols) == 0 || length(n_cols) == 0)
    stop("Cannot detect Beta/SE/N columns")

  traits <- sub("\\.Beta$", "", beta_cols)

  # ---- Chunk index ----
  idx <- split(seq_len(nrow(asset_input_df)),
               ceiling(seq_len(nrow(asset_input_df)) / chunk_size))

  message("▶ fastASSET stage started")
  message("   Total SNPs: ", nrow(asset_input_df))
  message("   Total chunks: ", length(idx))

  for (c in seq_along(idx)) {

    message("\n▶ Processing chunk ", c, "/", length(idx))

    out1 <- file.path(output_dir, paste0("fastasset_chunk_", c, "_1sided.tsv.gz"))
    out2 <- file.path(output_dir, paste0("fastasset_chunk_", c, "_2sided.tsv.gz"))

    if (file.exists(out1) && file.exists(out2)) {
      message("   Chunk already completed → skipping")
      next
    }

    rows <- idx[[c]]

    chunk_res <- mclapply(rows, function(i) {

      snp_id <- asset_input_df$ID[i]

      # ⭐ STREAMING extraction (no duplication)
      beta_vec <- as.numeric(unlist(asset_input_df[i, ..beta_cols]))
      se_vec   <- as.numeric(unlist(asset_input_df[i, ..se_cols]))
      n_vec    <- as.numeric(unlist(asset_input_df[i, ..n_cols]))

      if (anyNA(beta_vec) || anyNA(se_vec) || anyNA(n_vec) ||
          any(!is.finite(beta_vec)) || any(!is.finite(se_vec)) || any(!is.finite(n_vec)) ||
          any(se_vec <= 0) || any(n_vec <= 0)) {
        return(NULL)
      }

      res_dt <- tryCatch({

        test <- fast_asset(
          snp        = snp_id,
          traits.lab = traits,
          beta.hat   = beta_vec,
          sigma.hat  = se_vec,
          Neff       = n_vec,
          cor        = LDSC_cor,
          block      = block,
          scr_pthr   = scr_pthr
        )

        hs <- h.summary(test)

        out <- list()

        if (!is.null(hs$Subset.1sided)) {
          dt <- as.data.table(hs$Subset.1sided)
          dt[, ID := snp_id]
          dt[, analysis_type := "1sided"]
          out[[length(out)+1]] <- dt
        }

        if (!is.null(hs$Subset.2sided)) {
          dt <- as.data.table(hs$Subset.2sided)
          dt[, ID := snp_id]
          dt[, analysis_type := "2sided"]
          out[[length(out)+1]] <- dt
        }

        if (length(out) == 0) return(NULL)
        rbindlist(out, fill = TRUE)

      }, error = function(e) NULL)

      if (is.null(res_dt) || !is.data.frame(res_dt)) return(NULL)

      res_dt

    }, mc.cores = ncores)

    chunk_res <- Filter(is.data.frame, chunk_res)

    if (length(chunk_res) == 0) next

    chunk_df <- rbindlist(chunk_res, fill = TRUE)

    fwrite(chunk_df[analysis_type=="1sided"], out1)
    fwrite(chunk_df[analysis_type=="2sided"], out2)

    message("   Chunk ", c, " completed")
  }

  message("\n✅ fastASSET stage finished")
}




run_htraits_stage <- function(asset_input_df,
                              LDSC_cor,
                              trait_names,
                              output_dir,
                              chunk_size = 100000,
                              meth_pval = "DLM",
                              ncores = 5,
                              run_name) {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  message("▶ h.traits parallel chunk stage started")

  # ---- Detect columns ----
  beta_cols <- paste0(trait_names, ".Beta")
  se_cols   <- paste0(trait_names, ".SE")
  n_cols    <- paste0(trait_names, ".N")

  if (!all(beta_cols %in% names(asset_input_df)))
    stop("Missing Beta columns")
  if (!all(se_cols %in% names(asset_input_df)))
    stop("Missing SE columns")
  if (!all(n_cols %in% names(asset_input_df)))
    stop("Missing N columns")

  # ---- Chunk index ----
  idx <- split(
    seq_len(nrow(asset_input_df)),
    ceiling(seq_len(nrow(asset_input_df)) / chunk_size)
  )

  message("Total chunks: ", length(idx))

  # ---- Worker ----
  run_chunk <- function(chunk_id) {

    rows <- idx[[chunk_id]]

    outfile <- file.path(output_dir, paste0("htraits_chunk_", chunk_id, ".tsv.gz"))
    if (file.exists(outfile)) return(TRUE)

    # ⭐ MEMORY-SAFE extraction
    chunk_beta <- as.matrix(asset_input_df[rows, ..beta_cols])
    chunk_se   <- as.matrix(asset_input_df[rows, ..se_cols])
    chunk_n    <- as.matrix(asset_input_df[rows, ..n_cols])

    storage.mode(chunk_beta) <- "double"
    storage.mode(chunk_se)   <- "double"
    storage.mode(chunk_n)    <- "double"

    snp_ids <- asset_input_df$ID[rows]

    res <- tryCatch({

      h.traits(
        snp.vars   = snp_ids,
        traits.lab = trait_names,
        beta.hat   = chunk_beta,
        sigma.hat  = chunk_se,
        ncase      = chunk_n,
        ncntl      = chunk_n,
        cor        = LDSC_cor,
        meta       = TRUE,
        meth.pval  = meth_pval
      )

    }, error = function(e) NULL)

    if (is.null(res)) return(FALSE)

    out <- list()

    if (!is.null(res$Meta)) {
      dt <- as.data.table(res$Meta)
      dt[, ID := snp_ids]
      dt[, analysis_type := "meta"]
      out[[length(out)+1]] <- dt
    }

    if (!is.null(res$Subset.1sided)) {
      dt <- as.data.table(res$Subset.1sided)
      dt[, ID := snp_ids]
      dt[, analysis_type := "1sided"]
      out[[length(out)+1]] <- dt
    }

    if (!is.null(res$Subset.2sided)) {
      dt <- as.data.table(res$Subset.2sided)
      dt[, ID := snp_ids]
      dt[, analysis_type := "2sided"]
      out[[length(out)+1]] <- dt
    }

    if (length(out) > 0)
      fwrite(rbindlist(out, fill = TRUE), outfile)

    TRUE
  }

  mclapply(seq_along(idx), run_chunk, mc.cores = ncores)

  message("✅ h.traits parallel chunk stage finished")

  # ---- Combine ----
  chunk_files <- sort(list.files(output_dir,
                                 pattern="^htraits_chunk_.*\\.tsv\\.gz$",
                                 full.names=TRUE))

  if (length(chunk_files) == 0) return(invisible(NULL))

  meta_file <- file.path(output_dir, paste0(run_name, "_htraits_meta.tsv.gz"))
  one_file  <- file.path(output_dir, paste0(run_name, "_htraits_1sided.tsv.gz"))
  two_file  <- file.path(output_dir, paste0(run_name, "_htraits_2sided.tsv.gz"))

  first_meta <- first_one <- first_two <- TRUE

  for (f in chunk_files) {

    dt <- fread(f)
    if (nrow(dt) == 0) next

    if (any(dt$analysis_type=="meta")) {
      fwrite(dt[analysis_type=="meta"], meta_file, append=!first_meta)
      first_meta <- FALSE
    }

    if (any(dt$analysis_type=="1sided")) {
      fwrite(dt[analysis_type=="1sided"], one_file, append=!first_one)
      first_one <- FALSE
    }

    if (any(dt$analysis_type=="2sided")) {
      fwrite(dt[analysis_type=="2sided"], two_file, append=!first_two)
      first_two <- FALSE
    }
  }

  message("✅ Combined files created")
}





run_asset_pipeline <- function(params) {

  suppressPackageStartupMessages(library(data.table))

  # ------------------------------------------------
  # 0) Directory structure
  # ------------------------------------------------
  munge_dir     <- if (!is.null(params$munge_dir))     params$munge_dir     else file.path(params$output_dir, "3_munge_output")
  ldsc_dir      <- if (!is.null(params$ldsc_dir))      params$ldsc_dir      else file.path(params$output_dir, "4_ldsc_output")
  fastasset_dir <- if (!is.null(params$fastasset_dir)) params$fastasset_dir else file.path(params$output_dir, "5_fastasset_output")
  htraits_dir   <- if (!is.null(params$htraits_dir))   params$htraits_dir   else file.path(params$output_dir, "6_htraits_output")

  dir.create(munge_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(ldsc_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fastasset_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(htraits_dir, recursive = TRUE, showWarnings = FALSE)

  # ------------------------------------------------
  # 1) Report env
  # ------------------------------------------------
  report_env <- init_pipeline_report()

  # ------------------------------------------------
  # 2) Munge
  # ------------------------------------------------
  df <- safe_stage(
    "munge",
    function() run_munge_stage(
      input_manifest = params$input_manifest,
      hm3            = params$hm3,
      output_dir     = munge_dir,
      info_filter    = params$info_filter,
      maf_filter     = params$maf_filter,
      run_name       = params$run_name
    ),
    report_env
  )

  # ------------------------------------------------
  # 3) LDSC
  # ------------------------------------------------
  LDSCoutput <- safe_stage(
    "ldsc",
    function() run_ldsc_stage(
      df         = df,
      ld_ref     = params$ld_ref,
      output_dir = ldsc_dir,
      run_name   = params$run_name
    ),
    report_env
  )

  # ------------------------------------------------
  # 4) Load fastASSET
  # ------------------------------------------------
  asset_input_df <- safe_stage(
    "load_fastasset_input",
    function() data.table::fread(params$fastasset_input),
    report_env
  )

  # ------------------------------------------------
  # 5) Trait list
  # ------------------------------------------------
  beta_cols <- grep("\\.Beta$", names(asset_input_df), value = TRUE)
  if (length(beta_cols) == 0) stop("No .Beta columns found")

  traits <- sub("\\.Beta$", "", beta_cols)

  if (!all(traits %in% colnames(LDSCoutput$I))) {
    missing <- setdiff(traits, colnames(LDSCoutput$I))
    stop("Traits missing from LDSCoutput$I: ", paste(missing, collapse = ", "))
  }

  LDSC_cor <- LDSCoutput$I[traits, traits]
  block <- create_blocks(LDSC_cor)

  # ------------------------------------------------
  # 6) fastASSET
  # ------------------------------------------------
  safe_stage(
    "fastasset",
    function() run_fastasset_stage(
      asset_input_df = asset_input_df,
      LDSC_cor       = LDSC_cor,
      block          = block,
      output_dir     = fastasset_dir,
      chunk_size     = params$chunk_size,
      scr_pthr       = params$scr_pthr,
      ncores         = params$ncores,
      run_name       = params$run_name
    ),
    report_env
  )

  # ------------------------------------------------
  # 7) h.traits
  # ------------------------------------------------
  safe_stage(
    "htraits",
    function() run_htraits_stage(
      asset_input_df = asset_input_df,
      LDSC_cor       = LDSC_cor,
      trait_names    = traits,
      output_dir     = htraits_dir,
      chunk_size     = params$chunk_size,
      meth_pval      = params$meth_pval,
      ncores         = params$ncores,
      run_name       = params$run_name
    ),
    report_env
  )

  # ------------------------------------------------
  # 8) metadata
  # ------------------------------------------------
  safe_stage("environment_capture", function() capture_environment(params$output_dir), report_env)
  safe_stage("external_tools_capture", function() capture_external_tools(params$output_dir), report_env)

  message("✅ Pipeline completed successfully")

  invisible(finalize_pipeline_report(report_env))
}



