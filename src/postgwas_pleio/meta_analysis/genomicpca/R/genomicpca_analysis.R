############################################################
# 0) Libraries
############################################################
library(data.table)
library(glue)
library(parallel)
library(GenomicSEM)
library(ASSET)
library(checkmate)

options(expressions = 5e5)


  
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
      files       = df$ldsc_FILE,
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
        trait.names     = df$NAME,
        stand = T
    )

    # ---- Fix dimension names (GenomicSEM quirk) ----
    dimnames(LDSCoutput$I)[[1]] <- dimnames(LDSCoutput$S)[[2]]
    dimnames(LDSCoutput$I)[[2]] <- dimnames(LDSCoutput$S)[[2]]

    # ---- Save ----
    save(LDSCoutput, file = rdata)
    return(LDSCoutput)
}






run_genomicpca_gwama <- function(
    input_file_df,
    ldsc,
    output_dir,
    approach = c("both","correlation","covariance")
){
    ## Reference: https://annafurtjes.github.io/genomicPCA/25082021_geneticPCA_explanation.html
    library(data.table)
    library(glue)
    approach <- match.arg(approach)

    # ------------------------------------------------------------
    # 1. Load GWAS summary statistics listed in manifest
    # ------------------------------------------------------------
    dat <- vector("list")

    for(i in 1:nrow(input_file_df)){
        name <- input_file_df$NAME[i]
        file <- input_file_df$FILE[i]
        print(name)
        dat[[name]] <- fread(file, data.table = FALSE)
    }

    # ------------------------------------------------------------
    # Ensure LDSC intercept matrix matches GWAS ordering
    # ------------------------------------------------------------
    order <- names(dat)

    dimnames(ldsc$I)[[1]] <- dimnames(ldsc$S)[[2]]
    dimnames(ldsc$I)[[2]] <- dimnames(ldsc$S)[[2]]
    CTI <- as.matrix(ldsc$I[order, order])




    # ============================================================
    # CORRELATION APPROACH (recommended)
    # ============================================================
    if(approach %in% c("both","correlation")){

        dimnames(ldsc$S_Stand)[[1]] <- dimnames(ldsc$S)[[2]]
        dimnames(ldsc$S_Stand)[[2]] <- dimnames(ldsc$S)[[2]]
        cormatrix <- ldsc$S_Stand[order, order]

        # -------------------------------
        # Run GenomicPCA diagnostics
        # -------------------------------
        run_genomicpca_diagnostics(
            matrix_input = cormatrix,
            trait_names = order,
            output_prefix = paste0(output_dir,"/4_GenomicPCA_correlation/GenomicPCA_correlation")
        )

        # PCA decomposition
        eigenvectors <- eigen(cormatrix)$vectors
        eigenvalues  <- eigen(cormatrix)$values

        # Standardized PCA loadings
        loadings <- as.vector(eigenvectors %*% sqrt(diag(eigenvalues))[,1])

        # Create output directory
        gpca_cordir <- paste0(output_dir,"/4_GenomicPCA_correlation")
        dir.create(gpca_cordir, recursive = TRUE, showWarnings = FALSE)
        setwd(gpca_cordir)

        # Run GWAMA
        multivariate_GWAMA(
            x = dat,
            cov_Z = CTI,
            h2 = loadings,
            out = gpca_cordir,
            name = "GenomicPCA_correlation_meta",
            output_gz = FALSE,
            check_columns = FALSE
        )
    }


    # ============================================================
    # COVARIANCE APPROACH (alternative)
    # ============================================================
    if(approach %in% c("both","covariance")){

        dimnames(ldsc$S)[[1]] <- dimnames(ldsc$S)[[2]]
        covmatrix <- ldsc$S[order, order]

        # -------------------------------
        # Run GenomicPCA diagnostics
        # -------------------------------
        run_genomicpca_diagnostics(
            matrix_input = covmatrix,
            trait_names = order,
            output_prefix = paste0(output_dir,"/5_GenomicPCA_alternative_covariates/GenomicPCA_covariance")
        )
        
        # PCA decomposition
        eigenvectors <- eigen(covmatrix)$vectors
        eigenvalues  <- eigen(covmatrix)$values

        # Standardized PCA loadings
        loadings <- as.vector(eigenvectors %*% sqrt(diag(eigenvalues))[,1])

        # Flip direction if median negative
        median_loadings <- median(loadings)
        mean_direction  <- sign(median_loadings)

        if(mean_direction == -1){
            loadings <- loadings * (-1)
        }

        # Create output directory
        gpca_covdir <- paste0(output_dir,"/5_GenomicPCA_alternative_covariates")
        dir.create(gpca_covdir, recursive = TRUE, showWarnings = FALSE)
        setwd(gpca_covdir)

        # Run GWAMA
        multivariate_GWAMA(
            x = dat,
            cov_Z = CTI,
            h2 = loadings,
            out = gpca_covdir,
            name = "GenomicPCA_covariates_meta",
            output_gz = FALSE,
            check_columns = FALSE
        )
    }

}

run_genomicpca_pipeline <- function(params) {

  suppressPackageStartupMessages(library(data.table))

  # ------------------------------------------------
  # 0) Directory structure
  # ------------------------------------------------
  ldsc_dir <- if (!is.null(params$ldsc_dir)) {
    params$ldsc_dir
  } else {
    file.path(params$output_dir, "3_ldsc_analysis")
  }

  gpca_cor_dir <- if (!is.null(params$gpca_cor_dir)) {
    params$gpca_cor_dir
  } else {
    file.path(params$output_dir, "4_GenomicPCA_correlation")
  }

  gpca_cov_dir <- if (!is.null(params$gpca_cov_dir)) {
    params$gpca_cov_dir
  } else {
    file.path(params$output_dir, "5_GenomicPCA_alternative_covariates")
  }

  ldsc_input_dir  <- file.path(ldsc_dir, "ldsc_input")
  ldsc_output_dir <- file.path(ldsc_dir, "ldsc_output")

  dir.create(ldsc_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(ldsc_input_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(ldsc_output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(gpca_cor_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(gpca_cov_dir, recursive = TRUE, showWarnings = FALSE)

  # ------------------------------------------------
  # 1) Report environment
  # ------------------------------------------------
  report_env <- init_pipeline_report()

  # ------------------------------------------------
  # 2) Munge stage
  # ------------------------------------------------
  munge_df <- safe_stage(
    "munge",
    function() run_munge_stage(
      input_manifest = params$input_manifest,
      hm3            = params$hm3,
      output_dir     = ldsc_input_dir,
      info_filter    = params$info_filter,
      maf_filter     = params$maf_filter,
      run_name       = params$run_name
    ),
    report_env
  )

  # ------------------------------------------------
  # 3) LDSC stage
  # ------------------------------------------------
  ldsc <- safe_stage(
    "ldsc",
    function() run_ldsc_stage(
      df         = munge_df,
      ld_ref     = params$ld_ref,
      output_dir = ldsc_output_dir,
      run_name   = params$run_name
    ),
    report_env
  )

  # ------------------------------------------------
  # 4) Save LDSC diagnostics
  # ------------------------------------------------
  safe_stage(
    "ldsc_diagnostics",
    function() save_ldsc_outputs(
      ldsc_object  = ldsc,
      output_prefix = file.path(ldsc_output_dir, params$run_name)
    ),
    report_env
  )

  # ------------------------------------------------
  # 5) Load manifest
  # ------------------------------------------------
  input_file_df <- safe_stage(
    "load_manifest",
    function() fread(params$input_manifest),
    report_env
  )

  # ------------------------------------------------
  # 6) GenomicPCA + GWAMA
  # ------------------------------------------------
  safe_stage(
    "genomicpca_gwama",
    function() run_genomicpca_gwama(
      input_file_df = input_file_df,
      ldsc          = ldsc,
      output_dir    = params$output_dir,
      approach      = params$approach
    ),
    report_env
  )

  # ------------------------------------------------
  # 7) Metadata capture
  # ------------------------------------------------
  safe_stage(
    "environment_capture",
    function() capture_environment(params$output_dir),
    report_env
  )

  safe_stage(
    "external_tools_capture",
    function() capture_external_tools(params$output_dir),
    report_env
  )

  message("✅ GenomicPCA pipeline completed successfully")

  invisible(finalize_pipeline_report(report_env))
}