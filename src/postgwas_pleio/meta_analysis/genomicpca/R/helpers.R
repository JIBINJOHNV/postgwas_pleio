###############################################################################
# CLI args → pipeline params
###############################################################################
###############################################################################
# CLI args → pipeline params (GenomicPCA minimal version)
###############################################################################
args_to_params <- function(args) {

  suppressPackageStartupMessages(library(glue))

  message("▶ Converting CLI args → pipeline params")

  # ---- Normalize paths ----
  input_manifest <- normalizePath(args$input_manifest, mustWork = TRUE)
  hm3            <- normalizePath(args$hm3, mustWork = TRUE)
  ld_ref         <- normalizePath(args$ld_ref, mustWork = TRUE)

  if (!dir.exists(ld_ref))
    stop("ld_ref must be a directory: ", ld_ref)

  output_dir <- normalizePath(args$output_dir, mustWork = FALSE)

  if (!dir.exists(output_dir))
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # ---- Derived directories (used by pipeline) ----
  ldsc_dir      <- file.path(output_dir, "3_ldsc_analysis")
  ldsc_input    <- file.path(ldsc_dir, "ldsc_input")
  ldsc_output   <- file.path(ldsc_dir, "ldsc_output")

  gpca_cor_dir  <- file.path(output_dir, "4_GenomicPCA_correlation")
  gpca_cov_dir  <- file.path(output_dir, "5_GenomicPCA_alternative_covariates")

  dir.create(ldsc_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(ldsc_input, recursive = TRUE, showWarnings = FALSE)
  dir.create(ldsc_output, recursive = TRUE, showWarnings = FALSE)
  dir.create(gpca_cor_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(gpca_cov_dir, recursive = TRUE, showWarnings = FALSE)

  # ---- Parameter object ----
  params <- list(
    input_manifest = input_manifest,
    output_dir     = output_dir,
    hm3            = hm3,
    ld_ref         = ld_ref,
    run_name       = args$run_name,
    approach       = args$approach,
    info_filter    = args$info_filter,
    maf_filter     = args$maf_filter,

    # directories used in pipeline
    ldsc_dir       = ldsc_dir,
    ldsc_input     = ldsc_input,
    ldsc_output    = ldsc_output,
    gpca_cor_dir   = gpca_cor_dir,
    gpca_cov_dir   = gpca_cov_dir,

    # metadata
    run_time       = Sys.time(),
    working_dir    = getwd(),
    host           = Sys.info()[["nodename"]]
  )

  message("✅ Parameter object created")

  return(params)
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



save_ldsc_outputs <- function(
    ldsc_object,
    output_prefix
){
    library(ggplot2)
    library(reshape2)
    # ------------------------------
    # Extract matrices
    # ------------------------------
    S_Stand <- ldsc_object$S_Stand
    S <- ldsc_object$S
    I <- ldsc_object$I
    # ------------------------------
    # Save matrices
    # ------------------------------
    write.csv(S_Stand, paste0(output_prefix,"_genetic_correlation_matrix.csv"))
    write.csv(S, paste0(output_prefix,"_genetic_covariance_matrix.csv"))
    write.csv(I, paste0(output_prefix,"_intercept_matrix.csv"))
    # ------------------------------
    # Heatmap of genetic correlations
    # ------------------------------
    cor_df <- reshape2::melt(S_Stand)
    p <- ggplot(cor_df, aes(Var1, Var2, fill=value)) +
        geom_tile() +
        scale_fill_gradient2(
            low="blue",mid="white", high="red", midpoint=0
        ) +
        theme_bw() +
        theme(axis.text.x = element_text(angle=45,hjust=1)) +
        labs( title="LDSC Genetic Correlation Matrix",  x="Trait",
            y="Trait", fill="rg" )
    ggsave(
        paste0(output_prefix,"_genetic_correlation_heatmap.png"), p,  width=6,  height=5 )

}


run_genomicpca_diagnostics <- function(
    matrix_input,
    trait_names,
    output_prefix
){
      library(ggplot2)
      # --------------------------------------------------
      # Eigen decomposition
      # --------------------------------------------------
      eig <- eigen(matrix_input)
      eigenvectors <- eig$vectors
      eigenvalues  <- eig$values
      # --------------------------------------------------
      # Explained variance
      # --------------------------------------------------
      explained_variance <- eigenvalues / sum(eigenvalues) * 100
      # --------------------------------------------------
      # Standardised PCA loadings (author formula)
      # L = V * sqrt(Lambda)
      # --------------------------------------------------
      loadings <- as.vector(eigenvectors %*% sqrt(diag(eigenvalues))[,1])
      names(loadings) <- trait_names
      # --------------------------------------------------
      # Save eigenvalues
      # --------------------------------------------------
      eigen_df <- data.frame(
          PC = paste0("PC", seq_along(eigenvalues)),
          Eigenvalue = eigenvalues,
          ExplainedVariance = explained_variance
      )
      write.csv(  eigen_df,
          paste0(output_prefix,"_pca_eigenvalues.csv"),  row.names = FALSE )
      # --------------------------------------------------
      # Save loadings
      # --------------------------------------------------
      loading_df <- data.frame( Trait = trait_names, PC1_Loading = loadings )
      write.csv( loading_df, paste0(output_prefix,"_pc1_loadings.csv"), row.names = FALSE )
      # --------------------------------------------------
      # Scree plot
      # --------------------------------------------------
      p1 <- ggplot(eigen_df, aes(x=PC, y=Eigenvalue, group=1)) +
            geom_point(size=3) +
            geom_line() +
            theme_bw() +
            labs( title="Genomic PCA Scree Plot",  x="Principal Component",  y="Eigenvalue" )
      ggsave(  paste0(output_prefix,"_scree_plot.png"), p1,  width=6, height=4 )
      # --------------------------------------------------
      # Explained variance plot
      # --------------------------------------------------
      p2 <- ggplot(eigen_df, aes(x=PC, y=ExplainedVariance)) +
            geom_bar(stat="identity") +
            theme_bw() +
            labs( title="Explained Variance",  x="Principal Component", y="Variance Explained (%)" )
      ggsave(
          paste0(output_prefix,"_variance_explained.png"), p2, width=6, height=4 )
      # --------------------------------------------------
      # Trait loading plot
      # --------------------------------------------------
      p3 <- ggplot(loading_df, aes(x=Trait, y=PC1_Loading)) +
            geom_bar(stat="identity") +
            theme_bw() +
            theme(axis.text.x = element_text(angle=45,hjust=1)) +
            labs( title="PC1 Trait Loadings", x="Trait",  y="Loading"  )
      ggsave(
          paste0(output_prefix,"_pc1_loadings_plot.png"),
          p3, width=7, height=5 )
      return(loadings)

}