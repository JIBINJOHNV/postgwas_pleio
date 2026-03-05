run_placo_pipeline <- function(
    input,
    placo_script="/opt/placo/PLACO_v0.2.0.R",
    output_csv = "placo_results.csv",
    method = c("placo", "placo.plus"),
    id_col = "ID",
    z_suffix = "_Z",
    p_suffix = "_P",
    p.threshold = 1e-4,
    remove_na = TRUE,
    verbose = TRUE
){

    method <- match.arg(method)

    # ===============================
    # LOAD PLACO
    # ===============================
    source(placo_script)

    # ===============================
    # LOAD INPUT
    # ===============================
    if(is.character(input)){
      if(verbose) cat("Reading input file\n")
      df <- data.table::fread(input)
    } else {
      df <- input
    }

    # ===============================
    # DETECT TRAITS
    # ===============================
    z_cols <- grep(paste0(z_suffix, "$"), names(df), value = TRUE)

    if(length(z_cols) != 2){
      stop("Exactly 2 *_Z columns required")
    }

    traits <- sub(paste0(z_suffix, "$"), "", z_cols)
    p_cols <- paste0(traits, p_suffix)

    if(!all(p_cols %in% names(df))){
      stop("Matching *_P columns missing")
    }

    if(verbose){
      cat("Detected traits:\n")
      print(traits)
    }

    # ===============================
    # BUILD MATRICES
    # ===============================
    Z.matrix <- as.matrix(df[, z_cols, with = FALSE])
    P.matrix <- as.matrix(df[, p_cols, with = FALSE])
    ids <- df[[id_col]]

    # ===============================
    # NA FILTER
    # ===============================
    if(remove_na){
      keep <- complete.cases(Z.matrix) & complete.cases(P.matrix)
      Z.matrix <- Z.matrix[keep, , drop=FALSE]
      P.matrix <- P.matrix[keep, , drop=FALSE]
      ids <- ids[keep]
    }

    if(verbose) cat("Variants retained:", length(ids), "\n")

    # ===============================
    # PARAMETER ESTIMATION
    # ===============================
    if(verbose) cat("Estimating VarZ\n")
    VarZ <- var.placo(Z.matrix, P.matrix, p.threshold = p.threshold)

    CorZ <- NULL
    if(method == "placo.plus"){
      if(verbose) cat("Estimating CorZ\n")
      CorZ <- cor.pearson(Z.matrix, P.matrix,
                          p.threshold = p.threshold,
                          returnMatrix = FALSE)
    }

    # ===============================
    # GENOMEWIDE RUN
    # ===============================
    if(verbose) cat("Running", method, "\n")

    res <- t(
      apply(Z.matrix, 1, function(z){

        if(method == "placo"){
          out <- placo(Z = z, VarZ = VarZ)
          return(c(out$T.placo, out$p.placo))
        }

        if(method == "placo.plus"){
          out <- placo.plus(Z = z, VarZ = VarZ, CorZ = CorZ)
          return(c(out$T.placo.plus, out$p.placo.plus))
        }

      })
    )

    colnames(res) <- c("PLACO_stat", "PLACO_p")

    # ===============================
    # OUTPUT
    # ===============================
    out_df <- data.frame(
      ID = ids,
      PLACO_stat = res[,1],
      PLACO_p = res[,2]
    )

    if(!is.null(output_csv)){
      if(verbose) cat("Writing output:", output_csv, "\n")
      data.table::fwrite(out_df, output_csv)
    }

    if(verbose) cat("Done\n")


  # ------------------------------------------------
  # 8) metadata
  # ------------------------------------------------
  safe_stage("environment_capture", function() capture_environment(params$output_dir), report_env)
  safe_stage("external_tools_capture", function() capture_external_tools(params$output_dir), report_env)

  message("✅ Pipeline completed successfully")

  invisible(finalize_pipeline_report(report_env))


  return(out_df)

}