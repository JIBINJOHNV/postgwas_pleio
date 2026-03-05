# ============================================================
# Production-grade PLACO / PLACO+ pipeline (modular + parallel)
# - Stage 1 parallel: (total cores - 2) by default
# - Stage 2 retry parallel: 40% of total cores by default
# - Sources PLACO script inside EACH worker (critical on macOS/Ubuntu PSOCK)
# - No hardcoding: uses id_col, suffixes, etc.
# - Copy-paste friendly (single block)
#
# Requirements:
#   library(data.table)
#   library(parallel)
#   placo_script must define: var.placo, cor.pearson, placo, placo.plus
# ============================================================

# ✅ 2️⃣ Load Input
load_input <- function(input, verbose){
  if (is.character(input)) {
    if (isTRUE(verbose)) cat("Reading input file\n")
    return(data.table::fread(input))
  }
  return(input)
}

# ✅ 3️⃣ Detect Traits
detect_traits <- function(df, z_suffix, p_suffix, verbose){

  z_cols <- grep(paste0(z_suffix, "$"), names(df), value = TRUE)

  if (length(z_cols) != 2)
    stop("Exactly 2 *_Z columns required")

  traits <- sub(paste0(z_suffix, "$"), "", z_cols)
  p_cols <- paste0(traits, p_suffix)

  if (!all(p_cols %in% names(df)))
    stop("Matching *_P columns missing")

  if (isTRUE(verbose)) {
    cat("Detected traits:\n")
    print(traits)
  }

  list(z_cols = z_cols, p_cols = p_cols, traits = traits)
}

# ✅ 4️⃣ Build Matrices
build_matrices <- function(df, trait_info, id_col){

  if (!id_col %in% names(df))
    stop("id_col not found in input: ", id_col)

  list(
    Z = as.matrix(df[, trait_info$z_cols, with = FALSE]),
    P = as.matrix(df[, trait_info$p_cols, with = FALSE]),
    ids = df[[id_col]]
  )
}

# ✅ 5️⃣ NA Filter
filter_na <- function(matrices, remove_na, verbose){

  if (!isTRUE(remove_na)) {
    if (isTRUE(verbose)) cat("Variants retained:", length(matrices$ids), "\n")
    return(matrices)
  }

  keep <- complete.cases(matrices$Z) & complete.cases(matrices$P) & !is.na(matrices$ids)

  matrices$Z <- matrices$Z[keep, , drop = FALSE]
  matrices$P <- matrices$P[keep, , drop = FALSE]
  matrices$ids <- matrices$ids[keep]

  if (isTRUE(verbose))
    cat("Variants retained:", length(matrices$ids), "\n")

  if (nrow(matrices$Z) == 0)
    stop("No variants left after NA filtering.")

  matrices
}

# ✅ 6️⃣ Parameter Estimation
estimate_parameters <- function(Z, P, method, p.threshold, verbose){

  if (isTRUE(verbose)) cat("Estimating VarZ\n")
  VarZ <- var.placo(Z, P, p.threshold = p.threshold)

  CorZ <- NULL
  if (method == "placo.plus") {
    if (isTRUE(verbose)) cat("Estimating CorZ\n")
    CorZ <- cor.pearson(
      Z, P,
      p.threshold = p.threshold,
      returnMatrix = FALSE
    )
  }

  list(VarZ = VarZ, CorZ = CorZ)
}

# ✅ 7️⃣ Core Strategy
compute_core_strategy <- function(cores_stage1, cores_stage2, verbose){

  total <- parallel::detectCores()

  s1 <- if (is.null(cores_stage1)) {
    max(1, total - 2)
  } else {
    max(1, min(cores_stage1, total - 1))
  }

  s2 <- if (is.null(cores_stage2)) {
    max(1, floor(total * 0.4))
  } else {
    max(1, min(cores_stage2, total - 1))
  }

  s2 <- min(s2, s1)

  if (isTRUE(verbose)){
    cat("Stage 1 cores:", s1, "\n")
    cat("Stage 2 cores:", s2, "\n")
  }

  list(stage1 = s1, stage2 = s2, total = total)
}

# ✅ 8️⃣ Parallel Stage Runner (critical: source placo_script inside workers)
run_parallel_stage <- function(
    Z,
    params,
    method,
    cores,
    abs_tol,
    placo_script,
    verbose,
    stage_label = "Stage"
){

  if (isTRUE(verbose)) cat("Running", stage_label, "\n")

  safe_test <- function(z){
    tryCatch({
      if (method == "placo") {
        out <- placo(Z = z, VarZ = params$VarZ, AbsTol = abs_tol)
        c(out$T.placo, out$p.placo)
      } else {
        out <- placo.plus(
          Z = z,
          VarZ = params$VarZ,
          CorZ = params$CorZ,
          AbsTol = abs_tol
        )
        c(out$T.placo.plus, out$p.placo.plus)
      }
    }, error = function(e) {
      c(NA_real_, NA_real_)
    })
  }

  # Sequential fallback
  if (cores <= 1) {
    return(t(apply(Z, 1, safe_test)))
  }

  cl <- parallel::makeCluster(cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  # Make placo_script visible to workers, then source it there
  parallel::clusterExport(cl, varlist = c("placo_script"), envir = environment())

  parallel::clusterEvalQ(cl, {
    source(placo_script)
    NULL
  })

  # Export parameters + function after sourcing
  parallel::clusterExport(
    cl,
    varlist = c("params", "method", "abs_tol", "safe_test"),
    envir = environment()
  )

  res <- t(parallel::parApply(cl, Z, 1, safe_test))
  return(res)
}

# ✅ 9️⃣ Retry Stage (calls run_parallel_stage with fewer cores + relaxed AbsTol)
retry_failed_stage <- function(
    Z,
    res,
    params,
    method,
    cores,
    abs_tol_retry,
    placo_script,
    verbose
){

  fail_idx <- which(is.na(res[, 2]))

  if (length(fail_idx) == 0)
    return(res)

  if (isTRUE(verbose))
    cat("Retrying", length(fail_idx), "failed SNPs\n")

  Z_retry <- Z[fail_idx, , drop = FALSE]

  retry_res <- run_parallel_stage(
    Z = Z_retry,
    params = params,
    method = method,
    cores = cores,
    abs_tol = abs_tol_retry,
    placo_script = placo_script,
    verbose = verbose,
    stage_label = "Stage 2 (retry)"
  )

  res[fail_idx, ] <- retry_res
  res
}

# ✅ 1️⃣0️⃣ Assemble Output
assemble_output <- function(ids, res, id_col){

  out <- data.frame(
    ids,
    PLACO_stat = res[, 1],
    PLACO_p = res[, 2]
  )

  names(out)[1] <- id_col
  out
}

# ✅ 1️⃣1️⃣ Write Output
write_output <- function(out_df, output_csv, verbose){

  if (is.null(output_csv)) return(invisible(NULL))

  if (isTRUE(verbose))
    cat("Writing output:", output_csv, "\n")

  data.table::fwrite(out_df, output_csv)
  invisible(NULL)
}

# ✅ 1️⃣ Main Orchestrator
run_placo_pipeline <- function(
    input,
    placo_script = "/opt/placo/PLACO_v0.2.0.R",
    output_csv = "placo_results.csv",
    method = c("placo", "placo.plus"),
    id_col = "ID",
    z_suffix = "_Z",
    p_suffix = "_P",
    p.threshold = 1e-4,
    remove_na = TRUE,
    cores_stage1 = NULL,
    cores_stage2 = NULL,
    retry_failed = TRUE,
    abs_tol = .Machine$double.eps^0.8,
    abs_tol_retry = 1e-10,
    verbose = TRUE
){

  method <- match.arg(method)

  # Source once in master (workers will source again inside run_parallel_stage)
  source(placo_script)

  df <- load_input(input, verbose)

  trait_info <- detect_traits(df, z_suffix, p_suffix, verbose)

  matrices <- build_matrices(df, trait_info, id_col)

  matrices <- filter_na(matrices, remove_na, verbose)

  params <- estimate_parameters(
    matrices$Z, matrices$P, method, p.threshold, verbose
  )

  cores <- compute_core_strategy(
    cores_stage1, cores_stage2, verbose
  )

  res <- run_parallel_stage(
    Z = matrices$Z,
    params = params,
    method = method,
    cores = cores$stage1,
    abs_tol = abs_tol,
    placo_script = placo_script,
    verbose = verbose,
    stage_label = "Stage 1"
  )

  if (isTRUE(retry_failed)) {
    res <- retry_failed_stage(
      Z = matrices$Z,
      res = res,
      params = params,
      method = method,
      cores = cores$stage2,
      abs_tol_retry = abs_tol_retry,
      placo_script = placo_script,
      verbose = verbose
    )
  }

  out_df <- assemble_output(matrices$ids, res, id_col)

  write_output(out_df, output_csv, verbose)

  if (isTRUE(verbose)) {
    cat("Final failed SNPs:", sum(is.na(out_df$PLACO_p)), "\n")
    cat("Done\n")
  }

  return(out_df)
}