# 02_run_snapped_subset_sol4.R
for (p in c("phase7/R/helper-functions.R", "../phase7/R/helper-functions.R", "~/phase7/R/helper-functions.R")) if (file.exists(p)) { source(p); break }
for (p in c("phase7/R/nmf-functions_v2.R", "../phase7/R/nmf-functions_v2.R", "~/phase7/R/nmf-functions_v2.R")) if (file.exists(p)) { source(p); break }
for (p in c("phase7/R/simulation-functions.R", "../phase7/R/simulation-functions.R", "~/phase7/R/simulation-functions.R")) if (file.exists(p)) { source(p); break }

suppressPackageStartupMessages(library(quadprog))

# FCLS via Quadratic Programming
fit_fcls_qp <- function(P_obs, W) {
  K <- ncol(W); n <- ncol(P_obs)
  Dmat <- crossprod(W) + diag(1e-7, K)
  Amat <- cbind(1, diag(K))
  bvec <- c(1, rep(0, K))

  H <- vapply(seq_len(n), function(j) {
    qp <- tryCatch(
      solve.QP(Dmat, crossprod(W, P_obs[, j]), Amat, bvec, meq = 1),
      error = function(e) list(solution = rep(1 / K, K))
    )
    sol <- pmax(0, qp$solution)
    sol / sum(sol)
  }, numeric(K))

  if (K == 1) H <- matrix(H, nrow = 1)
  list(W = W, H = H, mse = mean((P_obs - W %*% H)^2))
}

# Vectorized 1-bit & 2-bit descent to hypercube facet
snap_candidate <- function(w, U, mu, M, h_diag) {
  w_curr <- as.integer(w)
  for (it in 1:10) {
    r <- (w_curr - mu) - U %*% crossprod(U, w_curr - mu)
    d1 <- 2 * (1L - 2L * w_curr) * as.vector(r) + h_diag
    m1 <- which.min(d1)
    if (d1[m1] < -1e-5) { w_curr[m1] <- 1L - w_curr[m1]; next }

    s <- (1L - 2L * w_curr)
    d2 <- outer(d1, d1, "+") + 2 * outer(s, s) * M
    diag(d2) <- Inf
    if (min(d2) < -1e-5) {
      idx <- which(d2 == min(d2), arr.ind = TRUE)[1, ]
      w_curr[idx[1]] <- 1L - w_curr[idx[1]]; w_curr[idx[2]] <- 1L - w_curr[idx[2]]
      next
    }
    break
  }
  w_curr
}

snap_pool_to_basis <- function(W_pool, U, mu) {
  M <- diag(1.0, nrow(W_pool)) - tcrossprod(U)
  h_diag <- diag(M)
  snapped <- apply(W_pool, 2, function(w) snap_candidate(w, U, mu, M, h_diag))
  W_all <- cbind(W_pool, snapped)
  W_all[, !duplicated(t(W_all)), drop = FALSE]
}

# Forward selection + 1-pass swap
select_best_k_subset <- function(P_obs, W_pool, K_target) {
  P_cands <- ncol(W_pool)
  if (P_cands <= K_target) return(fit_fcls_qp(P_obs, W_pool))

  sel <- integer(0)
  for (step in seq_len(K_target)) {
    rem <- setdiff(seq_len(P_cands), sel)
    mses <- vapply(rem, function(c) fit_fcls_qp(P_obs, W_pool[, c(sel, c), drop = FALSE])$mse, numeric(1))
    sel <- c(sel, rem[which.min(mses)])
  }

  for (i in seq_along(sel)) {
    rem <- setdiff(seq_len(P_cands), sel)
    curr_mse <- fit_fcls_qp(P_obs, W_pool[, sel, drop = FALSE])$mse
    mses <- vapply(rem, function(c) {
      test <- sel; test[i] <- c
      fit_fcls_qp(P_obs, W_pool[, test, drop = FALSE])$mse
    }, numeric(1))
    if (min(mses) < curr_mse - 1e-6) sel[i] <- rem[which.min(mses)]
  }

  fit_fcls_qp(P_obs, W_pool[, sel, drop = FALSE])
}

# Pipeline execution
run_pipeline <- function(manifest_file = "master_manifest.csv") {
  if (!file.exists(manifest_file)) return(invisible(NULL))
  manifest <- read.csv(manifest_file, stringsAsFactors = FALSE)
  cat(sprintf("[INFO] Processing %d tasks with sol4_snap_subset\n", nrow(manifest)))

  for (idx in seq_len(nrow(manifest))) {
    row <- manifest[idx, ]; out_dir <- row$out_dir; K <- as.integer(row$r_true)
    sum_p <- file.path(out_dir, "summary_hunter_sol4.csv")
    w_true_p <- file.path(out_dir, "W_true.csv")
    if (!file.exists(sum_p) || file.info(sum_p)$size == 0 || !file.exists(w_true_p)) next

    summary_df <- read.csv(sum_p)
    if (nrow(summary_df) == 0) next

    # Load sol4 pool
    w_files <- file.path(out_dir, sprintf("W_alt_sol4_%s_r%d.csv", summary_df$solution_id, K))
    w_files <- w_files[file.exists(w_files) & file.info(w_files)$size > 0]
    if (length(w_files) == 0) next

    W_pool <- do.call(cbind, lapply(w_files, function(f) as.matrix(read.csv(f))))
    W_pool <- W_pool[, !duplicated(t(W_pool)), drop = FALSE]

    # Reconstruct P_obs
    W_true <- as.matrix(read.csv(w_true_p))
    n <- as.integer(row$samples); set.seed(row$seed)
    alpha <- as.numeric(unlist(strsplit(row$alpha_str, ",")))
    H_true <- t(matrix(rgamma(n * K, alpha, 1), ncol = K, byrow = TRUE))
    H_true <- H_true / rep(colSums(H_true), each = K)
    X <- matrix(rbinom(365 * n, size = 100, prob = as.vector(W_true %*% H_true)), nrow = 365)
    P_obs <- X / 100

    # Fit canonical affine subspace using phase7 helper
    fit_base <- fit_affine_subspace(X, rep(100, n), K, transform = "none", eps_var = 1e-6)

    # Micro-snap and stepwise extract K
    W_snapped <- snap_pool_to_basis(W_pool, fit_base$U, fit_base$mu)
    fit <- select_best_k_subset(P_obs, W_snapped, K_target = K)

    # Save outputs
    write.csv(W_snapped, file.path(out_dir, sprintf("W_pool_snapped_r%d.csv", K)), row.names = FALSE)
    write.csv(fit$W, file.path(out_dir, sprintf("W_alt_sol4_snap_subset_1_r%d.csv", K)), row.names = FALSE)
    write.csv(data.frame(solution_id = 1, source = "sol4_snap_subset", mse_score = fit$mse,
                         avg_distance = 0.0, max_distance = 0.0, num_strains = ncol(fit$W)),
              file.path(out_dir, "summary_hunter_sol4_snap_subset.csv"), row.names = FALSE)

    if (idx %% 20 == 0) cat(sprintf("  Processed %d / %d tasks\n", idx, nrow(manifest)))
  }
}

if (sys.nframe() == 0) run_pipeline()
