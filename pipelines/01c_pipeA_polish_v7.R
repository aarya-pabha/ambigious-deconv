# 02_run_sol4_snap_polish_v7.R
# Calibrated Step-Drift Polish Engine (d_res <= 0.85)
# Abundance-Calibrated Screening, Pass A -> Pass B Order, Pareto Checkpointing
for (p in c("phase7/R/helper-functions.R", "../phase7/R/helper-functions.R", "~/phase7/R/helper-functions.R")) if (file.exists(p)) { source(p); break }
for (p in c("phase7/R/nmf-functions_v2.R", "../phase7/R/nmf-functions_v2.R", "~/phase7/R/nmf-functions_v2.R")) if (file.exists(p)) { source(p); break }
for (p in c("phase7/R/simulation-functions.R", "../phase7/R/simulation-functions.R", "~/phase7/R/simulation-functions.R")) if (file.exists(p)) { source(p); break }
suppressPackageStartupMessages(library(quadprog))

fit_fcls_h <- function(P_obs, W) {
  K <- ncol(W); n <- ncol(P_obs)
  if (K == 1) return(matrix(1.0, 1, n))
  Dmat <- crossprod(W) + diag(1e-6, K)
  Amat <- cbind(1, diag(K)); bvec <- c(1, rep(0, K))
  H <- vapply(seq_len(n), function(j) {
    qp <- tryCatch(solve.QP(Dmat, crossprod(W, P_obs[, j]), Amat, bvec, meq = 1),
                   error = function(e) list(solution = rep(1 / K, K)))
    sol <- pmax(0, qp$solution); sol / sum(sol)
  }, numeric(K))
  if (K == 1) matrix(H, 1) else H
}

score_binom <- function(X, W, H, eps = 1e-3) {
  P <- pmax(eps, pmin(1 - eps, W %*% H))
  sum(dbinom(X, size = 100, prob = P, log = TRUE))
}

compute_row_deviance <- function(x_row, p_row, eps = 1e-3) {
  p_clamped <- pmax(eps, pmin(1 - eps, p_row))
  mu_val <- 100 * p_clamped
  d <- numeric(length(x_row))
  
  pos_x <- x_row > 0
  d[pos_x] <- d[pos_x] + 2 * x_row[pos_x] * log(x_row[pos_x] / mu_val[pos_x])
  
  pos_nx <- (100 - x_row) > 0
  d[pos_nx] <- d[pos_nx] + 2 * (100 - x_row[pos_nx]) * log((100 - x_row[pos_nx]) / (100 - mu_val[pos_nx]))
  sum(d)
}

polish_bits_v7 <- function(X, P_obs, W_init, U_basis, mu_vec, eps = 1e-3, max_epochs = 3) {
  W <- W_init
  M <- nrow(W); K <- ncol(W); n <- ncol(P_obs)
  h_diag <- 1.0 - rowSums(U_basis^2)

  # Pre-clamp monomorphic loci
  inv_0 <- which(rowSums(X) == 0)
  inv_1 <- which(rowSums(X) == n * 100)
  W[inv_0, ] <- 0L; W[inv_1, ] <- 1L

  H <- fit_fcls_h(P_obs, W)
  init_ll <- score_binom(X, W, H, eps)
  init_mse <- mean((P_obs - W %*% H)^2)

  # Initial projection residual matrix
  res_mat <- (W - mu_vec) - U_basis %*% crossprod(U_basis, W - mu_vec)

  # Multi-epoch Pareto checkpoint tracking
  W_best <- W
  H_best <- H
  ll_best <- init_ll
  mse_best <- init_mse

  pairs <- if (K > 1) combn(K, 2, simplify = FALSE) else list()

  for (epoch in seq_len(max_epochs)) {
    flips_in_epoch <- 0
    P_hat <- W %*% H
    mean_h <- rowMeans(H)

    # -------------------------------------------------------------
    # PASS A: Abundance-Adaptive 1-Bit BCD (Anchors Genotypes)
    # -------------------------------------------------------------
    for (m in seq_len(M)) {
      if (m %in% inv_0 || m %in% inv_1) next
      x_m <- X[m, ]
      row_p <- P_hat[m, ]
      ll_base <- sum(dbinom(x_m, 100, pmax(eps, pmin(1 - eps, row_p)), log = TRUE))

      best_d1 <- 0; best_k1 <- 0
      for (k in seq_len(K)) {
        if (mean_h[k] < 0.01) next
        d_bit <- 1L - 2L * W[m, k]
        p_trial <- row_p + d_bit * H[k, ]
        dll <- sum(dbinom(x_m, 100, pmax(eps, pmin(1 - eps, p_trial)), log = TRUE)) - ll_base
        tau_k <- max(2.5, min(5.0, 4.5 * sqrt(mean_h[k] / 0.18)))

        if (dll > tau_k && dll > best_d1) {
          # Calibrated step-drift threshold: captures up to 92nd %ile of true mutations
          d_res <- 2 * d_bit * res_mat[m, k] + h_diag[m]
          if (d_res <= 0.85) {
            best_d1 <- dll; best_k1 <- k
          }
        }
      }

      if (best_k1 > 0) {
        d_bit <- 1L - 2L * W[m, best_k1]
        W[m, best_k1] <- 1L - W[m, best_k1]
        row_p <- row_p + d_bit * H[best_k1, ]
        P_hat[m, ] <- row_p
        res_mat[, best_k1] <- (W[, best_k1] - mu_vec) - U_basis %*% crossprod(U_basis, W[, best_k1] - mu_vec)
        flips_in_epoch <- flips_in_epoch + 1
      }
    }

    # -------------------------------------------------------------
    # PASS B: Deviance-Gated Pairwise Moves on Residual Saddles
    # -------------------------------------------------------------
    for (m in seq_len(M)) {
      if (m %in% inv_0 || m %in% inv_1) next
      x_m <- X[m, ]
      row_p <- P_hat[m, ]

      # Abundance-calibrated screening gate
      active_strains <- which(W[m, ] == 1L)
      max_h_active <- if (length(active_strains) > 0) max(mean_h[active_strains]) else 0.0
      gate_thresh <- n + 1.0 * sqrt(2 * n) * max(0.4, (1.0 - 0.5 * max_h_active))

      row_dev <- compute_row_deviance(x_m, row_p, eps)
      if (row_dev <= gate_thresh) next

      ll_base <- sum(dbinom(x_m, 100, pmax(eps, pmin(1 - eps, row_p)), log = TRUE))
      best_d2 <- 0; best_pair <- NULL; best_trans <- NULL

      for (pr in pairs) {
        k1 <- pr[1]; k2 <- pr[2]
        if (mean_h[k1] < 0.01 && mean_h[k2] < 0.01) next
        pair_h <- mean_h[k1] + mean_h[k2]

        if (W[m, k1] != W[m, k2]) {
          # Allele Swap: (1, 0) <-> (0, 1)
          dp <- (W[m, k2] - W[m, k1]) * H[k1, ] + (W[m, k1] - W[m, k2]) * H[k2, ]
          dll <- sum(dbinom(x_m, 100, pmax(eps, pmin(1 - eps, row_p + dp)), log = TRUE)) - ll_base
          tau_swap <- max(4.0, min(5.2, 4.2 * sqrt(pair_h / 0.36)))

          if (dll > tau_swap && dll > best_d2) {
            d_bit1 <- 1L - 2L * W[m, k1]
            d_bit2 <- 1L - 2L * W[m, k2]
            d_res1 <- 2 * d_bit1 * res_mat[m, k1] + h_diag[m]
            d_res2 <- 2 * d_bit2 * res_mat[m, k2] + h_diag[m]

            # Calibrated step-drift threshold
            if (d_res1 <= 0.85 && d_res2 <= 0.85) {
              best_d2 <- dll; best_pair <- pr; best_trans <- "swap"
            }
          }
        } else {
          # Joint Co-mutation: (0, 0) <-> (1, 1)
          s <- 1L - 2L * W[m, k1]
          dp <- s * (H[k1, ] + H[k2, ])
          dll <- sum(dbinom(x_m, 100, pmax(eps, pmin(1 - eps, row_p + dp)), log = TRUE)) - ll_base
          tau_joint <- max(6.0, min(8.0, 6.8 * sqrt(pair_h / 0.36)))

          if (dll > tau_joint && dll > best_d2) {
            d_res1 <- 2 * s * res_mat[m, k1] + h_diag[m]
            d_res2 <- 2 * s * res_mat[m, k2] + h_diag[m]

            # Calibrated step-drift threshold
            if (d_res1 <= 0.85 && d_res2 <= 0.85) {
              best_d2 <- dll; best_pair <- pr; best_trans <- "joint"
            }
          }
        }
      }

      if (!is.null(best_pair)) {
        k1 <- best_pair[1]; k2 <- best_pair[2]
        W[m, k1] <- 1L - W[m, k1]; W[m, k2] <- 1L - W[m, k2]
        row_p <- W[m, , drop = FALSE] %*% H
        P_hat[m, ] <- row_p
        res_mat[, k1] <- (W[, k1] - mu_vec) - U_basis %*% crossprod(U_basis, W[, k1] - mu_vec)
        res_mat[, k2] <- (W[, k2] - mu_vec) - U_basis %*% crossprod(U_basis, W[, k2] - mu_vec)
        flips_in_epoch <- flips_in_epoch + 2
      }
    }

    if (flips_in_epoch == 0) break

    # Re-solve H once per macro-epoch
    H <- fit_fcls_h(P_obs, W)
    curr_ll <- score_binom(X, W, H, eps)
    curr_mse <- mean((P_obs - W %*% H)^2)

    # Multi-epoch Pareto checkpointing: anchor to best running state
    if (curr_ll > ll_best && curr_mse <= (mse_best * 1.01)) {
      W_best <- W
      H_best <- H
      ll_best <- curr_ll
      mse_best <- curr_mse
    } else {
      W <- W_best
      H <- H_best
      break
    }
  }

  list(W = W_best, H = H_best, loglik = ll_best, mse = mse_best, status = "accepted")
}

run_manifest <- function(manifest_path = "master_manifest.csv") {
  if (!file.exists(manifest_path) && file.exists("../master_manifest.csv")) manifest_path <- "../master_manifest.csv"
  manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
  args <- commandArgs(trailingOnly = TRUE)
  tasks <- if (length(args) > 0) as.integer(args[1]) else seq_len(nrow(manifest))

  for (idx in tasks) {
    row <- manifest[idx, ]
    out_dir <- row$out_dir
    if (!dir.exists(out_dir) && dir.exists(file.path("..", out_dir))) out_dir <- file.path("..", out_dir)

    out_file <- file.path(out_dir, "summary_hunter_sol4_snap_polish_v7_v2.csv")
    if (file.exists(out_file)) next

    w_sub_p <- file.path(out_dir, sprintf("W_alt_sol4_snap_subset_1_r%s.csv", row$r_true))
    w_true_p <- file.path(out_dir, "W_true.csv")
    if (!file.exists(w_sub_p) || !file.exists(w_true_p)) next

    W_sub <- as.matrix(read.csv(w_sub_p))
    W_true <- as.matrix(read.csv(w_true_p))
    K <- as.integer(row$r_true); n <- as.integer(row$samples)

    set.seed(row$seed)
    alpha <- as.numeric(unlist(strsplit(row$alpha_str, ",")))
    H_true <- t(matrix(rgamma(n * K, alpha, 1), ncol = K, byrow = TRUE))
    H_true <- H_true / rep(colSums(H_true), each = K)
    X <- matrix(rbinom(365 * n, size = 100, prob = as.vector(W_true %*% H_true)), nrow = 365)
    P_obs <- X / 100

    fit_base <- fit_affine_subspace(X, rep(100, n), K, transform = "none", eps_var = 1e-6)
    fit_polished <- polish_bits_v7(X, P_obs, W_sub, fit_base$U, fit_base$mu)

    write.csv(fit_polished$W, file.path(out_dir, sprintf("W_alt_sol4_snap_polish_v7_1_r%d.csv", K)), row.names = FALSE)
    write.csv(data.frame(solution_id = 1, loglik = fit_polished$loglik,
                         mse_score = fit_polished$mse, num_strains = K,
                         status = fit_polished$status),
              out_file, row.names = FALSE)
    cat(sprintf("[DONE] Task %d: Status=%s, LL=%.1f, MSE=%.5f\n",
                idx, fit_polished$status, fit_polished$loglik, fit_polished$mse))
  }
}

if (sys.nframe() == 0) run_manifest()
