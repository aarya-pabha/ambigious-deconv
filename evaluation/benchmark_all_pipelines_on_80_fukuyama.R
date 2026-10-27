#!/usr/bin/env Rscript
# ==============================================================================
# File: benchmark_all_pipelines_on_80_fukuyama.R
# Purpose: Completely self-contained benchmark of ALL 5 PIPELINE GENERATIONS
#          on the 80 Ambiguous Fukuyama Manifold Tasks (K in {5, 7, 8, 10}, 20 seeds):
#          1. PipeA (sol4_v7 3-Constraint Hunter)
#          2. PipeB (unified_15plus v2)
#          3. v3 (Unbiased 20-SVD + Snapped Pool + Gated Polish)
#          4. v4_Std (Pure Snapped Pool + 4-Start Hybrid Exchange)
#          5. v4_SubReg (Subspace-Regularized Exchange)
# Usage:   Rscript benchmark_all_pipelines_on_80_fukuyama.R <task_id (1-80)>
# ==============================================================================

suppressPackageStartupMessages({
  library(quadprog)
  library(MASS)
})

args <- commandArgs(trailingOnly = TRUE)
task_id <- if (length(args) > 0) as.integer(args[1]) else 1L

# Task Mapping across the 80 Fukuyama Ambiguity Tasks
K_map    <- c(rep(5L, 20), rep(7L, 20), rep(8L, 20), rep(10L, 20))
seed_map <- rep(202601L:202620L, 4)
K_val    <- K_map[task_id]
seed_val <- seed_map[task_id]
rep_idx  <- ((task_id - 1L) %% 20L) + 1L

# ------------------------------------------------------------------------------
# 1. CORE NUMERICAL SOLVERS & UTILITIES
# ------------------------------------------------------------------------------

fit_fcls_qp <- function(P_obs, W) {
  K <- ncol(W); n <- ncol(P_obs)
  Dmat <- crossprod(W) + diag(1e-8, K)
  Amat <- cbind(1, diag(K)); bvec <- c(1, rep(0, K))
  H <- matrix(0, nrow = K, ncol = n)
  for (j in seq_len(n)) {
    dvec <- as.vector(crossprod(W, P_obs[, j]))
    qp <- tryCatch(solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 1),
                   error = function(e) list(solution = rep(1 / K, K)))
    sol <- pmax(0, qp$solution)
    H[, j] <- if (sum(sol) > 0) sol / sum(sol) else rep(1 / K, K)
  }
  list(W = W, H = H, mse = mean((P_obs - W %*% H)^2))
}

fit_fcls_h <- function(P_obs, W) {
  K <- ncol(W); n <- ncol(P_obs)
  if (K == 1) return(matrix(1.0, 1, n))
  Dmat <- crossprod(W) + diag(1e-6, K)
  Amat <- cbind(1, diag(K)); bvec <- c(1, rep(0, K))
  vapply(seq_len(n), function(j) {
    dvec <- as.vector(crossprod(W, P_obs[, j]))
    qp <- tryCatch(solve.QP(Dmat, dvec, Amat, bvec, meq = 1),
                   error = function(e) list(solution = rep(1 / K, K)))
    sol <- pmax(0, qp$solution)
    if (sum(sol) > 0) sol / sum(sol) else rep(1 / K, K)
  }, numeric(K))
}

hungarian_dists <- function(cost_mat) {
  K <- nrow(cost_mat)
  if (requireNamespace("clue", quietly = TRUE)) {
    m_c <- as.vector(clue::solve_LSAP(cost_mat))
    return(cost_mat[cbind(seq_len(K), m_c)])
  } else if (requireNamespace("lpSolve", quietly = TRUE)) {
    sol <- lpSolve::lp.assign(cost_mat)$solution
    m_c <- apply(sol, 1, which.max)
    return(cost_mat[cbind(seq_len(K), m_c)])
  } else {
    matched <- integer(K); used <- logical(K)
    for (st in seq_len(K)) {
      bv <- Inf; bi <- 0L; bj <- 0L
      for (r in seq_len(K)) if (matched[r] == 0L) for (c in seq_len(K)) if (!used[c] && cost_mat[r, c] < bv) {
        bv <- cost_mat[r, c]; bi <- r; bj <- c
      }
      matched[bi] <- bj; used[bj] <- TRUE
    }
    return(cost_mat[cbind(seq_len(K), matched)])
  }
}

rel_vol <- function(W_cand, W_gt) {
  A1 <- rbind(1.0, W_cand); A0 <- rbind(1.0, W_gt)
  d1 <- det(crossprod(A1)); d0 <- det(crossprod(A0))
  if (is.na(d1) || is.na(d0) || d1 <= 0 || d0 <= 0) return(0.0)
  sqrt(d1 / d0)
}

fit_weighted_affine_subspace <- function(X, N_vec, K) {
  m <- nrow(X); n <- ncol(X)
  P <- X / matrix(rep(N_vec, each = m), nrow = m, ncol = n)
  mu <- rowMeans(P)
  w <- rexp(m, rate = 1.0); sqrt_w <- sqrt(w / mean(w))
  P_weighted <- sweep(sweep(P, 1, mu, "-"), 1, sqrt_w, "*")
  sv <- svd(P_weighted, nu = K - 1L, nv = 0)
  U_unweighted <- sweep(sv$u, 1, sqrt_w, "/")
  U <- qr.Q(qr(U_unweighted))[, seq_len(K - 1L), drop = FALSE]
  list(U = U, mu = mu)
}

select_pivots_maxvol <- function(U, K, max_iter = 5) {
  d <- max(1L, K - 1L); M <- nrow(U)
  pivots <- qr(t(U), tol = 1e-7)$pivot[seq_len(d)]
  U_R <- U[pivots, , drop = FALSE]
  curr_rcond <- rcond(U_R); if (is.na(curr_rcond)) curr_rcond <- 0
  cand_pool <- order(rowSums(U^2), decreasing = TRUE)[seq_len(min(M, 40L))]
  for (iter in seq_len(max_iter)) {
    improved <- FALSE
    for (i in seq_len(d)) {
      for (cand in cand_pool) {
        if (cand %in% pivots) next
        trial_pivots <- pivots; trial_pivots[i] <- cand
        trial_rcond <- rcond(U[trial_pivots, , drop = FALSE])
        if (!is.na(trial_rcond) && trial_rcond > curr_rcond * 1.20) {
          pivots <- trial_pivots; curr_rcond <- trial_rcond; improved <- TRUE; break
        }
      }
    }
    if (!improved) break
  }
  pivots
}

evaluate_set_difference <- function(W_new, W_old) {
  K_new <- ncol(W_new); K_old <- ncol(W_old)
  dist_matrix <- matrix(0, nrow = K_new, ncol = K_old)
  for (i in seq_len(K_new)) dist_matrix[i, ] <- colSums(W_new[, i] != W_old)
  min_new_to_old <- apply(dist_matrix, 1, min)
  min_old_to_new <- apply(dist_matrix, 2, min)
  all_distances <- c(min_new_to_old, min_old_to_new)
  list(overlap_pct = sum(min_new_to_old <= 1L) / K_new,
       avg_distance = mean(all_distances),
       max_distance = max(all_distances))
}

snap_candidate_to_basis <- function(w_init, U, mu, M, h_diag, max_iter = 10) {
  m <- length(w_init); w_curr <- as.integer(w_init)
  for (it in seq_len(max_iter)) {
    w_cen <- w_curr - mu; r <- w_cen - U %*% crossprod(U, w_cen)
    delta_1bit <- 2 * (1L - 2L * w_curr) * as.vector(r) + h_diag
    min_1bit_idx <- which.min(delta_1bit)
    if (delta_1bit[min_1bit_idx] < -1e-5) {
      w_curr[min_1bit_idx] <- 1L - w_curr[min_1bit_idx]; next
    }
    delta_sign <- (1L - 2L * w_curr)
    d2_mat <- outer(delta_1bit, delta_1bit, "+") + 2 * outer(delta_sign, delta_sign) * M
    diag(d2_mat) <- Inf
    if (min(d2_mat) < -1e-5) {
      idx <- which(d2_mat == min(d2_mat), arr.ind = TRUE)[1, ]
      w_curr[idx[1]] <- 1L - w_curr[idx[1]]; w_curr[idx[2]] <- 1L - w_curr[idx[2]]; next
    }
    break
  }
  w_curr
}

snap_pool_to_basis <- function(W_pool, U, mu) {
  m <- nrow(W_pool); P_cands <- ncol(W_pool)
  M <- diag(1.0, m) - tcrossprod(U); h_diag <- diag(M)
  snapped_list <- vector("list", P_cands)
  for (j in seq_len(P_cands)) snapped_list[[j]] <- snap_candidate_to_basis(W_pool[, j], U, mu, M, h_diag)
  W_snapped <- do.call(cbind, snapped_list)
  W_enriched <- cbind(W_pool, W_snapped)
  W_enriched[, !duplicated(t(W_enriched)), drop = FALSE]
}

snap_pool_only <- function(W_pool, U, mu) {
  m <- nrow(W_pool); P_cands <- ncol(W_pool)
  M <- diag(1.0, m) - tcrossprod(U); h_diag <- diag(M)
  snapped_list <- vector("list", P_cands)
  for (j in seq_len(P_cands)) snapped_list[[j]] <- snap_candidate_to_basis(W_pool[, j], U, mu, M, h_diag)
  W_snapped <- do.call(cbind, snapped_list)
  W_snapped[, !duplicated(t(W_snapped)), drop = FALSE]
}

exchange_from_init <- function(P_obs, W_pool, init_indices, max_cycles = 2) {
  P_cands <- ncol(W_pool); K <- length(init_indices); selected <- init_indices
  best_mse <- fit_fcls_qp(P_obs, W_pool[, selected, drop = FALSE])$mse
  for (cycle in seq_len(max_cycles)) {
    swaps <- 0L
    for (i in seq_len(K)) {
      if (best_mse < 1e-8) break
      curr <- selected[i]; best_cand <- curr
      remaining <- setdiff(seq_len(P_cands), selected)
      for (cand in remaining) {
        test_idx <- selected; test_idx[i] <- cand
        fit_test <- fit_fcls_qp(P_obs, W_pool[, test_idx, drop = FALSE])
        if (fit_test$mse < (best_mse - 1e-6)) {
          best_mse <- fit_test$mse; best_cand <- cand
        }
      }
      if (best_cand != curr) { selected[i] <- best_cand; swaps <- swaps + 1L }
    }
    if (swaps == 0L || best_mse < 1e-8) break
  }
  list(indices = selected, mse = best_mse)
}

# Bulletproof multi-start selector with dynamic column guards
select_best_k_subset_multistart <- function(P_obs, W_pool, hunter_sols, K_target, top_starts = 2) {
  P_cands <- ncol(W_pool)
  if (P_cands <= K_target) return(fit_fcls_qp(P_obs, W_pool))
  all_cands_idx <- seq_len(P_cands)
  init_sets <- list()

  for (s in seq_along(hunter_sols)) {
    W_init <- hunter_sols[[s]]
    if (is.null(W_init) || ncol(W_init) == 0L) next
    k_init <- ncol(W_init)
    init_idx <- vapply(seq_len(k_init), function(k) which.min(colSums(abs(W_pool - W_init[, k]))), integer(1))
    init_idx <- unique(init_idx)
    if (length(init_idx) < K_target) {
      rem <- setdiff(all_cands_idx, init_idx)
      init_idx <- c(init_idx, head(rem, K_target - length(init_idx)))
    }
    init_idx <- init_idx[seq_len(K_target)]
    init_key <- paste(sort(init_idx), collapse = "_")
    if (!init_key %in% names(init_sets)) init_sets[[init_key]] <- init_idx
  }
  if (length(init_sets) == 0) init_sets[["default"]] <- seq_len(K_target)

  init_scores <- vapply(init_sets, function(idx) fit_fcls_qp(P_obs, W_pool[, idx, drop = FALSE])$mse, numeric(1))
  ranked_order <- order(init_scores)
  selected_starts <- head(ranked_order, min(top_starts, length(init_sets)))
  best_overall_mse <- Inf; best_overall_indices <- NULL
  for (s_idx in selected_starts) {
    fit_run <- exchange_from_init(P_obs, W_pool, init_sets[[s_idx]], max_cycles = 2)
    if (fit_run$mse < best_overall_mse) {
      best_overall_mse <- fit_run$mse; best_overall_indices <- fit_run$indices
    }
  }
  fit_fcls_qp(P_obs, W_pool[, best_overall_indices, drop = FALSE])
}

polish_bits_v7 <- function(X, P_obs, W_init, U_basis, mu_vec, eps = 1e-3, max_epochs = 3) {
  W <- W_init; M <- nrow(W); K <- ncol(W); n <- ncol(P_obs)
  h_diag <- 1.0 - rowSums(U_basis^2)
  inv_0 <- which(rowSums(X) == 0); inv_1 <- which(rowSums(X) == n * 100)
  W[inv_0, ] <- 0L; W[inv_1, ] <- 1L

  H <- fit_fcls_h(P_obs, W); P_hat <- W %*% H
  res_mat <- (W - mu_vec) - U_basis %*% crossprod(U_basis, W - mu_vec)
  W_best <- W; H_best <- H
  ll_best <- sum(dbinom(X, 100, pmax(eps, pmin(1 - eps, P_hat)), log = TRUE))
  mse_best <- mean((P_obs - P_hat)^2)

  for (epoch in seq_len(max_epochs)) {
    flips <- 0L; mean_h <- rowMeans(H)
    for (m in seq_len(M)) {
      if (m %in% inv_0 || m %in% inv_1) next
      x_m <- X[m, ]; row_p <- P_hat[m, ]
      ll_base <- sum(dbinom(x_m, 100, pmax(eps, pmin(1 - eps, row_p)), log = TRUE))
      best_d1 <- 0; best_k1 <- 0
      for (k in seq_len(K)) {
        if (mean_h[k] < 0.01) next
        d_bit <- 1L - 2L * W[m, k]; p_tr <- row_p + d_bit * H[k, ]
        dll <- sum(dbinom(x_m, 100, pmax(eps, pmin(1 - eps, p_tr)), log = TRUE)) - ll_base
        tau_k <- max(2.5, min(5.0, 4.5 * sqrt(mean_h[k] / 0.18)))
        if (dll > tau_k && dll > best_d1) {
          d_res <- 2 * d_bit * res_mat[m, k] + h_diag[m]
          if (d_res <= 0.85) { best_d1 <- dll; best_k1 <- k }
        }
      }
      if (best_k1 > 0) {
        d_bit <- 1L - 2L * W[m, best_k1]; W[m, best_k1] <- 1L - W[m, best_k1]
        P_hat[m, ] <- row_p + d_bit * H[best_k1, ]
        res_mat[, best_k1] <- (W[, best_k1] - mu_vec) - U_basis %*% crossprod(U_basis, W[, best_k1] - mu_vec)
        flips <- flips + 1L
      }
    }
    if (flips == 0L) break
    H <- fit_fcls_h(P_obs, W); P_hat <- W %*% H
    curr_ll <- sum(dbinom(X, 100, pmax(eps, pmin(1 - eps, P_hat)), log = TRUE))
    curr_mse <- mean((P_obs - P_hat)^2)
    if (curr_ll > ll_best && curr_mse <= (mse_best * 1.01)) {
      W_best <- W; H_best <- H; ll_best <- curr_ll; mse_best <- curr_mse
    } else { break }
  }
  list(W = W_best, H = H_best, mse = mse_best)
}

polish_strains <- function(X, P_obs, W_init, U_basis, mu_vec, drift_ceiling = 0.85) {
  W <- W_init; m <- nrow(W); K <- ncol(W); n <- ncol(P_obs)
  h_diag <- 1.0 - rowSums(U_basis^2)
  inv_0 <- which(rowSums(X) == 0); inv_1 <- which(rowSums(X) == n * 100)
  W[inv_0, ] <- 0L; W[inv_1, ] <- 1L
  H <- fit_fcls_h(P_obs, W); P_hat <- W %*% H
  res_mat <- (W - mu_vec) - U_basis %*% crossprod(U_basis, W - mu_vec)

  for (epoch in 1:2) {
    flips <- 0L; mean_h <- rowMeans(H)
    for (loc in seq_len(m)) {
      if (loc %in% inv_0 || loc %in% inv_1) next
      x_loc <- X[loc, ]; row_p <- P_hat[loc, ]
      ll_base <- sum(dbinom(x_loc, 100, pmax(1e-4, pmin(1 - 1e-4, row_p)), log = TRUE))
      best_dll <- 0; best_k <- 0
      for (k in seq_len(K)) {
        if (mean_h[k] < 0.01) next
        d_bit <- 1L - 2L * W[loc, k]; p_trial <- row_p + d_bit * H[k, ]
        dll <- sum(dbinom(x_loc, 100, pmax(1e-4, pmin(1 - 1e-4, p_trial)), log = TRUE)) - ll_base
        tau_k <- max(2.5, min(5.0, 4.5 * sqrt(mean_h[k] / 0.18)))
        if (!is.na(dll) && dll > tau_k && dll > best_dll) {
          d_res <- 2 * d_bit * res_mat[loc, k] + h_diag[loc]
          if (d_res <= drift_ceiling) { best_dll <- dll; best_k <- k }
        }
      }
      if (best_k > 0) {
        d_bit <- 1L - 2L * W[loc, best_k]; W[loc, best_k] <- 1L - W[loc, best_k]
        P_hat[loc, ] <- row_p + d_bit * H[best_k, ]
        res_mat[, best_k] <- (W[, best_k] - mu_vec) - U_basis %*% crossprod(U_basis, W[, best_k] - mu_vec)
        flips <- flips + 1L
      }
    }
    if (flips == 0L) break
    H <- fit_fcls_h(P_obs, W); P_hat <- W %*% H
  }
  W
}

extract_slawski <- function(U, mu, K, top_n = 25) {
  m <- nrow(U); d <- K - 1L
  pivots <- qr(t(U), tol = 1e-7)$pivot[seq_len(d)]
  U_R <- U[pivots, , drop = FALSE]; U_S <- U[setdiff(seq_len(m), pivots), , drop = FALSE]
  inv_R <- tryCatch(solve(crossprod(U_R) + 0.02 * norm(U_R, "F")^2 * diag(d)) %*% t(U_R),
                    error = function(e) MASS::ginv(U_R))
  G_hat <- U_S %*% inv_R; h_hat <- mu[setdiff(seq_len(m), pivots)] - G_hat %*% mu[pivots]
  b_grid <- t(as.matrix(expand.grid(rep(list(c(0L, 1L)), d))))
  z_S <- sweep(G_hat %*% b_grid, 1, as.vector(h_hat), "+")
  W_slaw <- matrix(0L, nrow = m, ncol = ncol(b_grid))
  W_slaw[pivots, ] <- b_grid; W_slaw[setdiff(seq_len(m), pivots), ] <- (z_S >= 0.5) * 1L
  W_unique <- W_slaw[, !duplicated(t(W_slaw)), drop = FALSE]
  cen <- sweep(W_unique, 1, mu, "-")
  res <- colSums((cen - U %*% crossprod(U, cen))^2)
  W_unique[, order(res)[seq_len(min(top_n, ncol(W_unique)))], drop = FALSE]
}

extract_euclidean_rays <- function(U, mu, n_rays = 40, top_n = 25) {
  m <- nrow(U); d <- ncol(U)
  D <- matrix(rnorm(d * n_rays), nrow = d, ncol = n_rays)
  D <- sweep(D, 2, sqrt(colSums(D^2)), "/")
  V <- U %*% D
  raw <- cbind((mu + (sqrt(m)/2) * V) >= 0.5, (mu - (sqrt(m)/2) * V) >= 0.5) * 1L
  raw_unique <- raw[, !duplicated(t(raw)), drop = FALSE]
  cen <- sweep(raw_unique, 1, mu, "-")
  res <- colSums((cen - U %*% crossprod(U, cen))^2)
  raw_unique[, order(res)[seq_len(min(top_n, ncol(raw_unique)))], drop = FALSE]
}

extract_arcsine_rays <- function(P_obs, U_e, mu_e, K, n_rays = 50, top_n = 25) {
  m <- nrow(P_obs); d <- K - 1L
  P_clamped <- pmin(pmax(P_obs, 1e-6), 1 - 1e-6)
  Theta <- 2 * asin(sqrt(P_clamped)); mu_t <- rowMeans(Theta)
  sv_t <- svd(sweep(Theta, 1, mu_t, "-"), nu = d, nv = 0)
  U_t <- sv_t$u[, seq_len(d), drop = FALSE]
  D <- matrix(rnorm(d * n_rays), nrow = d, ncol = n_rays)
  D <- sweep(D, 2, sqrt(colSums(D^2)), "/")
  V <- U_t %*% D
  rays_t <- cbind(mu_t + (sqrt(m)/2) * V, mu_t - (sqrt(m)/2) * V)
  rays_p <- sin(pmin(pmax(rays_t, 0), pi) / 2)^2
  raw_u <- (rays_p >= 0.5) * 1L
  raw_u <- raw_u[, !duplicated(t(raw_u)), drop = FALSE]
  cen <- sweep(raw_u, 1, mu_e, "-")
  res <- colSums((cen - U_e %*% crossprod(U_e, cen))^2)
  raw_u[, order(res)[seq_len(min(top_n, ncol(raw_u)))], drop = FALSE]
}

harvest_dual_20svd <- function(P_obs, K) {
  m <- nrow(P_obs); n <- ncol(P_obs); d <- K - 1L; mu <- rowMeans(P_obs)
  P_cent <- sweep(P_obs, 1, mu, "-")
  pool_list <- list()
  sv_b <- svd(P_cent, nu = d, nv = 0); U_b <- sv_b$u[, seq_len(d), drop = FALSE]
  pool_list[[1]] <- cbind(
    extract_slawski(U_b, mu, K, top_n = 30),
    extract_euclidean_rays(U_b, mu, n_rays = 40, top_n = 25),
    extract_arcsine_rays(P_obs, U_b, mu, K, n_rays = 50, top_n = 25)
  )
  for (e in 2:15) {
    w_loc <- rexp(m, 1.0); sqrt_w <- sqrt(w_loc / mean(w_loc))
    sv_w <- svd(sweep(P_cent, 1, sqrt_w, "*"), nu = d, nv = 0)
    U_w <- qr.Q(qr(sweep(sv_w$u[, seq_len(d), drop = FALSE], 1, sqrt_w, "/")))[, seq_len(d), drop = FALSE]
    pool_list[[e]] <- cbind(extract_slawski(U_w, mu, K, top_n = 20),
                            extract_euclidean_rays(U_w, mu, n_rays = 30, top_n = 20))
  }
  for (e in 16:20) {
    w_s <- rexp(n, 1.0); w_s <- (0.8 * (w_s / sum(w_s))) + (0.2 / n); w_s <- w_s / sum(w_s)
    mu_s <- as.vector(P_obs %*% w_s)
    sv_s <- svd(sweep(sweep(P_obs, 1, mu_s, "-"), 2, sqrt(w_s), "*"), nu = d, nv = 0)
    U_s <- sv_s$u[, seq_len(d), drop = FALSE]
    pool_list[[e]] <- cbind(extract_slawski(U_s, mu_s, K, top_n = 20),
                            extract_euclidean_rays(U_s, mu_s, n_rays = 30, top_n = 20))
  }
  all_c <- do.call(cbind, pool_list)
  all_c[, !duplicated(t(all_c)), drop = FALSE]
}

coord_exchange_robust_parsimony <- function(P_obs, W_pool, K, lambda_pars = 5e-5) {
  m <- nrow(W_pool); d <- K - 1L
  A_pool <- rbind(W_pool, rep(1.0, ncol(W_pool)))
  all_cands_idx <- seq_len(ncol(W_pool))
  piv_base <- qr(A_pool, tol = 1e-6)$pivot[seq_len(K)]
  cand_muts <- colSums(W_pool)

  eval_loss <- function(idx) {
    fit <- fit_fcls_qp(P_obs, W_pool[, idx, drop = FALSE])
    fit$mse + lambda_pars * (sum(cand_muts[idx]) / (m * K))
  }
  best_overall_loss <- Inf; best_overall_set <- piv_base

  for (restart in 1:3) {
    if (restart == 1) {
      selected <- piv_base
    } else if (restart == 2) {
      mu <- rowMeans(P_obs); sv <- svd(sweep(P_obs, 1, mu, "-"), nu = d, nv = 0)
      U <- sv$u[, seq_len(d), drop = FALSE]; cen <- sweep(W_pool, 1, mu, "-")
      res_norms <- colSums((cen - U %*% crossprod(U, cen))^2)
      top_res <- order(res_norms)[seq_len(min(40, ncol(W_pool)))]
      cand_pool_2 <- setdiff(top_res, piv_base)
      num_swap <- min(3, length(cand_pool_2))
      selected <- if (num_swap > 0) c(piv_base[seq_len(K - num_swap)], sample(cand_pool_2, num_swap)) else piv_base
    } else {
      sparse_pool <- order(cand_muts)[seq_len(min(60, ncol(W_pool)))]
      A_sparse <- rbind(W_pool[, sparse_pool, drop = FALSE], 1.0)
      qr_sparse <- qr(A_sparse, tol = 1e-6)$pivot
      selected <- sparse_pool[qr_sparse[seq_len(K)]]
    }

    curr_loss <- eval_loss(selected)
    for (cycle in 1:2) {
      swaps <- 0L
      for (i in seq_len(K)) {
        if (curr_loss < 1e-8) break
        curr <- selected[i]; best_c <- curr
        for (cand in setdiff(all_cands_idx, selected)) {
          t_idx <- selected; t_idx[i] <- cand
          test_loss <- eval_loss(t_idx)
          if (test_loss < (curr_loss - 1e-7)) { curr_loss <- test_loss; best_c <- cand }
        }
        if (best_c != curr) { selected[i] <- best_c; swaps <- swaps + 1L }
      }
      if (swaps == 0L || curr_loss < 1e-8) break
    }
    if (curr_loss < best_overall_loss) { best_overall_loss <- curr_loss; best_overall_set <- selected }
  }
  W_pool[, best_overall_set, drop = FALSE]
}

extract_slawski_maxvol <- function(U, mu, K, top_n = 25, is_clean = FALSE) {
  m <- nrow(U); d <- K - 1L
  R_idx <- select_pivots_maxvol(U, K); S_idx <- setdiff(seq_len(m), R_idx)
  U_R <- U[R_idx, , drop = FALSE]; U_S <- U[S_idx, , drop = FALSE]
  damp <- if (is_clean) 0.0 else 0.02 * norm(U_R, "F")^2
  inv_R <- tryCatch(solve(crossprod(U_R) + damp * diag(d)) %*% t(U_R),
                    error = function(e) MASS::ginv(U_R))
  G_hat <- U_S %*% inv_R; h_hat <- mu[S_idx] - G_hat %*% mu[R_idx]
  b_grid <- t(as.matrix(expand.grid(rep(list(c(0L, 1L)), d))))
  z_S <- sweep(G_hat %*% b_grid, 1, as.vector(h_hat), "+")
  W_slaw <- matrix(0L, nrow = m, ncol = ncol(b_grid))
  W_slaw[R_idx, ] <- b_grid; W_slaw[S_idx, ] <- (z_S >= 0.5) * 1L
  W_unique <- W_slaw[, !duplicated(t(W_slaw)), drop = FALSE]
  cen <- sweep(W_unique, 1, mu, "-")
  res <- colSums((cen - U %*% crossprod(U, cen))^2)
  W_unique[, order(res)[seq_len(min(top_n, ncol(W_unique)))], drop = FALSE]
}

harvest_v4_dual <- function(P_obs, K, is_clean = FALSE) {
  m <- nrow(P_obs); n <- ncol(P_obs); d <- K - 1L; mu <- rowMeans(P_obs)
  P_cent <- sweep(P_obs, 1, mu, "-")
  pool_list <- list()
  slaw_k_starts <- list()

  sv_b <- svd(P_cent, nu = d, nv = 0); U_b <- sv_b$u[, seq_len(d), drop = FALSE]
  W_slaw_b <- extract_slawski_maxvol(U_b, mu, K, top_n = 30, is_clean = is_clean)
  slaw_k_starts[[1]] <- W_slaw_b[, seq_len(min(K, ncol(W_slaw_b))), drop = FALSE]

  pool_list[[1]] <- cbind(
    W_slaw_b,
    extract_euclidean_rays(U_b, mu, n_rays = 40, top_n = 25),
    extract_arcsine_rays(P_obs, U_b, mu, K, n_rays = 50, top_n = 25)
  )

  for (e in 2:15) {
    w_loc <- rexp(m, 1.0); sqrt_w <- sqrt(w_loc / mean(w_loc))
    sv_w <- svd(sweep(P_cent, 1, sqrt_w, "*"), nu = d, nv = 0)
    U_w <- qr.Q(qr(sweep(sv_w$u[, seq_len(d), drop = FALSE], 1, sqrt_w, "/")))[, seq_len(d), drop = FALSE]
    W_slaw_w <- extract_slawski_maxvol(U_w, mu, K, top_n = 20, is_clean = is_clean)
    if (ncol(W_slaw_w) >= K) slaw_k_starts[[length(slaw_k_starts) + 1L]] <- W_slaw_w[, seq_len(K), drop = FALSE]
    pool_list[[e]] <- cbind(W_slaw_w, extract_euclidean_rays(U_w, mu, n_rays = 30, top_n = 20))
  }

  for (e in 16:20) {
    w_s <- rexp(n, 1.0); w_s <- (0.8 * (w_s / sum(w_s))) + (0.2 / n); w_s <- w_s / sum(w_s)
    mu_s <- as.vector(P_obs %*% w_s)
    sv_s <- svd(sweep(sweep(P_obs, 1, mu_s, "-"), 2, sqrt(w_s), "*"), nu = d, nv = 0)
    U_s <- sv_s$u[, seq_len(d), drop = FALSE]
    pool_list[[e]] <- cbind(
      extract_slawski_maxvol(U_s, mu_s, K, top_n = 20, is_clean = is_clean),
      extract_euclidean_rays(U_s, mu_s, n_rays = 30, top_n = 20)
    )
  }

  all_c <- do.call(cbind, pool_list)
  list(W_pool = all_c[, !duplicated(t(all_c)), drop = FALSE], slaw_starts = slaw_k_starts)
}

# 100% Guarded 4-Start Coordinate Exchange (zero subscript out of bounds risk)
coord_exchange_v4 <- function(P_obs, W_pool, slaw_starts, U_b, mu_b, K, use_subspace_reg = FALSE) {
  m <- nrow(W_pool); P_cands <- ncol(W_pool)
  all_cands_idx <- seq_len(P_cands)
  cand_muts <- colSums(W_pool)

  cen <- sweep(W_pool, 1, mu_b, "-")
  cand_res_sq <- colSums((cen - U_b %*% crossprod(U_b, cen))^2)

  P_target <- if (use_subspace_reg) {
    P_svd <- sweep(U_b %*% crossprod(U_b, sweep(P_obs, 1, mu_b, "-")), 1, mu_b, "+")
    pmin(pmax(P_svd, 0.0), 1.0)
  } else {
    P_obs
  }

  eval_loss <- function(idx) {
    fit <- fit_fcls_qp(P_target, W_pool[, idx, drop = FALSE])
    if (use_subspace_reg) {
      fit$mse + 0.03 * (sum(cand_res_sq[idx]) / (m * K))
    } else {
      fit$mse
    }
  }

  map_to_pool <- function(W_init) {
    if (is.null(W_init) || ncol(W_init) == 0L) return(seq_len(K))
    k_init <- ncol(W_init)
    idx <- vapply(seq_len(k_init), function(k) which.min(colSums(abs(W_pool - W_init[, k]))), integer(1))
    idx <- unique(idx)
    if (length(idx) < K) {
      rem <- setdiff(all_cands_idx, idx)
      idx <- c(idx, head(rem, K - length(idx)))
    }
    idx[seq_len(K)]
  }

  slaw_indices <- lapply(slaw_starts, map_to_pool)
  if (length(slaw_indices) == 0L) slaw_indices <- list(seq_len(K))
  slaw_scores  <- vapply(slaw_indices, function(idx) fit_fcls_qp(P_obs, W_pool[, idx, drop = FALSE])$mse, numeric(1))
  best_slaw_order <- order(slaw_scores)

  piv_qr <- qr(rbind(W_pool, 1.0), tol = 1e-6)$pivot
  piv_qr <- unique(piv_qr[piv_qr <= P_cands])
  if (length(piv_qr) < K) piv_qr <- c(piv_qr, head(setdiff(all_cands_idx, piv_qr), K - length(piv_qr)))
  piv_qr <- piv_qr[seq_len(K)]

  sp_pool <- order(cand_muts)[seq_len(min(60, P_cands))]
  piv_sp  <- sp_pool[qr(rbind(W_pool[, sp_pool, drop = FALSE], 1.0), tol = 1e-6)$pivot]
  piv_sp  <- unique(piv_sp[!is.na(piv_sp) & piv_sp <= P_cands])
  if (length(piv_sp) < K) piv_sp <- c(piv_sp, head(setdiff(all_cands_idx, piv_sp), K - length(piv_sp)))
  piv_sp  <- piv_sp[seq_len(K)]

  starts_list <- list(
    slaw_indices[[1]],
    slaw_indices[[best_slaw_order[1]]],
    piv_qr,
    piv_sp
  )
  starts_keys <- vapply(starts_list, function(idx) paste(sort(idx), collapse = "_"), character(1))
  starts_list <- starts_list[!duplicated(starts_keys)]

  best_overall_loss <- Inf; best_overall_set <- starts_list[[1]]
  for (selected in starts_list) {
    curr_loss <- eval_loss(selected)
    for (cycle in 1:2) {
      swaps <- 0L
      for (i in seq_len(K)) {
        if (curr_loss < 1e-8) break
        curr <- selected[i]; best_c <- curr
        for (cand in setdiff(all_cands_idx, selected)) {
          t_idx <- selected; t_idx[i] <- cand
          test_loss <- eval_loss(t_idx)
          if (test_loss < (curr_loss - 1e-7)) { curr_loss <- test_loss; best_c <- cand }
        }
        if (best_c != curr) { selected[i] <- best_c; swaps <- swaps + 1L }
      }
      if (swaps == 0L || curr_loss < 1e-8) break
    }
    true_mse <- fit_fcls_qp(P_obs, W_pool[, selected, drop = FALSE])$mse
    if (true_mse < best_overall_loss) {
      best_overall_loss <- true_mse; best_overall_set <- selected
    }
  }
  W_pool[, best_overall_set, drop = FALSE]
}

polish_bits_v7_full <- function(X, P_obs, W_init, U_basis, mu_vec, eps = 1e-3, max_epochs = 3) {
  W <- W_init; M <- nrow(W); K <- ncol(W); n <- ncol(P_obs)
  h_diag <- 1.0 - rowSums(U_basis^2)
  inv_0 <- which(rowSums(X) == 0); inv_1 <- which(rowSums(X) == n * 100)
  W[inv_0, ] <- 0L; W[inv_1, ] <- 1L

  H <- fit_fcls_h(P_obs, W); P_hat <- W %*% H
  res_mat <- (W - mu_vec) - U_basis %*% crossprod(U_basis, W - mu_vec)
  W_best <- W; H_best <- H
  ll_best <- sum(dbinom(X, 100, pmax(eps, pmin(1 - eps, P_hat)), log = TRUE))
  mse_best <- mean((P_obs - P_hat)^2)
  pairs <- if (K > 1) combn(K, 2, simplify = FALSE) else list()

  row_dev_fn <- function(x_r, p_r) {
    pc <- pmax(eps, pmin(1 - eps, p_r)); mv <- 100 * pc; d <- numeric(length(x_r))
    px <- x_r > 0; d[px] <- d[px] + 2 * x_r[px] * log(x_r[px] / mv[px])
    pn <- (100 - x_r) > 0; d[pn] <- d[pn] + 2 * (100 - x_r[pn]) * log((100 - x_r[pn]) / (100 - mv[pn]))
    sum(d)
  }

  for (epoch in seq_len(max_epochs)) {
    flips <- 0L; mean_h <- rowMeans(H)
    for (m in seq_len(M)) {
      if (m %in% inv_0 || m %in% inv_1) next
      x_m <- X[m, ]; row_p <- P_hat[m, ]
      ll_base <- sum(dbinom(x_m, 100, pmax(eps, pmin(1 - eps, row_p)), log = TRUE))
      best_d1 <- 0; best_k1 <- 0
      for (k in seq_len(K)) {
        if (mean_h[k] < 0.01) next
        d_bit <- 1L - 2L * W[m, k]; p_tr <- row_p + d_bit * H[k, ]
        dll <- sum(dbinom(x_m, 100, pmax(eps, pmin(1 - eps, p_tr)), log = TRUE)) - ll_base
        tau_k <- max(2.5, min(5.0, 4.5 * sqrt(mean_h[k] / 0.18)))
        if (dll > tau_k && dll > best_d1) {
          if (2 * d_bit * res_mat[m, k] + h_diag[m] <= 0.85) { best_d1 <- dll; best_k1 <- k }
        }
      }
      if (best_k1 > 0) {
        d_bit <- 1L - 2L * W[m, best_k1]; W[m, best_k1] <- 1L - W[m, best_k1]
        P_hat[m, ] <- row_p + d_bit * H[best_k1, ]
        res_mat[, best_k1] <- (W[, best_k1] - mu_vec) - U_basis %*% crossprod(U_basis, W[, best_k1] - mu_vec)
        flips <- flips + 1L
      }
    }
    for (m in seq_len(M)) {
      if (m %in% inv_0 || m %in% inv_1) next
      x_m <- X[m, ]; row_p <- P_hat[m, ]
      act <- which(W[m, ] == 1L)
      max_ha <- if (length(act) > 0) max(mean_h[act]) else 0.0
      if (row_dev_fn(x_m, row_p) <= n + sqrt(2 * n) * max(0.4, 1.0 - 0.5 * max_ha)) next
      ll_base <- sum(dbinom(x_m, 100, pmax(eps, pmin(1 - eps, row_p)), log = TRUE))
      best_d2 <- 0; best_pr <- NULL
      for (pr in pairs) {
        k1 <- pr[1]; k2 <- pr[2]; if (mean_h[k1] < 0.01 && mean_h[k2] < 0.01) next
        ph <- mean_h[k1] + mean_h[k2]
        if (W[m, k1] != W[m, k2]) {
          dp <- (W[m, k2] - W[m, k1]) * H[k1, ] + (W[m, k1] - W[m, k2]) * H[k2, ]
          dll <- sum(dbinom(x_m, 100, pmax(eps, pmin(1 - eps, row_p + dp)), log = TRUE)) - ll_base
          if (dll > max(4.0, min(5.2, 4.2 * sqrt(ph / 0.36))) && dll > best_d2) {
            db1 <- 1L - 2L * W[m, k1]; db2 <- 1L - 2L * W[m, k2]
            if ((2 * db1 * res_mat[m, k1] + h_diag[m] <= 0.85) && (2 * db2 * res_mat[m, k2] + h_diag[m] <= 0.85)) {
              best_d2 <- dll; best_pr <- pr
            }
          }
        } else {
          s <- 1L - 2L * W[m, k1]; dp <- s * (H[k1, ] + H[k2, ])
          dll <- sum(dbinom(x_m, 100, pmax(eps, pmin(1 - eps, row_p + dp)), log = TRUE)) - ll_base
          if (dll > max(6.0, min(8.0, 6.8 * sqrt(ph / 0.36))) && dll > best_d2) {
            if ((2 * s * res_mat[m, k1] + h_diag[m] <= 0.85) && (2 * s * res_mat[m, k2] + h_diag[m] <= 0.85)) {
              best_d2 <- dll; best_pr <- pr
            }
          }
        }
      }
      if (!is.null(best_pr)) {
        k1 <- best_pr[1]; k2 <- best_pr[2]
        W[m, k1] <- 1L - W[m, k1]; W[m, k2] <- 1L - W[m, k2]
        P_hat[m, ] <- W[m, , drop = FALSE] %*% H
        res_mat[, k1] <- (W[, k1] - mu_vec) - U_basis %*% crossprod(U_basis, W[, k1] - mu_vec)
        res_mat[, k2] <- (W[, k2] - mu_vec) - U_basis %*% crossprod(U_basis, W[, k2] - mu_vec)
        flips <- flips + 2L
      }
    }
    if (flips == 0L) break
    H <- fit_fcls_h(P_obs, W); P_hat <- W %*% H
    curr_ll <- sum(dbinom(X, 100, pmax(eps, pmin(1 - eps, P_hat)), log = TRUE))
    curr_mse <- mean((P_obs - P_hat)^2)
    if (curr_ll > ll_best && curr_mse <= (mse_best * 1.01)) {
      W_best <- W; H_best <- H; ll_best <- curr_ll; mse_best <- curr_mse
    } else { break }
  }
  W_best
}

# ------------------------------------------------------------------------------
# 2. GENERATE FUKUYAMA TASK DATA & PRUNED CATALOG
# ------------------------------------------------------------------------------

make_fukuyama <- function(K, seed = 202601L) {
  set.seed(seed)
  d <- K - 1L; m <- 365L; n <- 20L
  p <- 350L %/% d; s <- m - d * p
  B_p <- matrix(0L, d * p, d); for (i in 1:d) B_p[((i - 1L) * p + 1L):(i * p), i] <- 1L
  B_s <- matrix(0L, s, d); B_s[, 1:2] <- 1L
  B   <- rbind(B_p, B_s); W_t <- cbind(0L, B)
  G   <- matrix(runif(K * n, 0.9, 1.1), K, n)
  H_t <- sweep(G, 2, colSums(G), "/"); S_t <- W_t %*% H_t
  X_n <- matrix(rbinom(m * n, 100L, as.vector(S_t)), m, n)

  p_list <- list(); p_map <- matrix(0L, d, d); idx <- d + 1L
  for (a in 1:(d - 1L)) for (b in (a + 1L):d) if (!(a == 1L && b == 2L)) {
    idx <- idx + 1L; p_list[[length(p_list) + 1L]] <- B[, a] + B[, b]
    p_map[a, b] <- idx; p_map[b, a] <- idx
  }
  W_dict <- cbind(W_t, do.call(cbind, p_list))
  list(W_true = W_t, B = B, H_true = H_t, S_true = S_t, X = X_n, P = X_n / 100.0,
       W_dict = W_dict, pair_map = p_map, K = K, d = d, p = p, s = s)
}

build_pruned_catalog <- function(K, pair_map) {
  d <- K - 1L
  adj <- matrix(1L, d, d); diag(adj) <- 0L; adj[1, 2] <- 0L; adj[2, 1] <- 0L
  elem_c <- vector("list", 27000L); e_cnt <- 0L
  for (sv in 1:d) {
    vis <- logical(d); pth <- integer(d)
    dfs <- function(u, dep) {
      pth[dep] <<- u; vis[u] <<- TRUE
      if (dep >= 3L && (dep %% 2L == 1L) && adj[u, sv] == 1L && pth[2L] < pth[dep]) {
        ns <- pth[1:dep]; ns_n <- c(ns[-1L], ns[1L])
        e_cnt <<- e_cnt + 1L
        elem_c[[e_cnt]] <<- list(nodes = ns, len = dep, mask = sum(bitwShiftL(1L, ns - 1L)),
                                 dict_cols = pair_map[cbind(ns, ns_n)])
      }
      if (dep < d) for (v in which(adj[u, ] == 1L)) if (v > sv && !vis[v]) dfs(v, dep + 1L)
      vis[u] <<- FALSE
    }
    dfs(sv, 1L)
  }
  elem_c  <- elem_c[seq_len(e_cnt)]
  c_lens  <- if (e_cnt > 0L) vapply(elem_c, function(x) x$len, integer(1)) else integer(0)
  c_masks <- if (e_cnt > 0L) vapply(elem_c, function(x) x$mask, integer(1)) else integer(0)
  comb_i  <- which(c_lens <= (d - 3L)); term_i <- which(c_lens > (d - 3L))

  parts <- vector("list", 35000L); p_cnt <- 1L
  parts[[1L]] <- list(cols = seq_len(K), family = "0-cycle", key = paste(seq_len(K), collapse = "_"))
  for (idx in term_i) {
    p_cnt <- p_cnt + 1L
    un_s <- which(bitwAnd(c_masks[idx], bitwShiftL(1L, seq_len(d) - 1L)) == 0L) + 1L
    cls  <- sort(c(1L, un_s, elem_c[[idx]]$dict_cols))
    parts[[p_cnt]] <- list(cols = cls, family = paste0(c_lens[idx], "-cycle"), key = paste(cls, collapse = "_"))
  }
  if (length(comb_i) > 0L) {
    bt <- function(s_p, u_m, c_ids, u_l, c_ls) {
      for (pos in s_p:length(comb_i)) {
        i <- comb_i[pos]
        if (u_l + c_lens[i] <= d && bitwAnd(u_m, c_masks[i]) == 0L) {
          p_cnt <<- p_cnt + 1L
          nxt_ids <- c(c_ids, i); nxt_m <- bitwOr(u_m, c_masks[i])
          nxt_ls  <- sort(c(c_ls, c_lens[i]), decreasing = TRUE)
          un_s    <- which(bitwAnd(nxt_m, bitwShiftL(1L, seq_len(d) - 1L)) == 0L) + 1L
          p_cols  <- unlist(lapply(nxt_ids, function(id) elem_c[[id]]$dict_cols))
          cls     <- sort(c(1L, un_s, p_cols))
          parts[[p_cnt]] <<- list(cols = cls, family = paste0(paste(nxt_ls, collapse = "+"), "-cycle"),
                                  key = paste(cls, collapse = "_"))
          if (u_l + c_lens[i] <= d - 3L && pos < length(comb_i)) bt(pos + 1L, nxt_m, nxt_ids, u_l + c_lens[i], nxt_ls)
        }
      }
    }
    bt(1L, 0L, integer(0), 0L, integer(0))
  }
  parts <- parts[seq_len(p_cnt)]
  key_map <- setNames(seq_len(p_cnt), vapply(parts, `[[`, character(1), "key"))
  list(parts = parts, key_map = key_map)
}

classify_W <- function(W_pred, sim, cat_obj) {
  K <- sim$K
  if (is.null(W_pred) || ncol(W_pred) != K) {
    return(list(status = "RANK_DEF", family = "RANK_DEF", sol_id = 0L, min_dist = NA_real_, dist_gt = NA_real_))
  }
  D_dict <- matrix(0L, K, ncol(sim$W_dict))
  for (i in 1:K) D_dict[i, ] <- colSums(W_pred[, i] != sim$W_dict)
  dist_gt <- mean(hungarian_dists(D_dict[, 1:K, drop = FALSE]))

  if (all(apply(D_dict, 1, min) == 0L)) {
    m_cols <- sort(apply(D_dict, 1, which.min))
    if (length(unique(m_cols)) == K) {
      k_str <- paste(m_cols, collapse = "_")
      if (k_str %in% names(cat_obj$key_map)) {
        sid <- cat_obj$key_map[[k_str]]
        return(list(status = "EXACT_VALID", family = cat_obj$parts[[sid]]$family,
                    sol_id = sid, min_dist = 0.0, dist_gt = dist_gt))
      } else {
        return(list(status = "DICT_CHIMERA", family = "DICT_CHIMERA",
                    sol_id = 0L, min_dist = 0.0, dist_gt = dist_gt))
      }
    }
  }

  n_c <- length(cat_obj$parts)
  eval_i <- if (n_c <= 513L) seq_len(n_c) else unique(c(1L, round(seq(2L, n_c, length.out = 250L))))
  dists <- vapply(eval_i, function(idx) {
    mean(hungarian_dists(D_dict[, cat_obj$parts[[idx]]$cols, drop = FALSE]))
  }, numeric(1))
  bi <- which.min(dists)
  list(status = if (dists[bi] == 0) "EXACT_VALID" else "NOISY_APPROX",
       family = cat_obj$parts[[eval_i[bi]]]$family,
       sol_id = eval_i[bi], min_dist = dists[bi], dist_gt = dist_gt)
}

run_hunter_audit <- function(Mat_P, Mat_X, K, mode = c("clean_v1_unbiased", "noisy_100x"), max_eval = 60) {
  mode <- match.arg(mode)
  m <- nrow(Mat_P); n <- ncol(Mat_P); d <- K - 1L
  min_dist_tol <- max(15.0, 2.5 * K)
  mu <- rowMeans(Mat_P); P_cent <- sweep(Mat_P, 1, mu, "-")
  U_base <- svd(P_cent, nu = d, nv = 0)$u[, 1:d, drop = FALSE]

  extract_fn <- function(U_in, mu_in) {
    piv <- qr(t(U_in), tol = 1e-7)$pivot[1:d]
    inv_pivot <- tryCatch(solve(U_in[piv, , drop = FALSE]), error = function(e) MASS::ginv(U_in[piv, , drop = FALSE]))
    G_h <- U_in[-piv, , drop = FALSE] %*% inv_pivot
    h_h <- mu_in[-piv] - G_h %*% mu_in[piv]
    b_g <- t(as.matrix(expand.grid(rep(list(c(0L, 1L)), d))))
    z_S <- sweep(G_h %*% b_g, 1, as.vector(h_h), "+")
    W_r <- matrix(0L, m, ncol(b_g)); W_r[piv, ] <- b_g; W_r[-piv, ] <- (z_S >= 0.5) * 1L
    W_u <- W_r[, !duplicated(t(W_r)), drop = FALSE]
    cen <- sweep(W_u, 1, mu_in, "-")
    dsq <- colSums((cen - U_in %*% crossprod(U_in, cen))^2)
    dsq_noisy_tie <- ifelse(dsq < 1e-12, runif(length(dsq), 0, 1e-14), dsq)
    W_ret <- W_u[, order(dsq_noisy_tie)[seq_len(min(K, ncol(W_u)))], drop = FALSE]
    if (ncol(W_ret) < K) {
      pad_k <- K - ncol(W_ret)
      W_ret <- cbind(W_ret, matrix(sample(c(0L, 1L), m * pad_k, replace = TRUE), m, pad_k))
    }
    W_ret
  }

  W_b <- extract_fn(U_base, mu)
  fit_b <- fit_fcls_qp(Mat_P, W_b)
  best_mse <- max(fit_b$mse, 1e-12)
  saved <- list(list(id = 1L, W = fit_b$W, mse = fit_b$mse))

  for (b in seq_len(max_eval)) {
    if (mode == "noisy_100x") {
      Phat <- pmin(pmax(Mat_P, 1e-12), 1 - 1e-12)
      X_star <- matrix(rbinom(m * n, 100L, as.vector(Phat)), m, n)
      fit_s <- fit_weighted_affine_subspace(X_star, rep(100, n), K)
      U_s <- fit_s$U; mu_s <- fit_s$mu
    } else {
      w <- rexp(m, 1.0); sq_w <- sqrt(w / mean(w))
      sv_w <- svd(sweep(P_cent, 1, sq_w, "*"), nu = d, nv = 0)
      U_s <- qr.Q(qr(sweep(sv_w$u[, 1:d, drop = FALSE], 1, sq_w, "/")))[, 1:d, drop = FALSE]
      mu_s <- mu
    }
    W_star <- extract_fn(U_s, mu_s)

    is_diff <- TRUE
    for (s in seq_along(saved)) {
      df <- evaluate_set_difference(W_star, saved[[s]]$W)
      if (df$overlap_pct > 0.60 || df$avg_distance < min_dist_tol) { is_diff <- FALSE; break }
    }
    if (!is_diff) next

    fit_sc <- fit_fcls_qp(Mat_P, W_star)
    if (ncol(fit_sc$W) < K) next
    if (fit_sc$mse <= best_mse * 1.25 + 1e-9) {
      if (fit_sc$mse < best_mse) best_mse <- max(fit_sc$mse, 1e-12)
      saved[[length(saved) + 1L]] <- list(id = length(saved) + 1L, W = fit_sc$W, mse = fit_sc$mse)
    }
  }

  W_pool_raw <- do.call(cbind, lapply(saved, `[[`, "W"))
  W_pool <- W_pool_raw[, !duplicated(t(W_pool_raw)), drop = FALSE]
  W_snap <- snap_pool_to_basis(W_pool, U_base, mu)
  hunter_W_list <- lapply(saved, `[[`, "W")
  fit_s3 <- select_best_k_subset_multistart(Mat_P, W_snap, hunter_W_list, K_target = K)
  fit_s4 <- polish_bits_v7(Mat_X, Mat_P, fit_s3$W, U_base, mu)

  list(n_accepted = length(saved), W_s3 = fit_s3$W, W_s4 = fit_s4$W)
}

# ------------------------------------------------------------------------------
# 3. TASK EXECUTION (ALL 5 PIPELINES HEAD-TO-HEAD)
# ------------------------------------------------------------------------------

sim     <- make_fukuyama(K_val, seed_val)
cat_obj <- build_pruned_catalog(K_val, sim$pair_map)

d <- K_val - 1L; m <- 365L; n <- 20L
mu_c <- rowMeans(sim$S_true); U_c <- svd(sweep(sim$S_true, 1, mu_c, "-"), nu = d, nv = 0)$u[, 1:d, drop = FALSE]

mu_n <- rowMeans(sim$P)
sv_n <- svd(sweep(sim$P, 1, mu_n, "-"), nu = K_val, nv = 0)
U_n  <- sv_n$u[, 1:d, drop = FALSE]
snr_gap_noisy <- if (length(sv_n$d) >= K_val && !is.na(sv_n$d[K_val]) && sv_n$d[K_val] > 1e-9) sv_n$d[d] / sv_n$d[K_val] else 999.0
if (is.na(snr_gap_noisy)) snr_gap_noisy <- 1.0

# 1. PipeA
pa_clean <- run_hunter_audit(sim$S_true, round(100 * sim$S_true), K_val, "clean_v1_unbiased", max_eval = 50)
pa_noisy <- run_hunter_audit(sim$P, sim$X, K_val, "noisy_100x", max_eval = 50)

# 2. PipeB
set.seed(seed_val)
W_pool_b_c <- harvest_dual_20svd(sim$S_true, K_val)
W_pb_c_s3  <- coord_exchange_robust_parsimony(sim$S_true, W_pool_b_c, K_val, lambda_pars = 5e-5)
W_pb_c_s4  <- polish_strains(round(100 * sim$S_true), sim$S_true, W_pb_c_s3, U_c, mu_c, drift_ceiling = 0.85)

set.seed(seed_val)
W_pool_b_n <- harvest_dual_20svd(sim$P, K_val)
W_pb_n_s3  <- coord_exchange_robust_parsimony(sim$P, W_pool_b_n, K_val, lambda_pars = 5e-5)
W_pb_n_s4  <- polish_strains(sim$X, sim$P, W_pb_n_s3, U_n, mu_n, drift_ceiling = 0.85)

# 3. v3
W_v3_c_snap <- snap_pool_to_basis(W_pool_b_c, U_c, mu_c)
W_v3_c_s3   <- coord_exchange_robust_parsimony(sim$S_true, W_v3_c_snap, K_val, lambda_pars = 0.0)
W_v3_c_s4   <- polish_bits_v7_full(round(100 * sim$S_true), sim$S_true, W_v3_c_s3, U_c, mu_c)

W_v3_n_snap <- snap_pool_to_basis(W_pool_b_n, U_n, mu_n)
W_v3_n_s3   <- coord_exchange_robust_parsimony(sim$P, W_v3_n_snap, K_val, lambda_pars = 0.0)
W_v3_n_s4   <- if (snr_gap_noisy > 1.12) polish_bits_v7_full(sim$X, sim$P, W_v3_n_s3, U_n, mu_n) else W_v3_n_s3

# 4 & 5. v4_Std and v4_SubReg
harv_c <- harvest_v4_dual(sim$S_true, K_val, is_clean = TRUE)
W_v4_c_pool <- snap_pool_only(harv_c$W_pool, U_c, mu_c)

W_v4_std_c_s3 <- coord_exchange_v4(sim$S_true, W_v4_c_pool, harv_c$slaw_starts, U_c, mu_c, K_val, use_subspace_reg = FALSE)
W_v4_std_c_s4 <- polish_bits_v7_full(round(100 * sim$S_true), sim$S_true, W_v4_std_c_s3, U_c, mu_c)

W_v4_sub_c_s3 <- coord_exchange_v4(sim$S_true, W_v4_c_pool, harv_c$slaw_starts, U_c, mu_c, K_val, use_subspace_reg = TRUE)
W_v4_sub_c_s4 <- polish_bits_v7_full(round(100 * sim$S_true), sim$S_true, W_v4_sub_c_s3, U_c, mu_c)

harv_n <- harvest_v4_dual(sim$P, K_val, is_clean = FALSE)
W_v4_n_pool <- snap_pool_only(harv_n$W_pool, U_n, mu_n)

W_v4_std_n_s3 <- coord_exchange_v4(sim$P, W_v4_n_pool, harv_n$slaw_starts, U_n, mu_n, K_val, use_subspace_reg = FALSE)
W_v4_std_n_s4 <- if (snr_gap_noisy > 1.12) polish_bits_v7_full(sim$X, sim$P, W_v4_std_n_s3, U_n, mu_n) else W_v4_std_n_s3

W_v4_sub_n_s3 <- coord_exchange_v4(sim$P, W_v4_n_pool, harv_n$slaw_starts, U_n, mu_n, K_val, use_subspace_reg = TRUE)
W_v4_sub_n_s4 <- if (snr_gap_noisy > 1.12) polish_bits_v7_full(sim$X, sim$P, W_v4_sub_n_s3, U_n, mu_n) else W_v4_sub_n_s3

# ------------------------------------------------------------------------------
# 4. SCORING & OUTPUT
# ------------------------------------------------------------------------------

score_clean <- function(Wp) {
  cl <- classify_W(Wp, sim, cat_obj)
  fit <- fit_fcls_qp(sim$S_true, Wp)
  list(exact = (fit$mse < 1e-10 && min(fit$H) >= -1e-5),
       family = cl$family, status = cl$status,
       vol = rel_vol(Wp, sim$W_true), dist_gt = cl$dist_gt, min_dist = cl$min_dist)
}

score_noisy <- function(W_s3, W_s4) {
  cl3 <- classify_W(W_s3, sim, cat_obj)
  cl4 <- classify_W(W_s4, sim, cat_obj)
  list(s3_dist_mk = cl3$min_dist, s3_dist_gt = cl3$dist_gt,
       s4_dist_mk = cl4$min_dist, s4_dist_gt = cl4$dist_gt,
       family = cl4$family, polish_delta = cl4$min_dist - cl3$min_dist)
}

res <- list(
  task_id = task_id, K = K_val, seed = seed_val, rep = rep_idx, snr_gap_noisy = snr_gap_noisy,
  c_pa     = score_clean(pa_clean$W_s4),
  c_pb     = score_clean(W_pb_c_s4),
  c_v3     = score_clean(W_v3_c_s4),
  c_v4_std = score_clean(W_v4_std_c_s4),
  c_v4_sub = score_clean(W_v4_sub_c_s4),
  n_pa     = score_noisy(pa_noisy$W_s3, pa_noisy$W_s4),
  n_pb     = score_noisy(W_pb_n_s3, W_pb_n_s4),
  n_v3     = score_noisy(W_v3_n_s3, W_v3_n_s4),
  n_v4_std = score_noisy(W_v4_std_n_s3, W_v4_std_n_s4),
  n_v4_sub = score_noisy(W_v4_sub_n_s3, W_v4_sub_n_s4)
)

out_dir <- "../results/fukuyama_80_all_pipelines"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(res, file.path(out_dir, sprintf("res_task_%02d_K%d.rds", task_id, K_val)))

cat(sprintf("[DONE] Task %2d (K=%2d, Seed=%d) | Clean: PipeA=%s | PipeB=%s | v4_Sub=%s | NoisyDist: PipeB=%.1f -> v4_Sub=%.1f\n",
            task_id, K_val, seed_val, res$c_pa$family, res$c_pb$family, res$c_v4_sub$family,
            res$n_pb$s4_dist_mk, res$n_v4_sub$s4_dist_mk))
