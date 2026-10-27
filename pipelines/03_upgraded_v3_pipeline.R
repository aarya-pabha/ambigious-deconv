#!/usr/bin/env Rscript
# ==============================================================================
# File: audit_both_pipelines_fukuyama.R
# Purpose: Fully self-contained audit of Pipeline A (sol4 -> snap -> polish_v7)
#          and Pipeline B (unified_15plus) on Generalized Fukuyama Manifolds.
# Usage:   Rscript audit_both_pipelines_fukuyama.R <task_id (1-80)>
# ==============================================================================

suppressPackageStartupMessages({
  library(quadprog)
  library(MASS)
})

# ------------------------------------------------------------------------------
# 1. CORE NUMERICAL SOLVERS & UTILITIES (INLINED, ZERO FILE DEPENDENCY)
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
fit_fcls_quadprog <- fit_fcls_qp
fit_fcls_quadprog_mse <- fit_fcls_qp

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

# ------------------------------------------------------------------------------
# 2. PIPELINE A FUNCTIONS (SOL4 HUNTER -> SNAP_SUBSET -> POLISH_V7)
# ------------------------------------------------------------------------------

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

extract_slawski_vertex_set <- function(U, mu, K) {
  m <- length(mu); d <- max(1L, K - 1L)
  R_idx <- select_pivots_maxvol(U, K); S_idx <- setdiff(seq_len(m), R_idx)
  U_R <- U[R_idx, , drop = FALSE]; U_S <- U[S_idx, , drop = FALSE]
  damp_val <- 0.02 * norm(U_R, "F")^2
  inv_damped <- tryCatch(solve(crossprod(U_R) + damp_val * diag(d)) %*% t(U_R),
                         error = function(e) MASS::ginv(U_R))
  G_hat <- U_S %*% inv_damped; h_hat <- mu[S_idx] - G_hat %*% mu[R_idx]
  b_grid <- t(as.matrix(expand.grid(rep(list(c(0L, 1L)), d))))
  z_S_hat <- sweep(G_hat %*% b_grid, 1, as.vector(h_hat), "+")
  W_rec <- matrix(0L, nrow = m, ncol = ncol(b_grid))
  W_rec[R_idx, ] <- b_grid; W_rec[S_idx, ] <- (z_S_hat >= 0.5) * 1L
  W_unique <- W_rec[, !duplicated(t(W_rec)), drop = FALSE]
  W_centered <- sweep(W_unique, 1, mu, "-")
  dist_sq <- colSums((W_centered - U %*% crossprod(U, W_centered))^2)
  W_unique[, order(dist_sq)[seq_len(min(K, ncol(W_unique)))], drop = FALSE]
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

exchange_from_init <- function(P_obs, W_pool, init_indices, max_cycles = 2) {
  P_cands <- ncol(W_pool); K <- length(init_indices); selected <- init_indices
  best_mse <- fit_fcls_qp(P_obs, W_pool[, selected, drop = FALSE])$mse
  for (cycle in seq_len(max_cycles)) {
    swaps <- 0L
    for (i in seq_len(K)) {
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
    if (swaps == 0L) break
  }
  list(indices = selected, mse = best_mse)
}

select_best_k_subset_multistart <- function(P_obs, W_pool, hunter_sols, K_target, top_starts = 2) {
  P_cands <- ncol(W_pool)
  if (P_cands <= K_target) return(fit_fcls_qp(P_obs, W_pool))
  init_sets <- list()
  for (s in seq_along(hunter_sols)) {
    W_init <- hunter_sols[[s]]
    if (ncol(W_init) != K_target) next
    init_idx <- vapply(seq_len(K_target), function(k) which.min(colSums(abs(W_pool - W_init[, k]))), integer(1))
    if (length(unique(init_idx)) < K_target) {
      rem <- setdiff(seq_len(P_cands), init_idx); dups <- which(duplicated(init_idx))
      for (d in seq_along(dups)) if (length(rem) >= d) init_idx[dups[d]] <- rem[d]
    }
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

  H <- fit_fcls_h(P_obs, W)
  P_hat <- W %*% H
  res_mat <- (W - mu_vec) - U_basis %*% crossprod(U_basis, W - mu_vec)
  W_best <- W; H_best <- H
  ll_best <- sum(dbinom(X, 100, pmax(eps, pmin(1 - eps, P_hat)), log = TRUE))
  mse_best <- mean((P_obs - P_hat)^2)

  pairs <- if (K > 1) combn(K, 2, simplify = FALSE) else list()
  for (epoch in seq_len(max_epochs)) {
    flips <- 0L
    mean_h <- rowMeans(H)
    for (m in seq_len(M)) {
      if (m %in% inv_0 || m %in% inv_1) next
      x_m <- X[m, ]; row_p <- P_hat[m, ]
      ll_base <- sum(dbinom(x_m, 100, pmax(eps, pmin(1 - eps, row_p)), log = TRUE))
      best_d1 <- 0; best_k1 <- 0
      for (k in seq_len(K)) {
        if (mean_h[k] < 0.01) next
        d_bit <- 1L - 2L * W[m, k]
        p_tr <- row_p + d_bit * H[k, ]
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

# ------------------------------------------------------------------------------
# 3. PIPELINE B FUNCTIONS (UNIFIED_15PLUS STEPS 2, 3, 4)
# ------------------------------------------------------------------------------

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
      for (i in seq_len(K)) {
        curr <- selected[i]; best_c <- curr
        for (cand in setdiff(all_cands_idx, selected)) {
          t_idx <- selected; t_idx[i] <- cand
          test_loss <- eval_loss(t_idx)
          if (test_loss < (curr_loss - 1e-7)) { curr_loss <- test_loss; best_c <- cand }
        }
        selected[i] <- best_c
      }
    }
    if (curr_loss < best_overall_loss) { best_overall_loss <- curr_loss; best_overall_set <- selected }
  }
  W_pool[, best_overall_set, drop = FALSE]
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
        d_bit <- 1L - 2L * W[loc, k]
        p_trial <- row_p + d_bit * H[k, ]
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

# ------------------------------------------------------------------------------

# 1. Adaptive-Damped Slawski Extractor (damp = 0 on clean flat, 0.02 on noisy)
extract_slawski_adaptive <- function(U, mu, K, top_n = 25, is_clean = FALSE) {
  m <- nrow(U); d <- K - 1L
  pivots <- qr(t(U), tol = 1e-7)$pivot[seq_len(d)]
  U_R <- U[pivots, , drop = FALSE]; U_S <- U[setdiff(seq_len(m), pivots), , drop = FALSE]
  damp <- if (is_clean) 0.0 else 0.02 * norm(U_R, "F")^2
  inv_R <- tryCatch(solve(crossprod(U_R) + damp * diag(d)) %*% t(U_R),
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

# 2. Complete Polish v7 (Pass A 1-bit + Pass B 2-bit pairwise swaps/joint moves)
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
    # Pass A: 1-bit flips
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
    # Pass B: 2-bit pairwise swaps & joint moves
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

# 3. Load Manifest Row & Task Data
manifest_p <- if (file.exists("master_manifest.csv")) "master_manifest.csv" else "../master_manifest.csv"
man <- read.csv(manifest_p, stringsAsFactors = FALSE)
row <- man[man$run_id == task_id, ]
od  <- row$out_dir
if (!dir.exists(od) && dir.exists(file.path("..", od))) od <- file.path("..", od)

K <- as.integer(row$r_true); n <- as.integer(row$samples); m <- 365L; d <- K - 1L
W_true <- as.matrix(read.csv(file.path(od, "W_true.csv")))

set.seed(row$seed)
if (file.exists(file.path(od, "X.csv"))) {
  X <- as.matrix(read.csv(file.path(od, "X.csv")))
} else {
  alpha_v <- as.numeric(unlist(strsplit(row$alpha_str, ",")))
  xg <- matrix(rgamma(n * K, alpha_v, 1), ncol = K, byrow = TRUE)
  H_true <- t(xg / rowSums(xg))
  X <- matrix(rbinom(m * n, 100, as.vector(W_true %*% H_true)), nrow = m)
}
P_obs <- X / 100.0

# 4. Run Upgraded Unified_15plus_v3 (Steps 1 -> 2 -> 3A Snap -> 3B Exchange -> 4 SNR-Gated Polish v7)
set.seed(row$seed)
mu_b <- rowMeans(P_obs)
sv_b <- svd(sweep(P_obs, 1, mu_b, "-"), nu = K, nv = 0)
U_b  <- sv_b$u[, seq_len(d), drop = FALSE]
snr_gap <- if (length(sv_b$d) >= K && sv_b$d[K] > 1e-9) sv_b$d[d] / sv_b$d[K] else 999.0

W_pool_raw <- harvest_dual_20svd(P_obs, K)
# Step 3A: Basis-Proximity Snapping (1-bit & 2-bit descent onto U_b)
W_pool_snap <- snap_pool_to_basis(W_pool_raw, U_b, mu_b)

# Step 3B: 3-Start Coordinate Exchange on Snapped Pool (unbiased lambda_pars = 0 vs 5e-5)
W_v3_s3_unb  <- coord_exchange_robust_parsimony(P_obs, W_pool_snap, K, lambda_pars = 0.0)
W_v3_s3_pars <- coord_exchange_robust_parsimony(P_obs, W_pool_snap, K, lambda_pars = 5e-5)

# Step 4: SNR-Gated Polish v7 (Pass A + Pass B; skips if snr_gap <= 1.12 i.e. sub-BBP flat noise)
W_v3_s4_unb  <- if (snr_gap > 1.12) polish_bits_v7_full(X, P_obs, W_v3_s3_unb,  U_b, mu_b) else W_v3_s3_unb
W_v3_s4_pars <- if (snr_gap > 1.12) polish_bits_v7_full(X, P_obs, W_v3_s3_pars, U_b, mu_b) else W_v3_s3_pars

eval_mat <- function(Wp) {
  if (is.null(Wp) || ncol(Wp) != K) return(c(err = NA_real_, rec = NA_real_))
  dm <- matrix(0, K, K); for (i in 1:K) dm[i, ] <- colSums(Wp[, i] != W_true)
  hd <- hungarian_dists(dm)
  c(err = mean(hd), rec = mean(hd == 0) * 100.0)
}

eval_file <- function(fp) if (file.exists(fp)) eval_mat(as.matrix(read.csv(fp))) else c(err = NA_real_, rec = NA_real_)

pa_s3 <- eval_file(file.path(od, sprintf("W_alt_sol4_snap_subset_1_r%d.csv", K)))
pa_s4 <- eval_file(file.path(od, sprintf("W_alt_sol4_snap_polish_v7_1_r%d.csv", K)))
pb_v2 <- eval_file(file.path(od, "W_pred_unified_15plus.csv"))

v3_s3_u <- eval_mat(W_v3_s3_unb);  v3_s4_u <- eval_mat(W_v3_s4_unb)
v3_s3_p <- eval_mat(W_v3_s3_pars); v3_s4_p <- eval_mat(W_v3_s4_pars)

out_dir <- "../results/unified_v3_benchmark"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(data.frame(
  task_id = task_id, K = K, regime = row$alpha_str, n = n, snr_gap = snr_gap,
  pa_s3_err = pa_s3["err"], pa_s4_err = pa_s4["err"], pa_s4_rec = pa_s4["rec"],
  pb_v2_err = pb_v2["err"], pb_v2_rec = pb_v2["rec"],
  v3_s3_unb_err = v3_s3_u["err"], v3_s4_unb_err = v3_s4_u["err"], v3_s4_unb_rec = v3_s4_u["rec"],
  v3_s3_prs_err = v3_s3_p["err"], v3_s4_prs_err = v3_s4_p["err"], v3_s4_prs_rec = v3_s4_p["rec"],
  stringsAsFactors = FALSE
), file.path(out_dir, sprintf("v3_task_%03d.rds", task_id)))

cat(sprintf("[DONE] Task %3d (K=%2d, n=%2d, SNR_Gap=%.2f) | PipeA: %.2f | PipeB_v2: %.2f | v3_Unb: %.2f | v3_Pars: %.2f\n",
            task_id, K, n, snr_gap, pa_s4["err"], pb_v2["err"], v3_s4_u["err"], v3_s4_p["err"]))
