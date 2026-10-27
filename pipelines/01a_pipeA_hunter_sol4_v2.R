for (p in c("phase7/R/helper-functions.R", "../phase7/R/helper-functions.R", "~/phase7/R/helper-functions.R")) if (file.exists(p)) { source(p); break }
for (p in c("phase7/R/nmf-functions_v2.R", "../phase7/R/nmf-functions_v2.R", "~/phase7/R/nmf-functions_v2.R")) if (file.exists(p)) { source(p); break }
for (p in c("phase7/R/simulation-functions.R", "../phase7/R/simulation-functions.R", "~/phase7/R/simulation-functions.R")) if (file.exists(p)) { source(p); break }

library(utils)
library(quadprog)
library(MASS)

fit_fcls_quadprog <- function(P_obs, W, prune_tol = 1e-4) {
  m <- nrow(W)
  K <- ncol(W)
  n <- ncol(P_obs)

  H <- matrix(0, nrow = K, ncol = n)
  Dmat <- crossprod(W) + diag(1e-8, K)
  Amat <- cbind(rep(1, K), diag(K))
  bvec <- c(1, rep(0, K))
  meq <- 1

  for (j in seq_len(n)) {
    dvec <- as.vector(crossprod(W, P_obs[, j]))
    qp_fit <- tryCatch({
      solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = meq)
    }, error = function(e) list(solution = rep(1 / K, K)))
    h_sol <- pmax(0, qp_fit$solution)
    H[, j] <- h_sol / sum(h_sol)
  }

  active_idx <- which(rowSums(H) > prune_tol)
  if (length(active_idx) == 0) active_idx <- seq_len(K)

  W_pruned <- W[, active_idx, drop = FALSE]
  H_pruned <- H[active_idx, , drop = FALSE]
  H_pruned <- sweep(H_pruned, 2, colSums(H_pruned), "/")

  residuals <- P_obs - (W_pruned %*% H_pruned)
  return(list(W = W_pruned, H = H_pruned, mse = mean(residuals^2)))
}

# Weighted affine subspace decomposition
fit_weighted_affine_subspace <- function(X, N_vec, K) {
  m <- nrow(X)
  n <- ncol(X)
  Nmat <- matrix(rep(N_vec, each = m), nrow = m, ncol = n)
  P <- X / Nmat

  mu <- rowMeans(P)
  P_centered <- sweep(P, 1, mu, "-")

  # Exponential leverage weights
  w <- rexp(m, rate = 1.0)
  w <- w / mean(w)
  sqrt_w <- sqrt(w)

  P_weighted <- sweep(P_centered, 1, sqrt_w, "*")
  sv <- svd(P_weighted, nu = K - 1, nv = 0)

  # Unweight basis coordinates back to Euclidean metric
  U_unweighted <- sweep(sv$u, 1, sqrt_w, "/")
  qU <- qr(U_unweighted)
  U <- qr.Q(qU)[, seq_len(K - 1), drop = FALSE]

  return(list(U = U, mu = mu))
}

# Condition-Gated MaxVol Pivot Selection (Prevents condition number explosion)
select_pivots_maxvol <- function(U, K, max_iter = 5) {
  d <- max(1L, K - 1L)
  M <- nrow(U)

  qr_res <- qr(t(U), tol = 1e-7)
  pivots <- qr_res$pivot[seq_len(d)]
  U_R <- U[pivots, , drop = FALSE]
  curr_rcond <- rcond(U_R)
  if (is.na(curr_rcond)) curr_rcond <- 0

  leverage_order <- order(rowSums(U^2), decreasing = TRUE)
  cand_pool <- leverage_order[seq_len(min(M, 40))]

  for (iter in seq_len(max_iter)) {
    improved <- FALSE
    for (i in seq_len(d)) {
      for (cand in cand_pool) {
        if (cand %in% pivots) next
        trial_pivots <- pivots
        trial_pivots[i] <- cand
        trial_UR <- U[trial_pivots, , drop = FALSE]
        trial_rcond <- rcond(trial_UR)

        if (!is.na(trial_rcond) && trial_rcond > curr_rcond * 1.20) {
          pivots <- trial_pivots
          curr_rcond <- trial_rcond
          improved <- TRUE
          break
        }
      }
    }
    if (!improved) break
  }

  return(pivots)
}

# Regularized Slawski Vertex Generation Engine
extract_slawski_vertex_set <- function(U, mu, K) {
  m <- length(mu)
  d <- max(1L, K - 1L)

  # 1. Condition-gated MaxVol pivot selection
  R_idx <- select_pivots_maxvol(U, K)
  S_idx <- setdiff(seq_len(m), R_idx)

  U_R <- U[R_idx, , drop = FALSE]
  U_S <- U[S_idx, , drop = FALSE]

  # 2. Adaptive Tikhonov pseudo-inverse (M1 - suppresses projection overshoot)
  damp_val <- 0.02 * norm(U_R, "F")^2
  inv_damped <- tryCatch({
    solve(crossprod(U_R) + damp_val * diag(d)) %*% t(U_R)
  }, error = function(e) {
    MASS::ginv(U_R)
  })

  G_hat <- U_S %*% inv_damped
  h_hat <- mu[S_idx] - G_hat %*% mu[R_idx]

  # 3. Boolean hypercube enumeration over well-conditioned pivots
  b_grid <- as.matrix(expand.grid(rep(list(c(0L, 1L)), d)))
  b_mat <- t(b_grid)
  z_S_hat <- sweep(G_hat %*% b_mat, 1, as.vector(h_hat), "+")

  # 4. Discrete hypercube reconstruction
  W_reconstructed <- matrix(0L, nrow = m, ncol = ncol(b_mat))
  W_reconstructed[R_idx, ] <- b_mat
  W_reconstructed[S_idx, ] <- (z_S_hat >= 0.5) * 1L

  unique_cols <- !duplicated(t(W_reconstructed))
  W_unique <- W_reconstructed[, unique_cols, drop = FALSE]

  # 5. Order candidates by proximity to estimated affine subspace
  W_centered <- sweep(W_unique, 1, mu, "-")
  W_proj <- U %*% (crossprod(U, W_centered))
  dist_sq <- colSums((W_centered - W_proj)^2)

  best_idx <- order(dist_sq)[seq_len(min(K, ncol(W_unique)))]
  return(W_unique[, best_idx, drop = FALSE])
}

evaluate_set_difference <- function(W_new, W_old) {
  K_new <- ncol(W_new)
  K_old <- ncol(W_old)
  dist_matrix <- matrix(0, nrow = K_new, ncol = K_old)
  for (i in seq_len(K_new)) {
    dist_matrix[i, ] <- colSums(W_new[, i] != W_old)
  }
  min_new_to_old <- apply(dist_matrix, 1, min)
  min_old_to_new <- apply(dist_matrix, 2, min)
  all_distances <- c(min_new_to_old, min_old_to_new)
  return(list(
    overlap_pct = sum(min_new_to_old <= 1L) / K_new,
    avg_distance = mean(all_distances),
    max_distance = max(all_distances)
  ))
}

run_hunter_sol4 <- function(X, N_vec, K, N_consecutive_stop = 100, max_eval = 1000,
                            max_overlap_tol = 0.60, mse_margin = 1.25) {
  m <- nrow(X)
  n <- ncol(X)
  Nmat <- matrix(rep(N_vec, each = m), nrow = m, ncol = n)
  P_obs <- X / Nmat
  Phat <- pmin(pmax(P_obs, 1e-12), 1 - 1e-12)
  min_avg_dist_tol <- max(15.0, 2.5 * K)

  fit_base <- fit_affine_subspace(X, N_vec, K, transform = "none", eps_var = 1e-6)
  W_base <- extract_slawski_vertex_set(fit_base$U, fit_base$mu, K)
  fit_score_base <- fit_fcls_quadprog(P_obs, W_base)

  best_mse <- fit_score_base$mse
  saved_solutions <- list(list(
    id = 1, W = fit_score_base$W, H = fit_score_base$H, mse = best_mse,
    overlap_pct = 0.0, avg_distance = 0.0, max_distance = 0.0, source = "Original_Baseline"
  ))

  consecutive_no_sol <- 0
  total_bootstraps <- 0

  while (consecutive_no_sol < N_consecutive_stop && total_bootstraps < max_eval) {
    total_bootstraps <- total_bootstraps + 1
    X_star <- matrix(rbinom(m * n, size = as.vector(Nmat), prob = as.vector(Phat)), nrow = m, ncol = n)

    # Generate weighted affine plane
    fit_star <- fit_weighted_affine_subspace(X_star, N_vec, K)
    if (is.null(fit_star$U) || ncol(fit_star$U) < (K - 1)) next

    W_star <- extract_slawski_vertex_set(fit_star$U, fit_star$mu, K)
    is_different <- TRUE
    last_diff <- NULL

    for (s in seq_along(saved_solutions)) {
      diff <- evaluate_set_difference(W_star, saved_solutions[[s]]$W)
      last_diff <- diff
      if (diff$overlap_pct > max_overlap_tol || diff$avg_distance < min_avg_dist_tol) {
        is_different <- FALSE
        break
      }
    }

    if (!is_different) {
      consecutive_no_sol <- consecutive_no_sol + 1
      next
    }

    fit_score_star <- fit_fcls_quadprog(P_obs, W_star)
    if (fit_score_star$mse <= (best_mse * mse_margin)) {
      if (fit_score_star$mse < best_mse) best_mse <- fit_score_star$mse
      sol_idx <- length(saved_solutions) + 1
      saved_solutions[[sol_idx]] <- list(
        id = sol_idx, W = fit_score_star$W, H = fit_score_star$H, mse = fit_score_star$mse,
        overlap_pct = last_diff$overlap_pct, avg_distance = last_diff$avg_distance,
        max_distance = last_diff$max_distance, source = paste0("Weighted_Boot_", total_bootstraps)
      )
      consecutive_no_sol <- 0
    } else {
      consecutive_no_sol <- consecutive_no_sol + 1
    }
  }
  return(saved_solutions)
}

# Main Execution
args <- commandArgs(trailingOnly = TRUE)
task_id <- as.integer(args[1])
manifest <- read.csv("master_manifest.csv", stringsAsFactors = FALSE)
row <- manifest[task_id, ]
set.seed(row$seed)

draw_dirichlet <- function(n, alpha) {
  k <- length(alpha)
  x <- matrix(rgamma(n * k, alpha, 1), ncol = k, byrow = TRUE)
  return(x / as.vector(x %*% rep(1, k)))
}

alpha_vec <- as.numeric(unlist(strsplit(row$alpha_str, ",")))
W_true <- make_endmembers(365, row$r_true)
H_true <- t(draw_dirichlet(row$samples, alpha_vec))
Nmat <- matrix(100, nrow = 365, ncol = row$samples)
X <- matrix(rbinom(365 * row$samples, size = as.vector(Nmat), prob = as.vector(W_true %*% H_true)), nrow = 365)

solutions <- run_hunter_sol4(X, rep(100, row$samples), row$r_true)

summary_records <- data.frame()
for (s in seq_along(solutions)) {
  sol <- solutions[[s]]
  write.csv(sol$W, file.path(row$out_dir, sprintf("W_alt_sol4_%d_r%d.csv", sol$id, row$r_true)), row.names = FALSE)
  summary_records <- rbind(summary_records, data.frame(
    solution_id = sol$id, source = sol$source, mse_score = sol$mse,
    avg_distance = sol$avg_distance, max_distance = sol$max_distance, num_strains = ncol(sol$W)
  ))
}
write.csv(summary_records, file.path(row$out_dir, "summary_hunter_sol4.csv"), row.names = FALSE)
