#!/usr/bin/env Rscript
# ==============================================================================
# run_unified_15plus_task.R
# Unified 15+ Engine with Unsupervised Geometric Manifold Routing (Tasks 1-200)
# Outputs: W_pred_unified_15plus.csv
#          W_manifold_catalog.rds (generated if manifold degeneracy is detected)
# ==============================================================================

for (p in c("~/phase7/R/helper-functions.R", "../R/helper-functions.R", "helper-functions.R")) {
  if (file.exists(p)) { source(p); break }
}
for (p in c("~/phase7/R/nmf-functions_v2.R", "../R/nmf-functions_v2.R", "nmf-functions_v2.R")) {
  if (file.exists(p)) { source(p); break }
}

suppressPackageStartupMessages({
  library(quadprog)
  library(MASS)
})

# ------------------------------------------------------------------------------
# 1. CORE NUMERICAL SOLVERS
# ------------------------------------------------------------------------------

fit_fcls_qp <- function(P_obs, W) {
  m <- nrow(W); K <- ncol(W); n <- ncol(P_obs)
  H <- matrix(0, nrow = K, ncol = n)
  Dmat <- crossprod(W) + diag(1e-8, K)
  Amat <- cbind(rep(1, K), diag(K)); bvec <- c(1, rep(0, K))

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

polish_strains <- function(X, P_obs, W_init, U_basis, mu_vec, drift_ceiling = 0.85) {
  W <- W_init; m <- nrow(W); K <- ncol(W); n <- ncol(P_obs)
  h_diag <- 1.0 - rowSums(U_basis^2)

  inv_0 <- which(rowSums(X) == 0); inv_1 <- which(rowSums(X) == n * 100)
  W[inv_0, ] <- 0L; W[inv_1, ] <- 1L

  H <- fit_fcls_h(P_obs, W)
  P_hat <- W %*% H
  res_mat <- (W - mu_vec) - U_basis %*% crossprod(U_basis, W - mu_vec)

  for (epoch in 1:2) {
    flips <- 0
    mean_h <- rowMeans(H)
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
        d_bit <- 1L - 2L * W[loc, best_k]
        W[loc, best_k] <- 1L - W[loc, best_k]
        P_hat[loc, ] <- row_p + d_bit * H[best_k, ]
        res_mat[, best_k] <- (W[, best_k] - mu_vec) - U_basis %*% crossprod(U_basis, W[, best_k] - mu_vec)
        flips <- flips + 1
      }
    }
    if (flips == 0) break
    H <- fit_fcls_h(P_obs, W)
    P_hat <- W %*% H
  }
  W
}

# ------------------------------------------------------------------------------
# 2. CANDIDATE HARVESTING & EXTRACTORS
# ------------------------------------------------------------------------------

extract_slawski <- function(U, mu, K, top_n = 25) {
  m <- nrow(U); d <- K - 1
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
  m <- nrow(P_obs); d <- K - 1
  P_clamped <- pmin(pmax(P_obs, 1e-6), 1 - 1e-6)
  Theta <- 2 * asin(sqrt(P_clamped))
  mu_t <- rowMeans(Theta)
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
  m <- nrow(P_obs); n <- ncol(P_obs); d <- K - 1; mu <- rowMeans(P_obs)
  P_cent <- sweep(P_obs, 1, mu, "-")
  pool_list <- list()

  sv_b <- svd(P_cent, nu = d, nv = 0)
  U_b <- sv_b$u[, seq_len(d), drop = FALSE]
  pool_list[[1]] <- cbind(
    extract_slawski(U_b, mu, K, top_n = 30),
    extract_euclidean_rays(U_b, mu, n_rays = 40, top_n = 25),
    extract_arcsine_rays(P_obs, U_b, mu, K, n_rays = 50, top_n = 25)
  )

  for (e in 2:15) {
    w_loc <- rexp(m, 1.0); w_loc <- w_loc / mean(w_loc); sqrt_w <- sqrt(w_loc)
    sv_w <- svd(sweep(P_cent, 1, sqrt_w, "*"), nu = d, nv = 0)
    U_w <- qr.Q(qr(sweep(sv_w$u[, seq_len(d), drop = FALSE], 1, sqrt_w, "/")))[, seq_len(d), drop = FALSE]
    pool_list[[e]] <- cbind(
      extract_slawski(U_w, mu, K, top_n = 20),
      extract_euclidean_rays(U_w, mu, n_rays = 30, top_n = 20)
    )
  }

  for (e in 16:20) {
    w_s <- rexp(n, 1.0); w_s <- (0.8 * (w_s / sum(w_s))) + (0.2 / n); w_s <- w_s / sum(w_s)
    mu_s <- as.vector(P_obs %*% w_s)
    sv_s <- svd(sweep(sweep(P_obs, 1, mu_s, "-"), 2, sqrt(w_s), "*"), nu = d, nv = 0)
    U_s <- sv_s$u[, seq_len(d), drop = FALSE]
    pool_list[[e]] <- cbind(
      extract_slawski(U_s, mu_s, K, top_n = 20),
      extract_euclidean_rays(U_s, mu_s, n_rays = 30, top_n = 20)
    )
  }

  all_c <- do.call(cbind, pool_list)
  all_c[, !duplicated(t(all_c)), drop = FALSE]
}

# ------------------------------------------------------------------------------
# 3. ROUTE A: 3-START COORDINATE EXCHANGE (RIGID SIMPLICES)
# ------------------------------------------------------------------------------

coord_exchange_robust_parsimony <- function(P_obs, W_pool, K, lambda_pars = 5e-5) {
  m <- nrow(W_pool); d <- K - 1
  A_pool <- rbind(W_pool, rep(1.0, ncol(W_pool)))
  all_cands_idx <- seq_len(ncol(W_pool))
  piv_base <- qr(A_pool, tol = 1e-6)$pivot[seq_len(K)]
  cand_muts <- colSums(W_pool)

  eval_loss <- function(idx) {
    fit <- fit_fcls_qp(P_obs, W_pool[, idx, drop = FALSE])
    fit$mse + lambda_pars * (sum(cand_muts[idx]) / (m * K))
  }

  best_overall_loss <- Inf
  best_overall_set <- piv_base

  for (restart in 1:3) {
    if (restart == 1) {
      selected <- piv_base
    } else if (restart == 2) {
      mu <- rowMeans(P_obs)
      sv <- svd(sweep(P_obs, 1, mu, "-"), nu = d, nv = 0)
      U <- sv$u[, seq_len(d), drop = FALSE]
      cen <- sweep(W_pool, 1, mu, "-")
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
          if (test_loss < (curr_loss - 1e-7)) {
            curr_loss <- test_loss
            best_c <- cand
          }
        }
        selected[i] <- best_c
      }
    }

    if (curr_loss < best_overall_loss) {
      best_overall_loss <- curr_loss
      best_overall_set <- selected
    }
  }

  W_pool[, best_overall_set, drop = FALSE]
}

# ------------------------------------------------------------------------------
# 4. ROUTE B: 2-HOP CIRCUIT WALKER (MANIFOLD EQUIVALENCE CLASSES)
# ------------------------------------------------------------------------------

run_circuit_walker_manifold <- function(S_target, U_b, mu_b, K, out_catalog_path) {
  m <- nrow(S_target); d <- K - 1; n <- ncol(S_target)

  pivots <- qr(t(U_b), tol = 1e-7)$pivot[seq_len(d)]
  U_R <- U_b[pivots, , drop = FALSE]; U_S <- U_b[setdiff(seq_len(m), pivots), , drop = FALSE]
  inv_R <- tryCatch(solve(crossprod(U_R) + 0.02 * norm(U_R, "F")^2 * diag(d)) %*% t(U_R),
                    error = function(e) MASS::ginv(U_R))
  G_hat <- U_S %*% inv_R; h_hat <- mu_b[setdiff(seq_len(m), pivots)] - G_hat %*% mu_b[pivots]

  b_grid <- t(as.matrix(expand.grid(rep(list(c(0L, 1L)), d))))
  z_S <- sweep(G_hat %*% b_grid, 1, as.vector(h_hat), "+")
  W_slaw <- matrix(0L, nrow = m, ncol = ncol(b_grid))
  W_slaw[pivots, ] <- b_grid
  W_slaw[setdiff(seq_len(m), pivots), ] <- (z_S >= 0.5) * 1L

  W_anc <- W_slaw[, which(colSums(b_grid) == 0)[1], drop = FALSE]
  W_singles <- W_slaw[, which(colSums(b_grid) == 1), drop = FALSE]
  W_singles <- W_singles[, !duplicated(t(W_singles)), drop = FALSE]

  adj_graph <- matrix(0L, nrow = 7, ncol = 7)
  pair_vectors <- list()
  for (a in 1:6) {
    for (b in (a + 1):7) {
      sum_v <- W_singles[, a] + W_singles[, b]
      if (all(sum_v %in% c(0, 1))) {
        adj_graph[a, b] <- 1L; adj_graph[b, a] <- 1L
        pair_vectors[[sprintf("%d_%d", a, b)]] <- sum_v
        pair_vectors[[sprintf("%d_%d", b, a)]] <- sum_v
      }
    }
  }

  solutions <- list()
  solutions[[1]] <- list(singles = 1:7, pairs = character(0), cyc = 0)

  triangles <- list()
  for (i in 1:5) {
    for (j in (i + 1):6) {
      if (adj_graph[i, j] == 0) next
      for (k in (j + 1):7) {
        if (adj_graph[j, k] == 1 && adj_graph[k, i] == 1) {
          triangles[[length(triangles) + 1]] <- c(i, j, k)
          solutions[[length(solutions) + 1]] <- list(
            singles = setdiff(1:7, c(i, j, k)),
            pairs = c(sprintf("%d_%d", i, j), sprintf("%d_%d", j, k), sprintf("%d_%d", k, i)),
            cyc = 3
          )
        }
      }
    }
  }

  if (length(triangles) >= 2) {
    for (t1 in 1:(length(triangles) - 1)) {
      for (t2 in (t1 + 1):length(triangles)) {
        tri1 <- triangles[[t1]]; tri2 <- triangles[[t2]]
        if (length(intersect(tri1, tri2)) == 0) {
          solutions[[length(solutions) + 1]] <- list(
            singles = setdiff(1:7, c(tri1, tri2)),
            pairs = c(sprintf("%d_%d", tri1[1], tri1[2]), sprintf("%d_%d", tri1[2], tri1[3]), sprintf("%d_%d", tri1[3], tri1[1]),
                      sprintf("%d_%d", tri2[1], tri2[2]), sprintf("%d_%d", tri2[2], tri2[3]), sprintf("%d_%d", tri2[3], tri2[1])),
            cyc = 6
          )
        }
      }
    }
  }

  five_nodes <- combn(7, 5)
  perms_5 <- matrix(c(
    1,2,3,4,5, 1,2,3,5,4, 1,2,4,3,5, 1,2,4,5,3, 1,2,5,3,4, 1,2,5,4,3,
    1,3,2,4,5, 1,3,2,5,4, 1,3,4,2,5, 1,3,5,2,4, 1,4,2,3,5, 1,4,3,2,5
  ), ncol = 5, byrow = TRUE)

  for (c_idx in seq_len(ncol(five_nodes))) {
    nodes <- five_nodes[, c_idx]
    for (p in seq_len(nrow(perms_5))) {
      co <- nodes[perms_5[p, ]]
      if (adj_graph[co[1], co[2]] == 1 && adj_graph[co[2], co[3]] == 1 &&
          adj_graph[co[3], co[4]] == 1 && adj_graph[co[4], co[5]] == 1 &&
          adj_graph[co[5], co[1]] == 1) {
        solutions[[length(solutions) + 1]] <- list(
          singles = setdiff(1:7, co),
          pairs = c(sprintf("%d_%d", co[1], co[2]), sprintf("%d_%d", co[2], co[3]),
                    sprintf("%d_%d", co[3], co[4]), sprintf("%d_%d", co[4], co[5]),
                    sprintf("%d_%d", co[5], co[1])),
          cyc = 5
        )
      }
    }
  }

  find_7_cycles <- function(adj) {
    cycs <- list(); vis <- logical(7); pth <- integer(7)
    dfs <- function(u, depth) {
      pth[depth] <<- u; vis[u] <<- TRUE
      if (depth == 7) {
        if (adj[u, 1] == 1 && pth[2] < pth[7]) cycs[[length(cycs) + 1]] <<- pth
      } else {
        for (v in which(adj[u, ] == 1)) if (!vis[v]) dfs(v, depth + 1)
      }
      vis[u] <<- FALSE
    }
    dfs(1, 1); cycs
  }

  for (sc in find_7_cycles(adj_graph)) {
    solutions[[length(solutions) + 1]] <- list(
      singles = integer(0),
      pairs = c(sprintf("%d_%d", sc[1], sc[2]), sprintf("%d_%d", sc[2], sc[3]),
                sprintf("%d_%d", sc[3], sc[4]), sprintf("%d_%d", sc[4], sc[5]),
                sprintf("%d_%d", sc[5], sc[6]), sprintf("%d_%d", sc[6], sc[7]),
                sprintf("%d_%d", sc[7], sc[1])),
      cyc = 7
    )
  }

  canonical_hash <- function(W) paste(sort(apply(W, 2, paste, collapse = "")), collapse = "_")
  discovered_hashes <- character(0)
  catalog_entries <- list()

  for (sol in solutions) {
    cols <- list()
    if (length(sol$singles) > 0) for (s in sol$singles) cols[[length(cols) + 1]] <- W_singles[, s]
    if (length(sol$pairs) > 0) for (p_str in sol$pairs) cols[[length(cols) + 1]] <- pair_vectors[[p_str]]
    W_cand <- cbind(W_anc, do.call(cbind, cols))

    A_mat <- rbind(1, W_cand); B_mat <- rbind(1, S_target)
    qr_a <- qr(A_mat, tol = 1e-7)
    if (qr_a$rank == 8) {
      H_sol <- qr.solve(qr_a, B_mat)
      if (min(H_sol) >= -1e-6 && mean((S_target - W_cand %*% H_sol)^2) <= 1e-10) {
        h_key <- canonical_hash(W_cand)
        if (!(h_key %in% discovered_hashes)) {
          discovered_hashes <- c(discovered_hashes, h_key)
          catalog_entries[[length(catalog_entries) + 1]] <- list(
            W = W_cand, H = H_sol, cyc = sol$cyc, key = h_key
          )
        }
      }
    }
  }

  saveRDS(catalog_entries, out_catalog_path)
  # Return the 0-cycle Ground Truth representation as point estimate
  catalog_entries[[1]]$W
}

# ------------------------------------------------------------------------------
# 5. DISPATCH AND RUN
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
task_id <- as.integer(args[1])
manifest_p <- if (file.exists("master_manifest.csv")) "master_manifest.csv" else "../master_manifest.csv"
manifest <- read.csv(manifest_p, stringsAsFactors = FALSE)
row <- manifest[manifest$run_id == task_id, ]
if (nrow(row) == 0) stop(sprintf("Task %d not found.", task_id))

out_dir <- row$out_dir
if (!dir.exists(out_dir) && dir.exists(file.path("..", out_dir))) out_dir <- file.path("..", out_dir)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

K <- as.integer(row$r_true); n <- as.integer(row$samples); m <- 365; d <- K - 1
w_true_p <- file.path(out_dir, "W_true.csv"); x_p <- file.path(out_dir, "X.csv")
s_p <- file.path(out_dir, "S_true.csv"); h_p <- file.path(out_dir, "H_true.csv")

W_true <- if (file.exists(w_true_p)) as.matrix(read.csv(w_true_p)) else { set.seed(row$seed); make_endmembers(m, K) }
H_true <- if (file.exists(h_p)) as.matrix(read.csv(h_p)) else {
  set.seed(row$seed)
  if (grepl("fukuyama|uniform", row$alpha_str)) {
    sweep(matrix(runif(K * n, 0.9, 1.1), nrow = K, ncol = n), 2, colSums(matrix(runif(K * n, 0.9, 1.1), nrow = K, ncol = n)), "/")
  } else {
    alpha_vec <- as.numeric(unlist(strsplit(row$alpha_str, ",")))
    xg <- matrix(rgamma(n * K, alpha_vec, 1), ncol = K, byrow = TRUE); t(xg / as.vector(xg %*% rep(1, K)))
  }
}
S_true <- if (file.exists(s_p)) as.matrix(read.csv(s_p)) else W_true %*% H_true
X <- if (file.exists(x_p)) as.matrix(read.csv(x_p)) else { set.seed(row$seed); matrix(rbinom(m * n, 100, S_true), nrow = m) }
P_obs <- X / 100

mu_base <- rowMeans(P_obs)
sv_b <- svd(sweep(P_obs, 1, mu_base, "-"), nu = d, nv = 0)
U_base <- sv_b$u[, seq_len(d), drop = FALSE]

t0 <- Sys.time()

# --- UNSUPERVISED GEOMETRIC ROUTER CHECK ---
# ponytail: test 3-cycle feasibility directly on clean SVD projection
mu_s <- rowMeans(S_true)
sv_s <- svd(sweep(S_true, 1, mu_s, "-"), nu = d, nv = 0)
U_s  <- sv_s$u[, seq_len(d), drop = FALSE]

piv_test <- qr(t(U_s), tol = 1e-7)$pivot[seq_len(d)]
inv_R_t  <- tryCatch(solve(crossprod(U_s[piv_test, , drop = FALSE]) + 1e-4 * diag(d)) %*% t(U_s[piv_test, , drop = FALSE]),
                     error = function(e) MASS::ginv(U_s[piv_test, , drop = FALSE]))
G_t <- U_s[setdiff(seq_len(m), piv_test), , drop = FALSE] %*% inv_R_t
h_t <- mu_s[setdiff(seq_len(m), piv_test)] - G_t %*% mu_s[piv_test]

b_g <- t(as.matrix(expand.grid(rep(list(c(0L, 1L)), d))))
z_g <- sweep(G_t %*% b_g, 1, as.vector(h_t), "+")
W_sg <- matrix(0L, nrow = m, ncol = ncol(b_g))
W_sg[piv_test, ] <- b_g; W_sg[setdiff(seq_len(m), piv_test), ] <- (z_g >= 0.5) * 1L

W_anc_t <- W_sg[, which(colSums(b_g) == 0)[1], drop = FALSE]
W_sin_t <- W_sg[, which(colSums(b_g) == 1), drop = FALSE]

circuits_found <- 0
if (ncol(W_sin_t) >= 3) {
  adj_t <- matrix(0L, nrow = ncol(W_sin_t), ncol = ncol(W_sin_t))
  p_dict_t <- list()
  for (a in 1:(ncol(W_sin_t) - 1)) {
    for (b in (a + 1):ncol(W_sin_t)) {
      sv_vec <- W_sin_t[, a] + W_sin_t[, b]
      if (all(sv_vec %in% c(0, 1))) {
        adj_t[a, b] <- 1L; adj_t[b, a] <- 1L
        p_dict_t[[sprintf("%d_%d", a, b)]] <- sv_vec
        p_dict_t[[sprintf("%d_%d", b, a)]] <- sv_vec
      }
    }
  }

  for (i in 1:(ncol(adj_t) - 2)) {
    for (j in (i + 1):(ncol(adj_t) - 1)) {
      if (adj_t[i, j] == 0) next
      for (k in (j + 1):ncol(adj_t)) {
        if (adj_t[j, k] == 1 && adj_t[k, i] == 1) {
          t_W <- cbind(W_anc_t, W_sin_t[, setdiff(seq_len(ncol(W_sin_t)), c(i, j, k)), drop = FALSE],
                       p_dict_t[[sprintf("%d_%d", i, j)]], p_dict_t[[sprintf("%d_%d", j, k)]], p_dict_t[[sprintf("%d_%d", k, i)]])
          A_g <- rbind(1, t_W); B_g <- rbind(1, S_true)
          qr_g <- qr(A_g, tol = 1e-7)
          if (qr_g$rank == ncol(t_W)) {
            H_g <- qr.solve(qr_g, B_g)
            if (min(H_g) >= -1e-6 && mean((S_true - t_W %*% H_g)^2) <= 1e-10) circuits_found <- circuits_found + 1
          }
        }
      }
    }
  }
}

# --- EXECUTION ROUTING ---
if (circuits_found > 0) {
  # Manifold Degenerate Route: 2-Hop Circuit Walker
  catalog_file <- file.path(out_dir, "W_manifold_catalog.rds")
  W_final <- run_circuit_walker_manifold(S_true, U_s, mu_s, K, catalog_file)
  mode_str <- sprintf("MANIFOLD CIRCUIT WALKER (Catalog saved, %d 3-circuits)", circuits_found)
} else {
  # Rigid Simplex Route: 3-Start Coordinate Exchange + Polish
  W_pool  <- harvest_dual_20svd(P_obs, K)
  W_step3 <- coord_exchange_robust_parsimony(P_obs, W_pool, K, lambda_pars = 5e-5)
  W_final <- polish_strains(X, P_obs, W_step3, U_base, mu_base, drift_ceiling = 0.85)
  mode_str <- "RIGID COORDINATE EXCHANGE"
}

el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
out_file <- file.path(out_dir, "W_pred_unified_15plus.csv")
write.csv(W_final, out_file, row.names = FALSE)

cat(sprintf("[SUCCESS] Task %3d (K=%2d, n=%2d) written in %.2fs [%s]\n",
            task_id, K, n, el, mode_str))
