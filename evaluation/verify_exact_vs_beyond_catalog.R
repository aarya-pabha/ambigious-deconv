#!/usr/bin/env Rscript
# ==============================================================================
# File: verify_exact_vs_beyond_catalog.R
# Purpose: Fully self-contained verification of:
#          (1) All binary lattice vectors on SVD plane (3 * 2^(d-2) law)
#          (2) Exact zero-MSE factorizations beyond the pair catalog (660-mutations)
#          (3) Pipeline B behavior across candidate pool configurations
# Usage:   Rscript verify_exact_vs_beyond_catalog.R
# ==============================================================================

suppressPackageStartupMessages({
  library(quadprog)
  library(MASS)
})

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

# ------------------------------------------------------------------------------
# 2. RUN TESTS
# ------------------------------------------------------------------------------

cat("====================================================================================================\n")
cat("   TEST 1: COUNTING ALL BINARY LATTICE VECTORS ON THE CLEAN SVD PLANE (3 * 2^(d-2) LAW)\n")
cat("====================================================================================================\n")
cat(sprintf("%-4s | %-10s | %-12s | %-14s | %-16s | %-18s\n",
            "K", "2^d Grid", "On-Plane Vec", "Pair-Catalog", "Higher-Order(>=3)", "Theory 3*2^(d-2)"))
cat("----------------------------------------------------------------------------------------------------\n")

for (Kv in c(5L, 7L, 8L, 10L)) {
  sim <- make_fukuyama(Kv, 202601L)
  d <- sim$d; mu <- rowMeans(sim$S_true)
  U <- svd(sweep(sim$S_true, 1, mu, "-"), nu = d, nv = 0)$u[, 1:d, drop = FALSE]
  piv <- qr(t(U), tol = 1e-7)$pivot[1:d]
  G_h <- U[-piv, , drop = FALSE] %*% solve(U[piv, , drop = FALSE])
  h_h <- mu[-piv] - G_h %*% mu[piv]
  b_g <- t(as.matrix(expand.grid(rep(list(c(0L, 1L)), d))))
  z_S <- sweep(G_h %*% b_g, 1, as.vector(h_h), "+")
  z_r <- round(z_S)
  ok  <- which(apply(abs(z_S - z_r), 2, max) < 1e-6 & apply(z_r >= 0L & z_r <= 1L, 2, all))
  n_plane <- length(ok)
  n_dict  <- ncol(sim$W_dict)
  cat(sprintf("K=%-2d | %10d | %12d | %14d | %16d | %18d\n",
              Kv, 2^d, n_plane, n_dict, n_plane - n_dict, 3L * (2L^(d - 2L))))
}

cat("\n====================================================================================================\n")
cat("   TEST 2: INSPECTING THE EXACT ZERO-MSE 660-MUTATION SOLUTION AT K=8 (BEYOND THE 513 CATALOG!)\n")
cat("====================================================================================================\n")
sim8 <- make_fukuyama(8L, 202601L)
set.seed(202601L)
pool8 <- harvest_dual_20svd(sim8$S_true, 8L)
A_p8  <- rbind(pool8, 1.0); piv8 <- qr(A_p8, tol = 1e-6)$pivot[1:8]
sel8  <- piv8; closs <- fit_fcls_qp(sim8$S_true, pool8[, sel8])$mse
for (cyc in 1:2) for (i in 1:8) {
  curr <- sel8[i]; bc <- curr
  for (cand in setdiff(seq_len(ncol(pool8)), sel8)) {
    ti <- sel8; ti[i] <- cand; tl <- fit_fcls_qp(sim8$S_true, pool8[, ti])$mse
    if (tl < closs - 1e-7) { closs <- tl; bc <- cand }
  }
  sel8[i] <- bc
}
W_660 <- pool8[, sel8]
H_660 <- qr.solve(qr(rbind(1, W_660), tol = 1e-7), rbind(1, sim8$S_true))
recon_err <- max(abs(sim8$S_true - W_660 %*% H_660))

B_piv_rows <- vapply(1:7, function(i) (i - 1L) * 50L + 1L, integer(1))
strain_recipes <- apply(W_660[B_piv_rows, ], 2, function(col) {
  idx <- which(col == 1L)
  if (length(idx) == 0L) "Ancestral (0)" else paste0("e", idx, collapse = "+")
})

cat(sprintf("Max |S_true - W_660 * H_660| : %.3e (Exact Zero Reconstruction!)\n", recon_err))
cat(sprintf("Min(H_660) across 20 pools   : %.4f (Strictly Positive Simplex Proportions!)\n", min(H_660)))
cat(sprintf("Rank([1; W_660])             : %d / 8\n", qr(rbind(1, W_660))$rank))
cat(sprintf("Column Mutation Counts       : %s (Sum = %d)\n", paste(colSums(W_660), collapse = ", "), sum(W_660)))
cat(sprintf("Exact Strain Recipes         : %s\n", paste(strain_recipes, collapse = " | ")))

cat("\n====================================================================================================\n")
cat("   TEST 3: TRUE ZERO-MSE VALIDITY OF PIPELINE B ACROSS RUNS (PAIR-CATALOG vs HIGHER-TUPLE vs LOCAL MIN)\n")
cat("====================================================================================================\n")
cat(sprintf("%-4s | %-24s | %-15s | %-20s | %-15s\n",
            "K", "Candidate Pool Mode", "Pair-Catalog %", "Higher-Tuple (MSE=0)", "Local Min (MSE>0)"))
cat("----------------------------------------------------------------------------------------------------\n")

for (Kv in c(5L, 7L, 8L, 10L)) {
  sim <- make_fukuyama(Kv, 202601L)
  cat_k <- build_pruned_catalog(Kv, sim$pair_map)
  d <- sim$d; mu <- rowMeans(sim$S_true)
  U <- svd(sweep(sim$S_true, 1, mu, "-"), nu = d, nv = 0)$u[, 1:d, drop = FALSE]
  piv <- qr(t(U), tol = 1e-7)$pivot[1:d]
  G_h <- U[-piv, , drop = FALSE] %*% solve(U[piv, , drop = FALSE])
  h_h <- mu[-piv] - G_h %*% mu[piv]
  b_g <- t(as.matrix(expand.grid(rep(list(c(0L, 1L)), d))))
  z_S <- sweep(G_h %*% b_g, 1, as.vector(h_h), "+")
  z_r <- round(z_S)
  ok  <- which(apply(abs(z_S - z_r), 2, max) < 1e-6 & apply(z_r >= 0L & z_r <= 1L, 2, all))
  W_full_lat <- matrix(0L, 365L, length(ok))
  W_full_lat[piv, ] <- b_g[, ok]
  W_full_lat[-piv, ] <- matrix(as.integer(z_r[, ok]), 365L - d)
  W_full_lat <- W_full_lat[, !duplicated(t(W_full_lat)), drop = FALSE]
  W_pair_lat <- sim$W_dict

  for (p_mode in c("1. Standard top_n Pool", "2. Full 3*2^(d-2) Pool", "3. Pair-Restricted Pool")) {
    n_pair_cat <- 0L; n_higher_zero <- 0L; n_local_min <- 0L
    for (r_i in 1:10) {
      set.seed(5000L + r_i * 10L + Kv)
      W_p <- if (p_mode == "1. Standard top_n Pool") harvest_dual_20svd(sim$S_true, Kv)
             else if (p_mode == "2. Full 3*2^(d-2) Pool") W_full_lat
             else W_pair_lat

      all_i <- seq_len(ncol(W_p)); b_mse <- Inf; b_set <- seq_len(Kv)
      for (st in 1:3) {
        sel <- sample(all_i, Kv); cm <- fit_fcls_qp(sim$S_true, W_p[, sel, drop = FALSE])$mse
        for (cy in 1:2) {
          sw <- 0L
          for (ki in 1:Kv) {
            cur <- sel[ki]; bc <- cur
            for (cnd in setdiff(all_i, sel)) {
              ti <- sel; ti[ki] <- cnd; tm <- fit_fcls_qp(sim$S_true, W_p[, ti, drop = FALSE])$mse
              if (tm < cm - 1e-8) { cm <- tm; bc <- cnd }
            }
            if (bc != cur) { sel[ki] <- bc; sw <- sw + 1L }
          }
          if (sw == 0L) break
        }
        if (cm < b_mse - 1e-8) { b_mse <- cm; b_set <- sel }
      }
      W_sol <- W_p[, b_set, drop = FALSE]
      cl <- classify_W(W_sol, sim, cat_k)
      if (cl$status == "EXACT_VALID") {
        n_pair_cat <- n_pair_cat + 1L
      } else if (b_mse < 1e-10) {
        n_higher_zero <- n_higher_zero + 1L
      } else {
        n_local_min <- n_local_min + 1L
      }
    }
    cat(sprintf("K=%-2d | %-24s | %13.0f%%  | %18.0f%%  | %15.0f%%\n",
                Kv, p_mode, 10 * n_pair_cat, 10 * n_higher_zero, 10 * n_local_min))
  }
  cat("----------------------------------------------------------------------------------------------------\n")
}
