#!/usr/bin/env Rscript
# ==============================================================================
# File: verify_volume_and_unsupervised_pair_filter.R
# Purpose: (1) Verifies the exact Simplex Volume Law |det(A)| = 2^c across all
#          cycle families at K=10, and (2) verifies 100% blind Pair-Catalog
#          sampling in Pipeline A & Pipeline B using unsupervised support-depth <= 2
#          filtering (no oracle W_dict) + zero Tikhonov damping on clean subspaces.
# Usage:   Rscript verify_volume_and_unsupervised_pair_filter.R
# ==============================================================================

suppressPackageStartupMessages({
  library(quadprog)
  library(MASS)
})

lines <- readLines("verify_exact_vs_beyond_catalog.R")
eval(parse(text = lines[1:(grep("# 2. RUN TESTS", lines)[1] - 1)]))

# Unsupervised Support-Depth <= 2 Extractor (zero Tikhonov damping on clean flat, no oracle)
extract_unsupervised_pair_pool <- function(S_in, K) {
  m <- nrow(S_in); d <- K - 1L; mu <- rowMeans(S_in)
  U <- svd(sweep(S_in, 1, mu, "-"), nu = d, nv = 0)$u[, 1:d, drop = FALSE]
  piv <- qr(t(U), tol = 1e-7)$pivot[1:d]
  G_h <- U[-piv, , drop = FALSE] %*% solve(U[piv, , drop = FALSE])
  h_h <- mu[-piv] - G_h %*% mu[piv]
  b_g <- t(as.matrix(expand.grid(rep(list(c(0L, 1L)), d))))
  z_S <- sweep(G_h %*% b_g, 1, as.vector(h_h), "+")
  z_r <- round(z_S)
  ok  <- which(apply(abs(z_S - z_r), 2, max) < 1e-5 & apply(z_r >= 0L & z_r <= 1L, 2, all))
  W_lat <- matrix(0L, m, length(ok))
  W_lat[piv, ] <- b_g[, ok]; W_lat[-piv, ] <- matrix(as.integer(z_r[, ok]), m - d)
  W_lat <- W_lat[, !duplicated(t(W_lat)), drop = FALSE]

  # Unsupervised Support-Inclusion Depth: keep Ancestral (0), Singletons (1), and Pairs (2)
  anc_i <- which.min(colSums(W_lat)); W_anc <- W_lat[, anc_i, drop = FALSE]
  W_rst <- W_lat[, -anc_i, drop = FALSE]
  rel_r <- abs(sweep(W_rst, 1, as.vector(W_anc), "-"))
  gram  <- crossprod(rel_r); c_n <- diag(gram); n_r <- ncol(rel_r)
  is_sng <- rep(TRUE, n_r)
  for (i in seq_len(n_r)) for (j in seq_len(n_r)) {
    if (i != j && c_n[j] < c_n[i] && gram[i, j] == c_n[j]) { is_sng[i] <- FALSE; break }
  }
  sng_idx <- which(is_sng)
  # Count how many singletons are subsets of each candidate j: keep if <= 2
  sng_contained <- vapply(seq_len(n_r), function(j) {
    sum(vapply(sng_idx, function(s_i) gram[s_i, j] == c_n[s_i], logical(1)))
  }, integer(1))
  cbind(W_anc, W_rst[, sng_contained <= 2L, drop = FALSE])
}

cat("========================================================================================\n")
cat("   VERIFICATION 1: EXACT SIMPLEX VOLUME LAW |det(A)| = 2^c AT K = 10 (32,554 SOLUTIONS)\n")
cat("========================================================================================\n")
sim10 <- make_fukuyama(10L, 202601L)
cat10 <- build_pruned_catalog(10L, sim10$pair_map)
W_aug_gt <- rbind(1.0, sim10$W_true)
vol_gt   <- sqrt(det(crossprod(W_aug_gt)))

fams_u <- unique(vapply(cat10$parts, `[[`, character(1), "family"))
cat(sprintf("%-16s | %-12s | %-18s | %-18s\n", "Cycle Family", "Count in M_10", "Relative Vol |det(A)|", "Min(H_alt) Weight"))
cat("----------------------------------------------------------------------------------------\n")
for (fm in fams_u) {
  idx_f <- which(vapply(cat10$parts, `[[`, character(1), "family") == fm)
  W_alt <- sim10$W_dict[, cat10$parts[[idx_f[1L]]]$cols]
  W_aug <- rbind(1.0, W_alt)
  rel_v <- sqrt(det(crossprod(W_aug))) / vol_gt
  H_alt <- qr.solve(qr(W_aug, tol = 1e-7), rbind(1.0, sim10$S_true))
  cat(sprintf("%-16s | %12d | %17.4fx | %18.4f\n", fm, length(idx_f), rel_v, min(H_alt)))
}

cat("\n========================================================================================\n")
cat("   VERIFICATION 2: UNSUPERVISED SUPPORT-DEPTH <= 2 POOL (NO ORACLE) ACROSS 20 SEEDS\n")
cat("========================================================================================\n")
cat(sprintf("%-4s | %-14s | %-16s | %-14s | %-s\n",
            "K", "UnsupPoolSize", "Exact Pair-Cat%", "LocalMin(MSE>0)", "Cycle Families Sampled Across 20 Seeds"))
cat("----------------------------------------------------------------------------------------\n")

for (Kv in c(5L, 7L, 8L, 10L)) {
  sim0 <- make_fukuyama(Kv, 202601L); cat_k <- build_pruned_catalog(Kv, sim0$pair_map)
  fams <- character(20); p_sizes <- integer(20); ex_ok <- logical(20)
  for (r_i in 1:20) {
    sim <- make_fukuyama(Kv, 202600L + r_i)
    set.seed(202600L + r_i)
    W_p <- extract_unsupervised_pair_pool(sim$S_true, Kv)
    p_sizes[r_i] <- ncol(W_p)

    # Unbiased 3-start random coordinate exchange on unsupervised pair pool
    all_i <- seq_len(ncol(W_p)); bm <- Inf; bs <- seq_len(Kv)
    for (st in 1:3) {
      sel <- sample(all_i, Kv); cm <- fit_fcls_qp(sim$S_true, W_p[, sel])$mse
      for (cy in 1:2) {
        sw <- 0L
        for (ki in 1:Kv) {
          cur <- sel[ki]; bc <- cur
          for (cnd in setdiff(all_i, sel)) {
            ti <- sel; ti[ki] <- cnd; tm <- fit_fcls_qp(sim$S_true, W_p[, ti])$mse
            if (tm < cm - 1e-8) { cm <- tm; bc <- cnd }
          }
          if (bc != cur) { sel[ki] <- bc; sw <- sw + 1L }
        }
        if (sw == 0L) break
      }
      if (cm < bm - 1e-8) { bm <- cm; bs <- sel }
    }
    cl <- classify_W(W_p[, bs], sim, cat_k)
    ex_ok[r_i] <- (cl$status == "EXACT_VALID") && (bm < 1e-10)
    fams[r_i]  <- cl$family
  }
  tb <- table(fams)
  f_str <- paste(sprintf("%s:%d", names(tb), as.integer(tb)), collapse = ", ")
  cat(sprintf("K=%-2d | %8.1f cols | %14.0f%%  | %13.0f%%  | %s\n",
              Kv, mean(p_sizes), 100 * mean(ex_ok), 100 * mean(!ex_ok), f_str))
}
cat("========================================================================================\n")
