#!/usr/bin/env Rscript
base_dir <- if (dir.exists("data")) "data" else "../data"
tasks_dir <- file.path(base_dir, "fukuyama_80_tasks")
dir.create(tasks_dir, recursive = TRUE, showWarnings = FALSE)

K_map    <- c(rep(5L, 20), rep(7L, 20), rep(8L, 20), rep(10L, 20))
seed_map <- rep(202601L:202620L, 4)
manifest_rows <- list()

for (task_id in 1:80) {
  K <- K_map[task_id]; seed <- seed_map[task_id]; rep_i <- ((task_id - 1L) %% 20L) + 1L
  set.seed(seed)
  d <- K - 1L; m <- 365L; n <- 20L
  p <- 350L %/% d; s <- m - d * p

  B_p <- matrix(0L, d * p, d)
  for (i in 1:d) B_p[((i - 1L) * p + 1L):(i * p), i] <- 1L
  B_s <- matrix(0L, s, d); B_s[, 1:2] <- 1L
  B   <- rbind(B_p, B_s)
  W_true <- cbind(0L, B)
  colnames(W_true) <- paste0("Strain_", seq_len(K))

  G <- matrix(runif(K * n, 0.9, 1.1), K, n)
  H_true <- sweep(G, 2, colSums(G), "/")
  colnames(H_true) <- paste0("Sample_", seq_len(n))

  S_true <- W_true %*% H_true
  X_counts <- matrix(rbinom(m * n, 100L, as.vector(S_true)), m, n)
  colnames(X_counts) <- colnames(H_true)
  P_obs <- X_counts / 100.0

  t_folder <- sprintf("task_%02d_K%02d_rep%02d", task_id, K, rep_i)
  out_p <- file.path(tasks_dir, t_folder)
  dir.create(out_p, recursive = TRUE, showWarnings = FALSE)

  write.csv(W_true,   file.path(out_p, "W_true.csv"),   row.names = FALSE)
  write.csv(H_true,   file.path(out_p, "H_true.csv"),   row.names = FALSE)
  write.csv(S_true,   file.path(out_p, "S_true.csv"),   row.names = FALSE)
  write.csv(X_counts, file.path(out_p, "X_counts.csv"), row.names = FALSE)
  write.csv(P_obs,    file.path(out_p, "P_obs.csv"),    row.names = FALSE)

  m_k <- switch(as.character(K), "5" = 3L, "7" = 71L, "8" = 513L, "10" = 32554L)
  manifest_rows[[task_id]] <- data.frame(
    task_id = task_id, K = K, d = d, private_snps_p = p, shared_snps_s = s,
    samples_n = n, seed = seed, rep = rep_i, certified_manifold_size = m_k,
    task_dir = file.path("data/fukuyama_80_tasks", t_folder),
    stringsAsFactors = FALSE
  )
}

manifest_df <- do.call(rbind, manifest_rows)
write.csv(manifest_df, file.path(base_dir, "fukuyama_80_manifest.csv"), row.names = FALSE)
cat(sprintf("Generated all 80 Fukuyama tasks and manifest in %s\n", base_dir))
