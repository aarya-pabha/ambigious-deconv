#!/usr/bin/env Rscript
# ==============================================================================
# File: summarize_80_fukuyama_all_pipelines.R
# Purpose: Compares all 5 pipeline generations head-to-head across the 80
#          ambiguous Fukuyama tasks (20 seeds per K in {5, 7, 8, 10}):
#          - Part A: Clean Data (Exact MSE=0 %, Mean Volume, Cycle Distribution)
#          - Part B: Noisy 100x Data (Step 3 vs Step 4 SNPs, Polish Delta, GT Dist)
# ==============================================================================

files <- list.files(if (dir.exists("results/fukuyama_80_all_pipelines")) "results/fukuyama_80_all_pipelines" else "../results/fukuyama_80_all_pipelines", pattern = "^res_task_.*\\.rds$", full.names = TRUE)
if (length(files) == 0) stop("No results found in ../results/fukuyama_80_all_pipelines")
recs <- lapply(files, readRDS)
recs <- recs[order(vapply(recs, `[[`, integer(1), "task_id"))]

fmt_fams <- function(f_vec) {
  tb <- table(f_vec)
  paste(sprintf("%s:%d", names(tb), as.integer(tb)), collapse = ", ")
}

cat("\n========================================================================================================================================\n")
cat(sprintf("   HEAD-TO-HEAD BENCHMARK OF ALL 5 PIPELINES ON THE 80 FUKUYAMA TASKS (COMPLETED: %d / 80 TASKS)\n", length(recs)))
cat("========================================================================================================================================\n")

for (Kv in c(5L, 7L, 8L, 10L)) {
  sub <- Filter(function(r) r$K == Kv, recs)
  n_tasks <- length(sub)
  m_k_size <- switch(as.character(Kv), "5" = 3, "7" = 71, "8" = 513, "10" = 32554)

  cat(sprintf("\n----------------------------------------------------------------------------------------------------------------------------------------\n"))
  cat(sprintf("   K = %2d (Ambiguity Manifold Size = %d Solutions | 20 Independent Seeds Evaluated)\n", Kv, m_k_size))
  cat(sprintf("----------------------------------------------------------------------------------------------------------------------------------------\n"))

  # Table 1: CLEAN DATA
  cat("PART 1: CLEAN DATA (DOES IT LAND ON VALID CYCLES OR GET TRAPPED?)\n")
  cat(sprintf("  %-22s | %-12s | %-14s | %-s\n", "Pipeline Version", "Exact MSE=0%", "Mean |det(A)|", "Cycle Family Distribution Across 20 Seeds"))
  cat(sprintf("  %s\n", paste(rep("-", 100), collapse = "")))

  pipe_keys <- c("c_pa", "c_pb", "c_v3", "c_v4_std", "c_v4_sub")
  pipe_names <- c("1. PipeA (sol4_v7)   ",
                  "2. PipeB (v2_Current)",
                  "3. Upgraded v3       ",
                  "4. Upgraded v4_Std   ",
                  "5. Upgraded v4_SubReg")

  for (pk_i in seq_along(pipe_keys)) {
    pk <- pipe_keys[pk_i]; pname <- pipe_names[pk_i]
    ex_rate <- 100 * mean(vapply(sub, function(r) r[[pk]]$exact, logical(1)))
    vols    <- mean(vapply(sub, function(r) r[[pk]]$vol, numeric(1)))
    f_str   <- fmt_fams(vapply(sub, function(r) r[[pk]]$family, character(1)))
    cat(sprintf("  %-22s | %10.0f%%  | %12.2fx | %s\n", pname, ex_rate, vols, f_str))
  }

  # Table 2: NOISY DATA (100x READS)
  cat("\nPART 2: NOISY DATA (100x READS: DOES THE SNR GATE STOP NOISE OVERFITTING?)\n")
  cat(sprintf("  %-22s | %-14s | %-14s | %-14s | %-14s\n",
              "Pipeline Version", "Step 3 (M_K)", "Step 4 Polish", "Polish Delta", "Dist to GT (0-cyc)"))
  cat(sprintf("  %s\n", paste(rep("-", 100), collapse = "")))

  n_pipe_keys <- c("n_pa", "n_pb", "n_v3", "n_v4_std", "n_v4_sub")
  for (pk_i in seq_along(n_pipe_keys)) {
    pk <- n_pipe_keys[pk_i]; pname <- pipe_names[pk_i]
    s3_mk <- mean(vapply(sub, function(r) r[[pk]]$s3_dist_mk, numeric(1)))
    s4_mk <- mean(vapply(sub, function(r) r[[pk]]$s4_dist_mk, numeric(1)))
    p_del <- mean(vapply(sub, function(r) r[[pk]]$polish_delta, numeric(1)))
    d_gt  <- mean(vapply(sub, function(r) r[[pk]]$s4_dist_gt, numeric(1)))
    cat(sprintf("  %-22s | %10.2f SNPs | %10.2f SNPs | %+10.2f SNPs | %10.2f SNPs\n",
                pname, s3_mk, s4_mk, p_del, d_gt))
  }
}
cat("========================================================================================================================================\n\n")
