#!/usr/bin/env Rscript
# ==============================================================================
# File: summarize_v4_200.R
# Purpose: Unaggregated cell-by-cell report (K x regime x n, averaged ONLY over
#          reps) comparing PipeA, PipeB_v2, v3, v4_Std, and v4_SubReg.
# ==============================================================================

files <- list.files(if (dir.exists("results/unified_v4_benchmark")) "results/unified_v4_benchmark" else "../results/unified_v4_benchmark", pattern = "^v4_task_.*\\.rds$", full.names = TRUE)
if (length(files) == 0) stop("No v4_task RDS files found.")
df <- do.call(rbind, lapply(files, readRDS))

agg <- aggregate(cbind(snr_gap, pa_s4_err, pa_s4_rec, pb_v2_err, pb_v2_rec,
                       v3_err, v3_rec, v4_std_s3_err, v4_std_s4_err, v4_std_s4_rec,
                       v4_sub_s3_err, v4_sub_s4_err, v4_sub_s4_rec) ~ K + regime + n,
                 data = df, FUN = function(x) mean(x, na.rm = TRUE), na.action = na.pass)
agg <- agg[order(agg$K, agg$regime, agg$n), ]

cat("\n========================================================================================================================================================================\n")
cat(sprintf("   UNAGGREGATED 200-TASK BENCHMARK: PIPE_A vs PIPE_B(v2) vs V3 vs V4_STD vs V4_SUBREG (COMPLETED: %d / 200 TASKS)\n", nrow(df)))
cat("========================================================================================================================================================================\n")
cat(sprintf("%-4s | %-20s | %-4s | %-6s | %-15s | %-15s | %-15s | %-15s | %-15s\n",
            "K", "Regime", "n", "SV_Gap", "PipeA (sol4_v7)", "PipeB (v2_Curr)", "v3 (UnsnapPool)", "v4_Std (4Start)", "v4_SubReg (SVD)"))
cat("------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n")
for (i in seq_len(nrow(agg))) {
  r <- agg[i, ]
  cat(sprintf("K=%-2d | %-20s | n=%-2d | %5.2fx | %6.2f (%5.1f%%) | %6.2f (%5.1f%%) | %6.2f (%5.1f%%) | %6.2f (%5.1f%%) | %6.2f (%5.1f%%)\n",
              r$K, substr(r$regime, 1, 20), r$n, r$snr_gap,
              r$pa_s4_err, r$pa_s4_rec, r$pb_v2_err, r$pb_v2_rec,
              r$v3_err, r$v3_rec,
              r$v4_std_s4_err, r$v4_std_s4_rec,
              r$v4_sub_s4_err, r$v4_sub_s4_rec))
}
cat("========================================================================================================================================================================\n\n")
