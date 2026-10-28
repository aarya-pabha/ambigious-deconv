# Where the sparse candidates come from

The strongest demonstrated mechanism in the noisy ambiguous experiments is the use of distance to an inaccurately estimated affine plane to rank, modify, and sometimes penalize binary candidates. For binary vectors with low observed allele frequencies, that distance contains a positive weighted mutation-count term. With the correct plane, the projection term cancels it for every exact equivalent vertex. With the noisy estimated plane, that cancellation fails.

This concerns sparsity of the proposed and selected columns, and hence their distances to the catalog. It does not mean the algorithm recovered the labeled truth, nor that its output has fewer mutations than the truth. For example, task 41's labeled truth has 380 total ones, while the three noisy fitted predictions have 493–575.

## The equation behind the effect

The generators and snapping routines score a binary column w using

\[
R(w)=\|(I-UU^T)(w-\mu)\|^2
=\|\mu\|^2+\sum_i(1-2\mu_i)w_i-\|U^T(w-\mu)\|^2.
\]

This is an exact identity for orthonormal U and binary w. At K=8, the private-SNP frequencies are near 1/8, so their mutation-count weights are about 0.75; the shared-SNP weights are about 0.5. A mutation must be compensated by sufficient captured projection energy to avoid increasing R.

On clean data, all equivalent vertices lie on the true affine plane and have R=0 up to numerical error. In the three noisy tasks examined here, mean squared subspace overlap, `sum(crossprod(U_true, U_est)^2)/(K-1)`, is only **8.1–10.2%**. The estimated plane therefore fails to represent most of the true subspace. Residual ranking then strongly tracks mutation count: correlations for v4 Slawski candidate sets are **0.940–0.985** across these tasks' weighted and unweighted subspaces.

The source is `evaluation/benchmark_all_pipelines_on_80_fukuyama.R`: generator ranking at lines 322–367 and 449–464; snapping at lines 138–176; SubReg scoring at lines 508–530. The separate pipeline files contain these same mechanisms, though their orchestration is not identical to the benchmark.

## What the measured stages do

Instrumented the benchmark on tasks 41–43: K=8, seeds 202601–202603. Logging consumes no random numbers and leaves the algorithm unchanged. All nine noisy output matrices match the preceding tie-audit reruns exactly, including their column order.

Representative task 41, mutations per candidate column:

| Stage | Before | After | Scope |
|---|---:|---:|---|
| v3 Slawski residual ranking | 97.2 | 59.1 | Mean of the 20 extraction-call means |
| v4 Slawski residual ranking | 82.5 | 44.2 | Mean of the 20 extraction-call means |
| v4 Euclidean-ray residual ranking | 83.2 | 77.7 | Mean of the 20 extraction-call means |
| v4 arcsine-ray residual ranking | 17.9 | 13.5 | Baseline extraction |
| v3 snapping | 67.2 | 57.4 | Same 834 columns before/after snapping |
| v4 snapping | 60.2 | 50.4 | Same 830 columns before/after snapping |

These rows describe different candidate subsets, so they should not be read as one sequence with matching intermediate totals. v3 retains both raw and snapped candidates; v4 retains only snapped candidates and removes duplicates. The final deduplicated v4 pool in this task averages 51.6 mutations per column.

Across all three tasks, **4,971 of 4,977 candidate-column instances lost mutations during noisy snapping; none gained mutations**. Mean losses are 9.82–9.89 mutations per column. Snapping is limited to ten iterations. On clean data, mean count changes are small (for task 41, -0.267 for v3 and -0.125 for v4), consistent with the correctly estimated plane behaving differently.

The pruning is strong enough to exclude entire catalog alternatives. In task 41, the largest selectable column contains 95 mutations in v3 and 93 in v4. Every pair vertex in the K=8 catalog requires at least 100. Therefore **no nontrivial pair-cycle catalog solution can be assembled from these pools**. In all three tasks, the v4 pool maximum is below 100 (93, 86, 92). The v3 maxima in tasks 42 and 43 are 101 and 106, so the count bound alone does not establish exclusion there.

## Other mechanisms and their importance

**Candidate proposals already depend on the low observed frequencies.** Euclidean rays are centered at mu, have a fixed radius sqrt(m)/2, and are thresholded at 0.5. Most row means are about 0.125. The positive and negative rays are symmetric around this low mean, not around the binary threshold. The arcsine generator uses the same numerical radius in angle space; its threshold corresponds to pi/2, well above the transformed mean. Its proposals are particularly sparse (about 18 mutations per column before residual filtering in these K=8 runs). This is measured proposal behavior; it does not imply every individual proposal is sparse or that changing the radius alone would solve the problem.

**Damped Slawski extrapolation is an additional plausible contributor.** Its non-pivot coordinates have the form `mu_S + G_damped %*% (b - mu_R)`, followed by thresholding at 0.5. Ridge damping suppresses extrapolation around a low mean. The diagnostic did not ablate damping, so its contribution cannot be separated from pivot choice and subsequent residual ranking. v3 applies damping even on clean data; v4 disables it on clean subspaces but retains it for noisy data.

**SubReg applies the residual preference again during coordinate exchange.** Its score adds `0.03 * sum(R(w_k)) / (m*K)`. Evaluating this penalty term on the known equivalent catalog gives:

| Task | Penalty on labeled truth | Smallest penalty on any alternative |
|---|---:|---:|
| 41 | 0.002971 | 0.003968 |
| 42 | 0.003061 | 0.004060 |
| 43 | 0.003054 | 0.004047 |

On the clean plane, these penalties are around 1e-31. On the noisy plane, the penalty term favors the labeled truth over every alternative in the catalog. These are isolated penalty values, not full optimization scores; the catalog sets need not be available in the noisy candidate pool. SubReg uses this penalty during swaps, then compares the resulting restarts by unpenalized observed-data MSE.

**Sparse starts remain, but they are not a sufficient explanation.** v3's third start and v4's fourth start select from the 60 columns with the fewest mutations. Setting `lambda_pars=0` does not remove them. However, they win only three of nine noisy method/task combinations in this diagnostic. In every one of those three wins, the final result is *denser* than the best result from the other starts:

| Task/method | Winning sparse-start result: total ones | Best non-sparse-start result: total ones |
|---|---:|---:|
| 42 / v3 | 569 | 524 |
| 42 / v4_Std | 538 | 498 |
| 43 / v4_Std | 516 | 494 |

Because starts are run independently, the logged minimum-MSE non-sparse start is the result of omitting the sparse start with the pool and other starts held fixed. Removing that start would not eliminate the observed preference for low-mutation candidates, and would actually reduce the final mutation count in these three cases.

**Candidate order influences the nominal QR start.** The code uses base R `qr(..., tol=...)` with its default `LAPACK=FALSE`. That routine does not perform largest-norm-first pivoting; it largely preserves input order and moves numerically dependent columns aside. The pool starts with residual-sorted Slawski candidates, so the QR start inherits their ordering. Calling this initializer "QR MaxVol" in the README overstates what the code does. A local check with `diag(c(1,100,10))` yields default pivots `1,2,3`, whereas `LAPACK=TRUE` yields `2,3,1`.

**Local optimization preserves path dependence.** Coordinate exchange runs only two sweeps and requires an improvement of more than 1e-7. The benchmark stops early below a 1e-8 score. These rules do not explore the set of equivalent solutions or move freely across equal-loss solutions. Final strict-less-than comparisons also allow numerical residual differences to choose between nearly exact fits. None of these rules implies uniform selection across equally plausible factorizations.

**Polishing is not responsible for the noisy cases tested.** The singular-value-gap gate skips it in these tasks. PipeA uses related residual ranking and snapping mechanisms as well, but its diversity constraints and bootstrap search differ; it was reviewed in the source, not rerun in this diagnostic. PipeB's deprecated explicit mutation penalty is not needed for any of the measured v3/v4 effects.

## Interpretation and limits

The chain supported by these measurements is: little variation between mixture proportions makes the noisy plane inaccurate; residual ranking and snapping then prefer fewer mutations and exclude dense catalog vertices; the selector operates within that restricted pool. SubReg adds a further preference through the same residual. The resulting sparse predictions can be uniquely closest in Hamming distance to the sparse labeled truth, including after SNP locations are shuffled as shown by the previous tie audit.

This explains a substantial algorithmic source of sparse candidates in the noisy experiments. It does not establish that every selection bias has been isolated, that the same magnitude applies to other mixture regimes, or that a particular algorithm samples the clean equivalence class uniformly. A useful future ablation would target residual-based pruning and snapping with the candidate pool and other initialization choices controlled; simply removing sparse starts is insufficient in the runs measured here.

## Reproduction

```sh
Rscript evaluation/diagnose_sparsity.R 41 42 43
```

Run from the repository root with `quadprog`, `MASS`, and `clue` installed. The script requires the corresponding saved predictions from `evaluation/audit_equivalence_ties.R` to verify that instrumentation leaves results unchanged. Each task's RDS retains the source fingerprint and R session information. Per-task CSVs and combined `generators.csv`, `snapping.csv`, `starts.csv`, `projection.csv`, and `pools.csv` contain the measurements described above.
