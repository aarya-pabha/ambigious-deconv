# Genomic Haplotype Deconvolution: Audit & Unification of 5 Pipeline Architectures on the Fukuyama Equivalence Manifold

This repository provides a self-contained codebase, pre-generated 80-task ambiguity benchmark datasets ($K \in \{5, 7, 8, 10\}$), **5 deconvolution pipeline implementations**, evaluation scripts, and pre-computed results tables across:
1. **The 80 Ambiguous Fukuyama Manifold Tasks** ($K \in \{5, 7, 8, 10\} \times 20\text{ independent seeds}$, evaluated on both Clean expected allele frequencies $S_{\text{true}}$ and Noisy $100\times$ Binomial read counts $X$).
2. **The 200-Task General Benchmark Suite** ($K \in \{5, 7, 10\}$ across Balanced, Imbalanced, and Flat Dirichlet abundance regimes for $n \in \{15, 20\}$ samples, plus $K=8$ Fukuyama tasks).

---

## 0. Two Modes of Operation in This Repository (Read First)

When evaluating deconvolution on the ambiguous Fukuyama tasks, this repository distinguishes between two fundamentally different algorithmic goals:

* **Mode A — Complete Equivalence-Class Enumeration (The Stage 4 Pruned Circuit Walker):**
  Given clean allele frequencies $S_{\text{true}}$ (or block-denoised frequencies at high sequencing depth), the **Pruned Circuit Walker** (`pipelines/02_pipeB_unified_15plus_v2.R` and `evaluation/verify_exact_vs_beyond_catalog.R`) extracts the support-minimal singleton basis and traverses all vertex-disjoint odd cycles ($\ell \le d - 3$ pruning) to output **100% of the certified equivalence class on a single task**:
  * **$K = 5$:** Recovers all **$3 / 3$** solutions in `0.08s`
  * **$K = 7$:** Recovers all **$71 / 71$** solutions in `0.09s`
  * **$K = 8$:** Recovers all **$513 / 513$** solutions in `0.25s`
  * **$K = 10$:** Recovers all **$32,554 / 32,554$** solutions in `1.36s` (a $203\times$ speedup over unpruned backtracking)
* **Mode B — Single Point-Estimation Pipeline Audit (Pipelines 1 to 5):**
  Standard NMF/deconvolution pipelines (`PipeA: sol4_v7`, `PipeB: unified_15plus_v2`, `v3`, `v4_Std`, `v4_SubReg`) are designed to output **one single $K$-strain matrix $W_{\text{pred}} \in \{0,1\}^{365 \times K}$ per run**. Our primary audit question is:
  > *When run on the 80 ambiguous Fukuyama tasks (Clean and Noisy $100\times$), do these point-estimation pipelines always collapse onto the `0-cycle` Ground Truth, get trapped in non-exact local minima, or land across different valid odd-cycle solutions ($3, 3+3, 5, 5+3, 7, 9$-cycles)?*

---

## 1. Repository Directory Tree

All core helper libraries (`phase7/R/`) and pre-computed `.rds`/`.csv` results (`results/`) are included so every script runs out-of-the-box:

    fukuyama-deconvolution-audit/
    ├── README.md                                         # Complete guide to all 5 pipelines, columns, & results
    ├── install_dependencies.R                            # Installs CRAN packages (quadprog, MASS, Matrix, clue, lpSolve)
    ├── master_manifest.csv                               # 200-task general benchmark manifest (Tasks 1-180 Dirichlet, 181-200 Fukuyama)
    ├── phase7/R/                                         # Core plateNMF mathematical library sourced by Pipeline A
    │   ├── helper-functions.R                            # Affine projection, Anscombe transform, & error metrics
    │   ├── nmf-functions.R                               # Base subspace & ray-casting functions
    │   ├── nmf-functions_v2.R                            # FCLS QP solver (fit_fcls_qp) & weighted SVD subspace fitters
    │   └── simulation-functions.R                        # Bernoulli endmember & Dirichlet generators
    ├── data/
    │   ├── fukuyama_80_manifest.csv                      # Lookup table for all 80 ambiguous tasks (K=5,7,8,10 x 20 seeds)
    │   └── fukuyama_80_tasks/                            # 80 task folders (task_01_K05_rep01 .. task_80_K10_rep20)
    ├── pipelines/
    │   ├── 01a_pipeA_hunter_sol4_v2.R                    # Pipeline 1 (Stage 1-2): 3-Constraint Exponential SVD Hunter
    │   ├── 01b_pipeA_snapped_subset_v2.R                 # Pipeline 1 (Stage 3): 1/2-Bit Basis Snapping + Warm Exchange
    │   ├── 01c_pipeA_polish_v7.R                         # Pipeline 1 (Stage 4): Binomial Polish v7 (Pass A + Pass B)
    │   ├── 02_pipeB_unified_15plus_v2.R                  # Pipeline 2: Original Unified 15+ (20-SVD + Parsimony Exchange + Circuit Walker)
    │   ├── 03_upgraded_v3_pipeline.R                     # Pipeline 3: Unbiased 20-SVD (lambda=0) + Snapped Pool + SV_Gap-Gated Polish v7
    │   └── 04_upgraded_v4_std_and_subreg.R               # Pipelines 4 & 5: MaxVol 20-SVD + Pure Snapped Pool + 4-Start + Subspace-Reg Exchange
    ├── evaluation/
    │   ├── 00_generate_80_fukuyama_tasks.R               # Deterministically generates all 80 ambiguous tasks from scratch
    │   ├── benchmark_all_pipelines_on_80_fukuyama.R      # Self-contained runner for all 5 pipelines on any task (1-80)
    │   ├── audit_both_pipelines_fukuyama.R               # Multi-run trajectory & 3-constraint rejection audit on the 80 tasks
    │   ├── verify_exact_vs_beyond_catalog.R              # Verifies the 3 * 2^(d-2) binary lattice law & 660-mutation triple solution
    │   ├── verify_volume_and_unsupervised_pair_filter.R  # Verifies Simplex Volume Law |det(A)| = 2^c & unsupervised pair filter
    │   ├── summarize_80_fukuyama_all_pipelines.R         # Prints the 80-task Clean & Noisy comparison tables from RDS files
    │   └── summarize_v4_200.R                            # Prints the 200-task unaggregated cell-by-cell comparison table
    └── results/
        ├── MEETING_REPORT_TABLES.txt                     # Formatted side-by-side console tables for all experiments
        ├── table1_80_fukuyama_head_to_head.csv           # Clean CSV of all 5 pipelines on the 80 ambiguous tasks
        ├── table2_200_benchmark_unaggregated_cells.csv   # Clean CSV of all 5 pipelines across the 200-task cells
        ├── fukuyama_80_all_pipelines/                    # Pre-computed RDS outputs for the 80 ambiguous tasks
        └── unified_v4_benchmark/                         # Pre-computed RDS outputs for the 200 benchmark tasks

---

## 2. Exact Anatomy of the 80 Ambiguous Tasks (`data/fukuyama_80_tasks/`)

The script `evaluation/00_generate_80_fukuyama_tasks.R` generates 80 task directories (`task_01_K05_rep01` through `task_80_K10_rep20`), indexed in `data/fukuyama_80_manifest.csv` (`K` $\in \{5, 7, 8, 10\}$, `d = K - 1`, `seed` $\in \{202601\dots 202620\}$). Inside each task folder are **5 CSV matrices**:

1. **`W_true.csv` ($365 \times K$ integer binary matrix $\in \{0,1\}^{365 \times K}$):**
   * **Column 1 (`Strain_1`):** The all-zero ancestral strain $\mathbf{0}_{365}$.
   * **Columns $2 \dots K$ (`Strain_2` $\dots$ `Strain_K`):** The $d = K - 1$ derived strains $B \in \{0,1\}^{365 \times d}$.
   * **Rows $1 \dots d \cdot p$ (Private SNPs):** Each derived strain $i \in \{1 \dots d\}$ has $p = \lfloor 350 / d \rfloor$ private SNPs ($p = 87, 58, 50, 38$ for $K = 5, 7, 8, 10$) where only column $i + 1$ equals `1`.
   * **Rows $d \cdot p + 1 \dots 365$ (Shared SNPs):** The remaining $s = 365 - d \cdot p$ SNPs ($s = 17, 17, 15, 23$) are shared exclusively by Derived Strains 1 and 2 (`Strain_2` and `Strain_3` both equal `1`; all other columns equal `0`).
2. **`H_true.csv` ($K \times 20$ numeric matrix $\in \Delta^{K \times 20}$):**
   * True strain proportions across $n = 20$ samples. Raw entries are drawn i.i.d. from $g_{kj} \sim \operatorname{Uniform}(0.9, 1.1)$ and normalized per column ($H_{kj} = g_{kj} / \sum_{i=1}^K g_{ij}$), keeping every sample strictly inside a tight ball around the simplex centroid $(1/K)\mathbf{1}$.
3. **`S_true.csv` ($365 \times 20$ numeric matrix $\in [0,1]^{365 \times 20}$):**
   * Noiseless expected allele frequencies $S_{\text{true}} = W_{\text{true}} H_{\text{true}}$. Used for all **Clean Data** evaluations.
4. **`X_counts.csv` ($365 \times 20$ integer matrix $\in \{0, 1, \dots, 100\}^{365 \times 20}$):**
   * Noisy sequencing alternative-allele read counts at $100\times$ total read depth: $X_{mj} \sim \operatorname{Binomial}(100, S_{\text{true}, mj})$.
5. **`P_obs.csv` ($365 \times 20$ numeric matrix $\in [0,1]^{365 \times 20}$):**
   * Observed noisy allele frequencies $P_{\text{obs}} = X_{\text{counts}} / 100$. Used together with `X_counts.csv` for all **Noisy $100\times$ Data** evaluations.

---

## 3. Mathematical Laws Verified in `evaluation/`

1. **Certified Pair-Cycle Equivalence Class Size ($M_K$):**
   Let $E_{\text{valid}} = \binom{d}{2} - 1$ be the pairs of derived strains that do not share SNPs (all pairs $(a,b)$ except $(1,2)$). Replacing any subset of singletons $\{B e_{v_1}, \dots, B e_{v_\ell}\}$ ($\ell \ge 3$ odd) with the cycle of pairwise sums $\{B(e_{v_1}+e_{v_2}), \dots, B(e_{v_\ell}+e_{v_1})\}$ yields a valid binary matrix $W_{\text{alt}} = W_{\text{true}} A \in \{0,1\}^{365 \times K}$ and mixture matrix $H_{\text{alt}} = A^{-1} H_{\text{true}} > 0$ such that $W_{\text{alt}} H_{\text{alt}} = S_{\text{true}}$ exactly ($\text{MSE} = 0$ and identical Binomial log-likelihood $\Delta \log L = 0$):
   * **$K = 5$ ($d=4$):** **$M_5 = 3$** solutions (`1` 0-cycle, `2` 3-cycles)
   * **$K = 7$ ($d=6$):** **$M_7 = 71$** solutions (`1` 0-cycle, `10` 3-cycles, `12` 3+3-cycles, `48` 5-cycles)
   * **$K = 8$ ($d=7$):** **$M_8 = 513$** solutions (`1` 0-cycle, `25` 3-cycles, `95` 3+3-cycles, `192` 5-cycles, `200` 7-cycles)
   * **$K = 10$ ($d=9$):** **$M_{10} = 32,554$** solutions (`1` 0-cycle, `77` 3-cycles, `700` 3+3-cycles, `210` 3+3+3-cycles, `1,302` 5-cycles, `4,704` 5+3-cycles, `10,440` 7-cycles, `15,120` 9-cycles)
2. **The $3 \times 2^{d-2}$ Binary Lattice Law & Higher-Order Tuples (`verify_exact_vs_beyond_catalog.R`):**
   Out of all $2^d$ binary combinations of the $d$ pivot rows, exactly **$3 \times 2^{d-2}$** ($12, 48, 96, 384$ for $K=5,7,8,10$) produce binary vectors in $\{0,1\}^{365}$ lying on the clean SVD hyperplane (all subsets of derived strains that do not simultaneously activate both Strain 1 and Strain 2). Beyond the $1 + d + \binom{d}{2} - 1$ pair-catalog vectors ($10, 21, 28, 45$), the remaining vectors are higher-order tuples (triples $e_a+e_b+e_c$, quads, etc.) that also form valid exact $\text{MSE}=0, H>0$ factorizations (e.g., `W_660` at $K=8$ with 660 mutations, $\max|S_{\text{true}} - W_{660}H_{660}| = 9.99 \times 10^{-16}$, and $\min H_{660} = 0.0365$), proving empirically that $M_K$ is a strict lower bound. Unsupervised **Boolean support-inclusion depth $\le 2$** (`verify_volume_and_unsupervised_pair_filter.R`) isolates the exact pair dictionary ($10, 21, 28, 45$ columns) without oracle knowledge.
3. **The Simplex Volume Law $|\det(A)| = 2^c$ (`verify_volume_and_unsupervised_pair_filter.R`):**
   Every odd cycle has an incidence block with determinant $2$. Therefore, any solution $W_{\text{alt}} = W_{\text{true}} A$ composed of $c$ disjoint odd cycles has relative simplex volume $|\det(A)| = \sqrt{\det(\widetilde{W}_{\text{alt}}^T \widetilde{W}_{\text{alt}}) / \det(\widetilde{W}_{\text{true}}^T \widetilde{W}_{\text{true}})} = 2^c$:
   * **0-cycle Ground Truth ($c=0$):** $|\det(A)| = \mathbf{1.0000\times}$ (unique Minimum-Volume Enclosing Simplex)
   * **Single odd cycles ($c=1$: $3, 5, 7, 9$-cycles):** $|\det(A)| = \mathbf{2.0000\times}$
   * **Double odd cycles ($c=2$: $3+3, 5+3$-cycles):** $|\det(A)| = \mathbf{4.0000\times}$
   * **Triple odd cycles ($c=3$: $3+3+3$-cycles):** $|\det(A)| = \mathbf{8.0000\times}$

---

## 4. Detailed Code Walkthrough of the 5 Pipelines (`pipelines/`)

All 5 point-estimation pipelines reconstruct $W \in \{0,1\}^{365 \times K}$ through a 4-stage architecture:
* **Stage 1 (Subspace Estimation):** Computes sample mean $\mu = \frac{1}{n}\sum_{j=1}^n P_{\cdot j}$ and top $d = K - 1$ left singular vectors $U_b$ of $P - \mu \mathbf{1}^T$.
* **Stage 2 (Candidate Harvesting):** Generates a candidate pool `W_pool` of binary columns via Slawski affine subspace inversion (`z_S = G_hat %*% b_grid + h_hat` rounded to $\{0,1\}$) and random ray projections.
* **Stage 3A (Basis Snapping) & Stage 3B (Coordinate Exchange):** Optionally snaps each candidate column onto $(U_b, \mu)$ via 1-bit/2-bit residual descent, then selects $K$ columns via cyclic 1-column swap coordinate exchange using Fully Constrained Least Squares (`quadprog::solve.QP`, $H \ge 0, \mathbf{1}^T H = \mathbf{1}$).
* **Stage 4 (Bit-Flipping Binomial Polish):** Tests 1-bit flips (Pass A) and 2-bit pairwise swaps/co-mutations (Pass B) on raw counts $X$, subject to the subspace step-drift ceiling $d_{\text{res}} = 2(1 - 2W_{mk})r_{mk} + M_{mm} \le 0.85$ (where $M = I - U_b U_b^T$).

### Pipeline 1: `PipeA` (`01a_pipeA_hunter_sol4_v2.R` -> `01b_pipeA_snapped_subset_v2.R` -> `01c_pipeA_polish_v7.R`)
* **Stage 2 (`run_hunter_sol4`):** Draws parametric bootstrap count matrices $X^* \sim \operatorname{Binomial}(100, \hat{P})$ and applies random exponential locus weights $w_i \sim \operatorname{Exp}(1)$ to tilt the SVD hyperplane (`fit_weighted_affine_subspace`). On each tilt, `extract_slawski_vertex_set` selects $d$ well-conditioned pivots (`select_pivots_maxvol`), applies a `2%` Tikhonov ridge pseudo-inverse (`0.02 * norm(U_R, "F")^2`), enumerates all $2^d$ binary anchor combinations, and keeps the top $K$ columns closest to $U$. A proposed $K$-strain set $W^*$ is accepted into the Hunter catalog only if it passes **3 Constraints** against all previously saved sets:
  1. **Overlap Gate:** Fraction of columns matching an existing solution within $\le 1$ SNP must be $\le 60\%$ (`max_overlap_tol = 0.60`).
  2. **Distance Gate:** Mean bidirectional nearest-neighbor Hamming distance must be $\ge \max(15.0, 2.5K)$ SNPs.
  3. **MSE Ratchet:** FCLS reconstruction error must satisfy $\text{MSE}(W^*) \le 1.25 \times \text{best\_mse}$.
* **Stage 3A & 3B (`snap_pool_to_basis` + `select_best_k_subset_multistart`):** Snaps all unique Hunter columns via 1-bit and 2-bit residual descent onto $U_b$, combines raw and snapped columns (`cbind(W_pool, W_snapped)`), and runs 1-column coordinate exchange from the top 2 Hunter warm starts with **zero parsimony penalty ($\lambda_{\text{pars}} = 0$)**.
* **Stage 4 (`polish_bits_v7`):** Runs up to 3 epochs of **Pass A (1-bit flips)** and **Pass B (2-bit allele swaps $(1,0) \leftrightarrow (0,1)$ and joint co-mutations $(0,0) \leftrightarrow (1,1)$)**.
* *Implementation Note on `benchmark_all_pipelines_on_80_fukuyama.R`:* When auditing `PipeA` on **noiseless** $S_{\text{true}}$ (`clean_v1_unbiased` mode), `run_hunter_audit` tilts $S_{\text{true}}$ directly with $w \sim \operatorname{Exp}(1)$ using exact `solve()` (`damp = 0`) and random tie-breaking among zero-residual lattice vectors (`runif(..., 0, 1e-14)`) so that artificial Binomial resampling noise is not injected into clean data; on noisy $100\times$ data (`noisy_100x` mode), it runs the full Binomial bootstrap.

### Pipeline 2: `PipeB` (`02_pipeB_unified_15plus_v2.R`)
* **Stage 2 (`harvest_dual_20svd`, Point-Estimation Route):** Removes the 3-constraint rejection gate and harvests ~600-800 unique candidates across **20 SVDs**: `1` baseline SVD + `14` locus-weighted SVDs ($w_{\text{loc}} \sim \operatorname{Exp}(1)$) + `5` sample-weighted SVDs ($w_{\text{samp}} \sim 0.8\operatorname{Dir}(1) + 0.2/n$, which amplify rare/imbalanced strains), combining `extract_slawski` (`2%` Tikhonov damping), `extract_euclidean_rays`, and `extract_arcsine_rays`.
* **Stage 3 (`coord_exchange_robust_parsimony`):** Skips Stage 3A Basis Snapping and runs 3-start cyclic coordinate exchange (`Restart 1: QR MaxVol`, `Restart 2: Leverage Swap`, `Restart 3: Sparse MinVol`) with an explicit L1 parsimony penalty:
  $$\mathcal{L}(W) = \text{MSE}(P_{\text{obs}}, WH) + \lambda_{\text{pars}} \frac{\|W\|_1}{m K}, \qquad \lambda_{\text{pars}} = 5 \times 10^{-5}$$
* **Stage 4 (`polish_strains`):** Runs Pass A (1-bit flips only; omits Pass B 2-bit pairwise swaps) unconditionally.

### Pipeline 3: `Upgraded v3` (`03_upgraded_v3_pipeline.R`)
* **What Changed from `PipeB`:**
  1. **Zero Parsimony Penalty ($\lambda_{\text{pars}} = 0.0$):** Removes the L1 mutation penalty in Stage 3B so all valid cycle solutions have identical loss ($\text{MSE} = 0$).
  2. **Adds Stage 3A Basis Snapping (`snap_pool_to_basis`):** Snaps the 20-SVD candidate pool onto $(U_b, \mu_b)$ via 1-bit/2-bit descent (`cbind(W_pool, W_snapped)`).
  3. **Upgrades Stage 4 to `polish_bits_v7_full` (Pass A + Pass B) with the `SV_Gap > 1.12` Gate:** Computes the empirical singular value gap $\text{SV\_Gap} = \sigma_{K-1}(P_{\text{cent}}) / \sigma_K(P_{\text{cent}})$. When $\text{SV\_Gap} > 1.12$ (well-conditioned Dirichlet data), full Pass A + Pass B Polish runs; when $\text{SV\_Gap} \le 1.12$ (sub-BBP noisy $100\times$ Fukuyama data where $\text{SV\_Gap} \approx 1.02\times$), Stage 4 Polish is automatically bypassed so it cannot overfit shot noise.

### Pipeline 4: `Upgraded v4_Std` (`04_upgraded_v4_std_and_subreg.R`, `use_subspace_reg = FALSE`)
* **What Changed from `v3`:**
  1. **MaxVol Pivot Conditioning + Adaptive Tikhonov Damping (`extract_slawski_maxvol`):** Integrates `PipeA`'s `select_pivots_maxvol` into all 20 SVDs and sets Tikhonov damping to **`0.0` on clean subspaces** ($\sigma_{K-1}/\sigma_K > 10^4$, preventing shared-SNP shrinkage) and `0.02` on noisy subspaces.
  2. **Pure Snapped Pool (`snap_pool_only`):** Replaces raw unsnapped random rays with their snapped facets on $(U_b, \mu_b)$ instead of keeping unsnapped noise columns in the pool.
  3. **4-Start Hybrid Coordinate Exchange (`coord_exchange_v4`):** Combines **both** `PipeA`'s top 2 Slawski $K$-simplex warm starts and `PipeB`'s `QR MaxVol` and `Sparse` starts in a single tournament with $\lambda_{\text{pars}} = 0.0$.

### Pipeline 5: `Upgraded v4_SubReg` (`04_upgraded_v4_std_and_subreg.R`, `use_subspace_reg = TRUE`)
* **What Changed from `v4_Std`:**
  * During Stage 3B candidate swaps, evaluates FCLS against the **rank-$d$ SVD-denoised frequency matrix** $\hat{P}_{\text{svd}} = \Pi_{[0,1]}(\mu_b \mathbf{1}^T + U_b U_b^T(P_{\text{obs}} - \mu_b \mathbf{1}^T))$ plus an off-plane subspace residual penalty $0.03 \frac{\|(I - U_b U_b^T)(W_{\text{sub}} - \mu_b)\|_F^2}{m K}$.
  * On noisy Flat Dirichlet mixtures, scoring against $\hat{P}_{\text{svd}}$ prevents candidate swaps from overfitting the $n - K + 1$ orthogonal dimensions of binomial shot noise. On clean Fukuyama manifold tasks, $\hat{P}_{\text{svd}} = S_{\text{true}}$ and all manifold vertices lie exactly on $U_b$ (`cand_res_sq = 0`), preserving **zero bias ($\lambda_{\text{pars}} = 0$) and 100% exact manifold recovery**.

---

## 5. Glossary of CSV Columns, Metrics, & FAQ for Reading the Tables

### A. Column Glossary for `results/table1_80_fukuyama_head_to_head.csv`
* **`K`:** True number of strains ($K \in \{5, 7, 8, 10\}$; certified pair-cycle manifold size $M_K \in \{3, 71, 513, 32554\}$).
* **`Pipeline`:** Pipeline version (`1_PipeA_sol4_v7`, `2_PipeB_v2_Curr`, `3_Upgraded_v3`, `4_Upgraded_v4_Std`, `5_Upgraded_v4_SubReg`).
* **`Clean_Exact_Pct` (`Exact MSE=0%` in console output):** Percentage of the 20 independent seeds where the pipeline's output $W_{\text{pred}}$ on noiseless $S_{\text{true}}$ achieved **exact zero reconstruction error** ($\text{MSE}(S_{\text{true}}, W_{\text{pred}}\hat{H}) < 10^{-10}$ and $\hat{H} \ge 0$, i.e., **0.00 SNPs wrong across all $K$ strains**). If a run has even 1–2 SNPs wrong ($\text{MSE} > 10^{-10}$), it is counted as a failure (`FALSE`) in this column.
* **`Clean_Mean_Vol` (`Mean |det(A)|` in console output):** Mean relative simplex volume $|\det(A)| = \sqrt{\det(\widetilde{W}_{\text{pred}}^T \widetilde{W}_{\text{pred}}) / \det(\widetilde{W}_{\text{true}}^T \widetilde{W}_{\text{true}})}$ across the 20 seeds (`1.00x` = 0-cycle Ground Truth, `2.00x` = single odd cycle, `4.00x` = double odd cycle).
* **`Clean_Families` (`Cycle Family Distribution` in console output):** Classification of each of the 20 predictions into its **nearest** catalog family (`0-cycle`, `3-cycle`, `3+3-cycle`, `5-cycle`, `5+3-cycle`, `7-cycle`, `9-cycle`) by Hungarian matching distance, or `DICT_CHIMERA` if all $K$ columns belong to the pair dictionary $W_{\text{dict}}$ (`min_dist_to_dict == 0`) but do not form a valid disjoint odd-cycle partition (e.g., an unclosed path or duplicate singleton, which gives $\text{MSE} > 0$).
* **`Noisy_Step3_SNP` (`Step 3 (M_K)` in console output):** On noisy $100\times$ read counts $X$, the mean per-strain Hungarian Hamming distance (in SNPs out of 365) between the **Stage 3** prediction $W_{\text{s3}}$ and the **closest valid solution anywhere in the $M_K$ catalog**: $\min_{m \in \{1\dots M_K\}} d_{\text{Hungarian}}(W_{\text{s3}}, W^{(m)})$.
* **`Noisy_Step4_SNP` (`Step 4 Polish` in console output):** Same minimum distance to the $M_K$ catalog measured **after Stage 4 Binomial Polish** ($W_{\text{s4}}$).
* **`Noisy_Pol_Delta` (`Polish Delta` in console output):** `Noisy_Step4_SNP - Noisy_Step3_SNP`. A positive value (`+2.5` to `+16.1 SNPs`) means Stage 4 Polish **worsened** accuracy by flipping bits to overfit $100\times$ binomial shot noise; `+0.00 SNPs` means the `SV_Gap > 1.12` gate successfully bypassed Polish.
* **`Noisy_GT_SNP` (`Dist to GT (0-cyc)` in console output):** Mean per-strain Hungarian Hamming distance from $W_{\text{s4}}$ specifically to the **`0-cycle` Ground Truth $W_{\text{true}}$**.

### B. Column Glossary for `results/table2_200_benchmark_unaggregated_cells.csv`
* **`K`, `regime`, `n`:** Strain count $K \in \{5, 7, 10\}$, Dirichlet concentration vector $\alpha$ (`1,1,..` = Balanced, `5,1,..` = Imbalanced, `5,5,..` = Flat), and sample count $n \in \{15, 20\}$. Each row averages **only over the 10 replicate seeds** within that exact `(K, regime, n)` cell.
* **`snr_gap` (`SV_Gap`):** Empirical ratio of the last signal singular value to the first noise singular value, $\sigma_{K-1}(P_{\text{cent}}) / \sigma_K(P_{\text{cent}})$.
* **`pa_s4_err`, `pb_v2_err`, `v3_err`, `v4_std_s4_err`, `v4_sub_s4_err`:** Mean per-strain Hungarian Hamming error (in SNPs out of 365) against `W_true.csv` for `PipeA`, `PipeB`, `v3`, `v4_Std`, and `v4_SubReg`.
* **`pa_s4_rec`, `pb_v2_rec`, `v3_rec`, `v4_std_s4_rec`, `v4_sub_s4_rec`:** Exact strain recall percentage (percentage of the $K$ strains recovered with **0 SNP errors**, $d_{\text{Hamming}} = 0$).

### C. Key Questions & Explanations When Interpreting Table 1
1. **Why does `2. PipeB (v2_Curr)` at $K=8$ Clean show `Cycle Family = 0-cycle: 20` while `Exact MSE=0% = 30%`?**
   * Because `PipeB`'s parsimony penalty ($\lambda_{\text{pars}} = 5 \times 10^{-5}$) and `Restart 3 (Sparse)` push all 20 runs closest to the `0-cycle` (`0-cycle: 20`). However, `PipeB`'s `extract_slawski` uses a hardcoded `2%` Tikhonov ridge penalty (`0.02 * norm(U_R, "F")^2`), which shrinks continuous coordinates toward $\mu$ and flips 2 to 12 shared SNPs on 14 of the 20 seeds (70%), leaving $\text{MSE} \approx 10^{-5} > 0$. Setting adaptive damping to `0.0` on clean subspaces in `v4` brings `Exact MSE=0%` to **100%**.
2. **Why does `1. PipeA (sol4_v7)` at $K=8$ and $K=10$ Clean never output a `0-cycle` or `3-cycle`?**
   * Because Constraint 2 (`overlap_pct <= 0.60`) rejects any proposed $K$-strain set that shares $>60\%$ of its columns with an earlier accepted set. A `3-cycle` replaces 3 singletons and keeps $K - 3$ columns identical to the `0-cycle` ($5/8 = 62.5\% > 60\%$ at $K=8$; $7/10 = 70.0\% > 60\%$ at $K=10$), so `PipeA` systematically rejects `3-cycles` at $K \ge 8$ and forces convergence into long cycles (`5, 3+3, 5+3, 7, 9`-cycles).
3. **Why does `v3` at $K=10$ Clean sample `9-cycles` and `5+3-cycles` (`Mean |det(A)| = 2.35x`) while `v4` samples `0, 3, 5, 7`-cycles (`Mean |det(A)| = 1.95x`)?**
   * Both `v3` and `v4` achieve **100% exact $\text{MSE}=0$ recovery** on $K=10$ Clean. `v3` starts Coordinate Exchange with `Restart 1 = QR MaxVol`, which initializes with high-norm pair columns and lands in `9-cycle` and `5+3-cycle` basins. `v4` places `PipeA`'s Baseline Slawski simplex as `Start 1` (ahead of `QR MaxVol` in `Start 3`); on clean data, `Start 1` initializes near the singletons and reaches $\text{MSE} = 0$ on a `0, 3, 5, or 7`-cycle in 1–2 swaps, winning the strict `<` tie-break before `Start 3` runs.
4. **Why are `Noisy_Step4_SNP` (distance to the closest of all $M_K$ solutions) and `Noisy_GT_SNP` (distance to the `0-cycle`) identical on Noisy $100\times$ data?**
   * Under $100\times$ noise on $\operatorname{Uniform}(0.9, 1.1)$, $\text{SV\_Gap} \approx 1.02\times$ ($\text{SNR} = 0.05\times\text{–}0.21\times$). Candidate binary columns generated by thresholding noisy projections contain ~50–65 ones per column, matching the sparsity of the `0-cycle` singletons ($p = 50$ ones at $K=8$) rather than the $2p = 100$ ones of cycle pairs. Consequently, a noisy prediction is ~79 SNPs away from the `0-cycle` and ~100–135 SNPs away from odd cycles, so the minimum distance across all $M_K$ solutions is achieved at the `0-cycle`.

---

## 6. Summary Tables (`results/MEETING_REPORT_TABLES.txt`)

### Table 1: Head-to-Head on the 80 Ambiguous Fukuyama Tasks (`results/table1_80_fukuyama_head_to_head.csv`)

| $K$ ($M_K$) | Pipeline Version | Clean `Exact MSE=0%` | Clean `Mean \|det(A)\|`| Clean Cycle Family Distribution (20 Seeds) | Noisy Step 3 | Noisy Step 4 | Polish Delta | Noisy Dist to GT |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **$K=5$** ($M_5=3$) | `1. PipeA (sol4_v7)` | 95% | `1.40x` | `0-cycle:11, 3-cycle:8, DICT_CHIMERA:1` | `104.35` | `106.87` | `+2.52 SNPs` | `107.89 SNPs` |
| | `2. PipeB (v2_Curr)` | 85% | `1.23x` | `0-cycle:17, 3-cycle:2, DICT_CHIMERA:1` | `97.61` | `102.05` | `+4.44 SNPs` | `103.07 SNPs` |
| | `3. Upgraded v3` | 90% | `1.25x` | `0-cycle:14, 3-cycle:5, DICT_CHIMERA:1` | `101.66` | `101.66` | **`+0.00 SNPs`** | `101.76 SNPs` |
| | `4. Upgraded v4_Std` | **100%** | `1.10x` | `0-cycle:18, 3-cycle:2` | `95.44` | `95.44` | **`+0.00 SNPs`** | `95.44 SNPs` |
| | `5. Upgraded v4_SubReg` | **100%** | `1.10x` | `0-cycle:18, 3-cycle:2` | **`92.59`** | **`92.59`** | **`+0.00 SNPs`** | **`92.59 SNPs`** |
| **$K=7$** ($M_7=71$) | `1. PipeA (sol4_v7)` | 95% | `3.30x` | `3+3-cycle:7, 3-cycle:1, 5-cycle:12` | `107.55` | `114.86` | `+7.31 SNPs` | `115.19 SNPs` |
| | `2. PipeB (v2_Curr)` | 45% | `2.48x` | `0-cycle:18, 3-cycle:2` | `92.00` | `107.65` | `+15.65 SNPs` | `107.65 SNPs` |
| | `3. Upgraded v3` | **100%** | `1.95x` | `0-cyc:6, 3+3-cyc:1, 3-cyc:4, 5-cyc:9` | `92.79` | `92.79` | **`+0.00 SNPs`** | `92.79 SNPs` |
| | `4. Upgraded v4_Std` | **100%** | `1.70x` | `0-cyc:9, 3+3-cyc:1, 3-cyc:5, 5-cyc:5` | `85.51` | `85.51` | **`+0.00 SNPs`** | `85.51 SNPs` |
| | `5. Upgraded v4_SubReg` | **100%** | `1.70x` | `0-cyc:9, 3+3-cyc:1, 3-cyc:5, 5-cyc:5` | **`84.35`** | **`84.35`** | **`+0.00 SNPs`** | **`84.35 SNPs`** |
| **$K=8$** ($M_8=513$) | `1. PipeA (sol4_v7)` | 90% | `4.35x` | `3+3-cycle:4, 5-cycle:4, 7-cycle:12` | `106.52` | `113.03` | `+6.51 SNPs` | `113.28 SNPs` |
| | `2. PipeB (v2_Curr)` | 30% | `2.63x` | `0-cycle:20` | `85.56` | `98.72` | `+13.16 SNPs` | `98.72 SNPs` |
| | `3. Upgraded v3` | **100%** | `2.60x` | `0:1, 3+3:2, 3:8, 5:7, 7:2` | `85.70` | `85.70` | **`+0.00 SNPs`** | `85.70 SNPs` |
| | `4. Upgraded v4_Std` | **100%** | `2.55x` | `0:1, 3+3:1, 3:8, 5:8, 7:2` | `80.30` | `80.30` | **`+0.00 SNPs`** | `80.30 SNPs` |
| | `5. Upgraded v4_SubReg` | **100%** | `2.55x` | `0:1, 3+3:1, 3:8, 5:8, 7:2` | **`79.01`** | **`79.01`** | **`+0.00 SNPs`** | **`79.01 SNPs`** |
| **$K=10$** ($M_{10}=32554$) | `1. PipeA (sol4_v7)` | 30% | `5.65x` | `5+3:2, 7:5, 9:12, DICT_CHIMERA:1` | `99.28` | `105.97` | `+6.68 SNPs` | `105.97 SNPs` |
| | `2. PipeB (v2_Curr)` | 55% | `2.18x` | `0-cycle:20` | `81.86` | `88.84` | `+6.99 SNPs` | `88.84 SNPs` |
| | `3. Upgraded v3` | **100%** | `2.35x` | `0:5, 3+3:1, 3:6, 5+3:2, 5:2, 7:3, 9:1` | `80.92` | `80.92` | **`+0.00 SNPs`** | `80.92 SNPs` |
| | `4. Upgraded v4_Std` | **100%** | `1.95x` | `0-cycle:9, 3-cycle:4, 5-cycle:3, 7-cycle:4` | `74.91` | `74.91` | **`+0.00 SNPs`** | `74.91 SNPs` |
| | `5. Upgraded v4_SubReg` | **100%** | `1.95x` | `0-cycle:9, 3-cycle:4, 5-cycle:3, 7-cycle:4` | **`71.72`** | **`71.72`** | **`+0.00 SNPs`** | **`71.72 SNPs`** |

---

### Table 2: Unaggregated Cell-by-Cell Benchmark Across the 18 Dirichlet Cells (`results/table2_200_benchmark_unaggregated_cells.csv`)

| $K$ | Dirichlet Regime ($\alpha$) | Samples ($n$) | `SV_Gap` ($\sigma_{K-1}/\sigma_K$) | `1. PipeA (sol4_v7)` | `2. PipeB (v2_Curr)` | `3. Upgraded v3` | `4. Upgraded v4_Std` | `5. Upgraded v4_SubReg` |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **$K=5$** | `Balanced (1,1,1,1,1)` | $n=15$ | `4.35x` | `0.00 (100.0%)` | `0.00 (100.0%)` | `0.00 (100.0%)` | **`0.00 (100.0%)`** | **`0.00 (100.0%)`** |
| **$K=5$** | `Balanced (1,1,1,1,1)` | $n=20$ | `4.63x` | `0.00 (100.0%)` | `0.00 (100.0%)` | `0.00 (100.0%)` | **`0.00 (100.0%)`** | **`0.00 (100.0%)`** |
| **$K=5$** | `Imbalanced (5,1,1,1,1)` | $n=15$ | `2.70x` | `0.02 (98.0%)` | `0.00 (100.0%)` | `0.00 (100.0%)` | **`0.00 (100.0%)`** | **`0.00 (100.0%)`** |
| **$K=5$** | `Imbalanced (5,1,1,1,1)` | $n=20$ | `3.39x` | `0.02 (98.0%)` | `0.00 (100.0%)` | `0.00 (100.0%)` | **`0.00 (100.0%)`** | **`0.00 (100.0%)`** |
| **$K=5$** | `Flat (5,5,5,5,5)` | $n=15$ | `2.06x` | `0.04 (96.0%)` | `0.28 (76.0%)` | `0.00 (100.0%)` | **`0.00 (100.0%)`** | **`0.00 (100.0%)`** |
| **$K=5$** | `Flat (5,5,5,5,5)` | $n=20$ | `2.79x` | `0.00 (100.0%)` | `0.12 (88.0%)` | `0.00 (100.0%)` | **`0.00 (100.0%)`** | **`0.00 (100.0%)`** |
| **$K=7$** | `Balanced (1,...,1)` | $n=15$ | `2.48x` | `0.00 (100.0%)` | `0.27 (85.7%)` | `0.00 (100.0%)` | **`0.00 (100.0%)`** | **`0.00 (100.0%)`** |
| **$K=7$** | `Balanced (1,...,1)` | $n=20$ | `3.03x` | `0.00 (100.0%)` | `0.00 (100.0%)` | `0.00 (100.0%)` | **`0.00 (100.0%)`** | **`0.00 (100.0%)`** |
| **$K=7$** | `Imbalanced (5,1,...,1)` | $n=15$ | `1.69x` | `2.34 (81.4%)` | `1.66 (57.1%)` | `0.50 (74.3%)` | **`0.39 (81.4%)`** | `0.47 (80.0%)` |
| **$K=7$** | `Imbalanced (5,1,...,1)` | $n=20$ | `2.35x` | `0.07 (95.7%)` | `0.11 (91.4%)` | `0.01 (98.6%)` | **`0.01 (98.6%)`** | **`0.01 (98.6%)`** |
| **$K=7$** | `Flat (5,...,5)` | $n=15$ | `1.34x` | `4.51 (44.3%)` | `14.76 (5.7%)` | `9.36 (22.9%)` | **`4.29 (22.9%)`** | **`4.29 (22.9%)`** |
| **$K=7$** | `Flat (5,...,5)` | $n=20$ | `1.65x` | **`0.14 (85.7%)`** | `2.46 (25.7%)` | `0.46 (67.1%)` | `0.37 (65.7%)` | `0.37 (65.7%)` |
| **$K=10$** | `Balanced (1,...,1)` | $n=15$ | `1.41x` | `9.87 (58.0%)` | `15.26 (18.0%)` | `4.03 (56.0%)` | `1.57 (65.0%)` | **`1.54 (65.0%)`** |
| **$K=10$** | `Balanced (1,...,1)` | $n=20$ | `1.78x` | `1.31 (85.0%)` | `2.75 (41.0%)` | `0.04 (96.0%)` | **`0.01 (99.0%)`** | **`0.01 (99.0%)`** |
| **$K=10$** | `Imbalanced (5,1,...,1)` | $n=15$ | `1.23x` | `49.67 (19.0%)` | `25.87 (13.0%)` | `8.69 (28.0%)` | **`7.18 (36.0%)`** | `7.39 (35.0%)` |
| **$K=10$** | `Imbalanced (5,1,...,1)` | $n=20$ | `1.54x` | `18.25 (41.0%)` | `9.68 (28.0%)` | `2.04 (70.0%)` | **`0.58 (72.0%)`** | `0.63 (69.0%)` |
| **$K=10$** | `Flat (5,...,5)` | $n=15$ | `1.08x` | `74.00 (0.0%)` | `75.48 (0.0%)` | `72.34 (0.0%)` | `66.87 (0.0%)` | **`62.69 (0.0%)`** |
| **$K=10$** | `Flat (5,...,5)` | $n=20$ | `1.15x` | `43.67 (0.0%)` | `69.62 (0.0%)` | `54.59 (0.0%)` | `47.31 (0.0%)` | **`37.62 (0.0%)`** |

---

## 7. Quick-Start Reproduction Commands

    # 1. Install R package dependencies (if not already installed)
    Rscript install_dependencies.R

    # 2. Print the pre-computed 80-task and 200-task summary tables (< 1 second)
    Rscript evaluation/summarize_80_fukuyama_all_pipelines.R
    Rscript evaluation/summarize_v4_200.R

    # 3. Re-generate all 80 ambiguous Fukuyama tasks from scratch (writes to data/fukuyama_80_tasks/)
    Rscript evaluation/00_generate_80_fukuyama_tasks.R

    # 4. Run all 5 pipelines head-to-head on a single ambiguous task (e.g., Task 41: K=8, Seed 202601)
    Rscript evaluation/benchmark_all_pipelines_on_80_fukuyama.R 41

    # 5. Verify the Simplex Volume Law |det(A)| = 2^c and Unsupervised Support-Depth <= 2 Pair Pool
    Rscript evaluation/verify_volume_and_unsupervised_pair_filter.R

    # 6. Verify the 3 * 2^(d-2) Binary Lattice Law and Higher-Tuple Exact Factorization (W_660)
    Rscript evaluation/verify_exact_vs_beyond_catalog.R
