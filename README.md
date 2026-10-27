# Genomic Haplotype Deconvolution: Audit & Unification of 5 Pipeline Architectures on the Fukuyama Equivalence Manifold

This repository contains the complete, self-contained codebase, pre-generated 80-task ambiguity benchmark datasets ($K \in \{5, 7, 8, 10\}$), **5 deconvolution pipeline implementations**, evaluation scripts, and benchmark results across both the **80 Ambiguous Fukuyama Manifold Tasks** and the **200-Task General Benchmark Suite**.

---

## 1. Repository Structure & Dependencies

All required core R libraries (`phase7/R/`) and pre-computed `.rds`/`.csv` evaluation tables (`results/`) are bundled directly in the repository so every script runs out-of-the-box without path modifications:

    fukuyama-deconvolution-audit/
    ├── README.md                                         # Complete documentation of all 5 pipelines & results
    ├── install_dependencies.R                            # Installs CRAN packages (quadprog, MASS, Matrix, clue, lpSolve)
    ├── master_manifest.csv                               # 200-task general benchmark registry (K=5,7,10 Dirichlet + K=8 Fukuyama)
    ├── phase7/R/                                         # Core plateNMF mathematical library
    │   ├── helper-functions.R                            # Affine projection, Anscombe transform, metric helpers
    │   ├── nmf-functions.R                               # Base subspace & ray-casting functions
    │   ├── nmf-functions_v2.R                            # FCLS QP solver & affine subspace estimators
    │   └── simulation-functions.R                        # Bernoulli endmember & Dirichlet generators
    ├── data/
    │   ├── fukuyama_80_manifest.csv                      # Manifest for the 80 ambiguous tasks (K=5,7,8,10 x 20 seeds)
    │   └── fukuyama_80_tasks/                            # Pre-generated W_true, H_true, S_true, X_counts, P_obs per task
    ├── pipelines/
    │   ├── 01a_pipeA_hunter_sol4_v2.R                    # Pipeline 1 (Stage 1-2): 3-Constraint Exponential SVD Hunter
    │   ├── 01b_pipeA_snapped_subset_v2.R                 # Pipeline 1 (Stage 3): 1/2-Bit Basis Snapping + Warm Exchange
    │   ├── 01c_pipeA_polish_v7.R                         # Pipeline 1 (Stage 4): Binomial Polish v7 (Pass A + Pass B)
    │   ├── 02_pipeB_unified_15plus_v2.R                  # Pipeline 2: Original Unified 15+ (20-SVD + Parsimony Exchange + Circuit Router)
    │   ├── 03_upgraded_v3_pipeline.R                     # Pipeline 3: Unbiased 20-SVD (lambda=0) + Snapped Pool + SV_Gap-Gated Polish v7
    │   └── 04_upgraded_v4_std_and_subreg.R               # Pipelines 4 & 5: MaxVol 20-SVD + Pure Snapped Pool + 4-Start + Subspace-Reg Exchange
    ├── evaluation/
    │   ├── 00_generate_80_fukuyama_tasks.R               # Replicates creation of the 80 ambiguous tasks from scratch
    │   ├── benchmark_all_pipelines_on_80_fukuyama.R      # Head-to-head runner for all 5 pipelines on any task (1-80)
    │   ├── audit_both_pipelines_fukuyama.R               # Multi-run trajectory & gate rejection audit on the 80 tasks
    │   ├── verify_exact_vs_beyond_catalog.R              # Verifies the 3 * 2^(d-2) binary lattice law & 660-mutation triple solution
    │   ├── verify_volume_and_unsupervised_pair_filter.R  # Verifies Simplex Volume Law |det(A)| = 2^c & unsupervised pair filter
    │   ├── summarize_80_fukuyama_all_pipelines.R         # Prints the 80-task Clean & Noisy comparison tables
    │   └── summarize_v4_200.R                            # Prints the 200-task unaggregated cell-by-cell comparison table
    └── results/
        ├── MEETING_REPORT_TABLES.txt                     # Formatted side-by-side console tables for all experiments
        ├── table1_80_fukuyama_head_to_head.csv           # CSV summary of all 5 pipelines on the 80 ambiguous tasks
        ├── table2_200_benchmark_unaggregated_cells.csv   # CSV summary of all 5 pipelines across the 200-task cells
        ├── fukuyama_80_all_pipelines/                    # Saved per-task RDS outputs for the 80 ambiguous tasks
        └── unified_v4_benchmark/                         # Saved per-task RDS outputs for the 200 benchmark tasks

---

## 2. Construction of the 80 Ambiguous Tasks & Theoretical Laws (`00_generate_80_fukuyama_tasks.R`)

For each rank $K \in \{5, 7, 8, 10\}$ (subspace dimension $d = K - 1$, $m = 365$ SNPs, $n = 20$ samples, $20$ seeds `202601..202620`):
1. **Genotype Matrix $W_{\text{true}} = [\mathbf{0}, B] \in \{0,1\}^{365 \times K}$:**
   * Column 1 is the all-zero ancestral strain $\mathbf{0}$.
   * Each of the $d = K - 1$ derived strains has $p = \lfloor 350 / d \rfloor$ private SNPs ($p = 87, 58, 50, 38$ for $K = 5, 7, 8, 10$).
   * Derived Strains 1 and 2 share the remaining $s = 365 - d \cdot p$ SNPs ($s = 17, 17, 15, 23$).
2. **Mixture Proportions $H_{\text{true}} \in \Delta^{K \times 20}$ & Observations:**
   * Unnormalized weights $g_{kj} \sim \operatorname{Uniform}(0.9, 1.1)$ are normalized per column ($H_{kj} = g_{kj} / \sum_{i=1}^K g_{ij}$), keeping every sample strictly inside the interior of the simplex around $(1/K)\mathbf{1}$.
   * Clean expected allele frequencies are $S_{\text{true}} = W_{\text{true}} H_{\text{true}}$, and $100\times$ noisy sequencing read counts are drawn as $X_{mj} \sim \operatorname{Binomial}(100, S_{\text{true}, mj})$ with $P_{\text{obs}} = X / 100$.
3. **Exact Pair-Cycle Equivalence Catalog ($M_K$):**
   Replacing any subset of singletons $\{B e_{v_1}, \dots, B e_{v_\ell}\}$ ($\ell \ge 3$ odd) with the pairwise sums $\{B(e_{v_1}+e_{v_2}), \dots, B(e_{v_\ell}+e_{v_1})\}$ on the non-overlapping graph ($E_{\text{valid}} = \binom{d}{2} - 1$, excluding edge $(1,2)$) yields $W_{\text{alt}} = W_{\text{true}} A$ and $H_{\text{alt}} = A^{-1} H_{\text{true}} > 0$ with $W_{\text{alt}} H_{\text{alt}} = S_{\text{true}}$ exactly:
   * **$K = 5$ ($d=4$):** **$M_5 = 3$** solutions (`1` 0-cycle, `2` 3-cycles)
   * **$K = 7$ ($d=6$):** **$M_7 = 71$** solutions (`1` 0-cycle, `10` 3-cycles, `12` 3+3-cycles, `48` 5-cycles)
   * **$K = 8$ ($d=7$):** **$M_8 = 513$** solutions (`1` 0-cycle, `25` 3-cycles, `95` 3+3-cycles, `192` 5-cycles, `200` 7-cycles)
   * **$K = 10$ ($d=9$):** **$M_{10} = 32,554$** solutions (`1` 0-cycle, `77` 3-cycles, `700` 3+3-cycles, `210` 3+3+3-cycles, `1,302` 5-cycles, `4,704` 5+3-cycles, `10,440` 7-cycles, `15,120` 9-cycles)

### Three Mathematical Laws Verified in Code
1. **The $3 \times 2^{d-2}$ Binary Lattice Law (`verify_exact_vs_beyond_catalog.R`):**
   Evaluating all $2^d$ binary combinations of the $d$ pivot rows on the clean SVD subspace produces exactly **$3 \times 2^{d-2}$** binary vectors in $\{0,1\}^{365}$ ($12, 48, 96, 384$ for $K = 5, 7, 8, 10$), corresponding to all subsets of derived strains that do not simultaneously activate both Strain 1 and Strain 2. Beyond the $1 + d + \binom{d}{2} - 1$ pair-catalog vectors ($10, 21, 28, 45$), the remaining vectors are higher-order tuples (triples $e_a+e_b+e_c$, quads, etc.) that also yield exact $\text{MSE} = 0, H > 0$ factorizations (e.g., the 660-mutation solution `W_660` at $K=8$ with $\max|S_{\text{true}} - W_{660}H_{660}| = 9.99 \times 10^{-16}$ and $\min H_{660} = 0.0365$), proving empirically that $M_K$ is a strict lower bound.
2. **Unsupervised Boolean Support-Inclusion Depth Filter (`verify_volume_and_unsupervised_pair_filter.R`):**
   Relative to the ancestral strain $\mathbf{w}_0$, minimal non-empty support vectors are singletons (depth 1), vectors containing exactly two minimal supports are pairs (depth 2), and vectors containing $\ge 3$ minimal supports are higher-order tuples. Filtering the lattice to support-inclusion depth $\le 2$ recovers the exact pair dictionary ($10, 21, 28, 45$ columns) 100% blindly without ground-truth oracle knowledge.
3. **The Simplex Volume Law $|\det(A)| = 2^c$ (`verify_volume_and_unsupervised_pair_filter.R`):**
   Every odd cycle has an incidence block with determinant $2$, so any solution $W_{\text{alt}} = W_{\text{true}} A$ formed by $c$ vertex-disjoint odd cycles has relative simplex volume $|\det(A)| = \sqrt{\det(\widetilde{W}_{\text{alt}}^T \widetilde{W}_{\text{alt}}) / \det(\widetilde{W}_{\text{true}}^T \widetilde{W}_{\text{true}})} = 2^c$:
   * **0-cycle Ground Truth ($c=0$):** $|\det(A)| = \mathbf{1.0000\times}$ (unique Minimum-Volume Enclosing Simplex)
   * **Single odd cycles ($c=1$: $3, 5, 7, 9$-cycles):** $|\det(A)| = \mathbf{2.0000\times}$
   * **Double odd cycles ($c=2$: $3+3, 5+3$-cycles):** $|\det(A)| = \mathbf{4.0000\times}$
   * **Triple odd cycles ($c=3$: $3+3+3$-cycles):** $|\det(A)| = \mathbf{8.0000\times}$

---

## 3. Detailed Code Walkthrough: How the 5 Pipelines Work & Where They Differ

Every point-estimation pipeline reconstructs $W \in \{0,1\}^{365 \times K}$ across 4 stages:
* **Stage 1 (Subspace Estimation):** Computes sample mean $\mu$ and top $d = K - 1$ left singular vectors $U_b$ of centered allele frequencies $P - \mu \mathbf{1}^T$.
* **Stage 2 (Candidate Harvesting):** Generates a dictionary `W_pool` of candidate binary strain columns.
* **Stage 3A (Basis Snapping) & Stage 3B (Coordinate Exchange):** Optionally snaps candidates onto $(U_b, \mu)$ via 1-bit/2-bit projection descent, then selects the best $K$ columns via cyclic coordinate exchange using Fully Constrained Least Squares (`quadprog::solve.QP`, $H \ge 0, \mathbf{1}^T H = \mathbf{1}$).
* **Stage 4 (Bit-Flipping Binomial Polish):** Tests 1-bit flips (Pass A) and 2-bit pairwise swaps/co-mutations (Pass B) on raw counts $X$, subject to the subspace step-drift ceiling $d_{\text{res}} = 2(1 - 2W_{mk})r_{mk} + M_{mm} \le 0.85$.

### Pipeline 1: `PipeA` (`01a_pipeA_hunter_sol4_v2.R` -> `01b_pipeA_snapped_subset_v2.R` -> `01c_pipeA_polish_v7.R`)
* **Stage 2 (`run_hunter_sol4`):** Draws parametric bootstrap count matrices $X^* \sim \operatorname{Binomial}(100, \hat{P})$ and applies random exponential locus weights $w_i \sim \operatorname{Exp}(1)$ to tilt the SVD hyperplane (`fit_weighted_affine_subspace`). On each tilt, `extract_slawski_vertex_set` selects $d$ well-conditioned pivots (`select_pivots_maxvol`), applies a `2%` Tikhonov ridge pseudo-inverse (`damp_val = 0.02 * norm(U_R, "F")^2`), enumerates all $2^d$ binary anchor combinations, and keeps the top $K$ columns closest to $U$. A proposed $K$-strain set $W^*$ is accepted into the Hunter catalog only if it passes **3 Constraints** against all previously saved sets:
  1. **Overlap Gate:** Fraction of columns matching an existing solution within $\le 1$ SNP must be $\le 60\%$ (`max_overlap_tol = 0.60`).
  2. **Distance Gate:** Mean bidirectional nearest-neighbor Hamming distance must be $\ge \max(15.0, 2.5K)$ SNPs.
  3. **MSE Ratchet:** FCLS reconstruction error must satisfy $\text{MSE}(W^*) \le 1.25 \times \text{best\_mse}$.
* **Stage 3A & 3B (`snap_pool_to_basis` + `select_best_k_subset_multistart`):** Snaps all unique Hunter columns via vectorized 1-bit and 2-bit residual descent onto $U_b$, combines raw and snapped columns (`cbind(W_pool, W_snapped)`), pre-scores all accepted Hunter $K$-sets, and runs 1-column cyclic coordinate exchange from the top 2 Hunter warm starts with **zero parsimony penalty ($\lambda_{\text{pars}} = 0$)**.
* **Stage 4 (`polish_bits_v7`):** Runs up to 3 epochs of **Pass A (abundance-adaptive 1-bit flips)** and **Pass B (deviance-gated 2-bit allele swaps $(1,0) \leftrightarrow (0,1)$ and joint co-mutations $(0,0) \leftrightarrow (1,1)$)** with Pareto checkpointing.

### Pipeline 2: `PipeB` (`02_pipeB_unified_15plus_v2.R`)
* **Geometric Manifold Router & Stage 4 Circuit Walker:** First tests whether the clean SVD subspace admits valid odd-cycle circuits (`circuits_found > 0`). When triggered on clean manifold data, it routes to the **Circuit Walker**, which enumerates all certified odd-cycle equivalence classes ($M_K = 3, 71, 513, 32554$).
* **Stage 2 (`harvest_dual_20svd`, Point-Estimation Route):** Removes the 3-constraint rejection gate and harvests ~600-800 unique candidates across **20 SVDs**: `1` baseline SVD + `14` locus-weighted SVDs ($w_{\text{loc}} \sim \operatorname{Exp}(1)$) + `5` sample-weighted SVDs ($w_{\text{samp}} \sim 0.8\operatorname{Dir}(1) + 0.2/n$, which amplify rare/imbalanced strains), combining `extract_slawski` (`2%` Tikhonov damping), `extract_euclidean_rays`, and `extract_arcsine_rays`.
* **Stage 3 (`coord_exchange_robust_parsimony`):** Skips Stage 3A Basis Snapping and runs 3-start cyclic coordinate exchange (`Restart 1: QR MaxVol`, `Restart 2: Leverage Swap`, `Restart 3: Sparse MinVol`) with explicit L1 parsimony penalty **$\lambda_{\text{pars}} = 5 \times 10^{-5}$**.
* **Stage 4 (`polish_strains`):** Runs Pass A (1-bit flips only; omits Pass B 2-bit pairwise swaps) unconditionally.

### Pipeline 3: `Upgraded v3` (`03_upgraded_v3_pipeline.R`)
* **What Changed from `PipeB`:**
  1. **Zero Parsimony Penalty ($\lambda_{\text{pars}} = 0.0$):** Removes the L1 mutation penalty in Stage 3B so all valid cycle solutions have identical loss ($\text{MSE} = 0$).
  2. **Adds Stage 3A Basis Snapping (`snap_pool_to_basis`):** Snaps the 20-SVD candidate pool onto $(U_b, \mu_b)$ via 1-bit/2-bit descent (`cbind(W_pool, W_snapped)`).
  3. **Upgrades Stage 4 to `polish_bits_v7_full` (Pass A + Pass B) with the `SV_Gap > 1.12` Gate:** Computes the empirical singular value gap $\text{SV\_Gap} = \sigma_{K-1}(P_{\text{cent}}) / \sigma_K(P_{\text{cent}})$. When $\text{SV\_Gap} > 1.12$, full Pass A + Pass B Polish runs; when $\text{SV\_Gap} \le 1.12$ (sub-BBP noisy $100\times$ Fukuyama data where $\text{SV\_Gap} \approx 1.02\times$), Stage 4 Polish is automatically bypassed so it cannot overfit shot noise.

### Pipeline 4: `Upgraded v4_Std` (`04_upgraded_v4_std_and_subreg.R`, `use_subspace_reg = FALSE`)
* **What Changed from `v3`:**
  1. **MaxVol Pivot Conditioning + Adaptive Tikhonov Damping (`extract_slawski_maxvol`):** Integrates `PipeA`'s `select_pivots_maxvol` into all 20 SVDs and sets Tikhonov damping to **`0.0` on clean subspaces** ($\sigma_{K-1}/\sigma_K > 10^4$, preventing shared-SNP shrinkage) and `0.02` on noisy subspaces.
  2. **Pure Snapped Pool (`snap_pool_only`):** Replaces raw unsnapped random rays with their snapped facets on $(U_b, \mu_b)$ instead of keeping unsnapped noise columns in the pool.
  3. **4-Start Hybrid Coordinate Exchange (`coord_exchange_v4`):** Combines **both** `PipeA`'s top 2 Slawski $K$-simplex warm starts and `PipeB`'s `QR MaxVol` and `Sparse` starts in a single tournament with $\lambda_{\text{pars}} = 0.0$.

### Pipeline 5: `Upgraded v4_SubReg` (`04_upgraded_v4_std_and_subreg.R`, `use_subspace_reg = TRUE`)
* **What Changed from `v4_Std`:**
  * During Stage 3B candidate swaps, evaluates FCLS against the **rank-$d$ SVD-denoised frequency matrix** $\hat{P}_{\text{svd}} = \Pi_{[0,1]}(\mu_b \mathbf{1}^T + U_b U_b^T(P_{\text{obs}} - \mu_b \mathbf{1}^T))$ plus an off-plane subspace residual penalty $0.03 \frac{\|(I - U_b U_b^T)(W_{\text{sub}} - \mu_b)\|_F^2}{m K}$.
  * On noisy Flat Dirichlet mixtures, scoring against $\hat{P}_{\text{svd}}$ prevents candidate swaps from overfitting orthogonal shot noise. On clean Fukuyama manifold tasks, $\hat{P}_{\text{svd}} = S_{\text{true}}$ and all manifold vertices lie exactly on $U_b$ (`cand_res_sq = 0`), preserving **zero bias ($\lambda_{\text{pars}} = 0$) and 100% exact manifold recovery**.

---

## 4. Complete Benchmark Results (`results/MEETING_REPORT_TABLES.txt`)

### Table 1: Head-to-Head on the 80 Ambiguous Fukuyama Tasks (`results/table1_80_fukuyama_head_to_head.csv`)

| $K$ ($M_K$) | Pipeline Version | Clean `Exact MSE=0%` | Clean `Mean \|det(A)\|`| Clean Cycle Family Distribution (20 Seeds) | Noisy Step 3 | Noisy Step 4 | Polish Delta | Noisy Dist to GT |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **$K=5$** ($M_5=3$) | `1. PipeA (sol4_v7)` | 95% | `1.40x` | `0-cycle:11, 3-cycle:8, CHIMERA:1` | `104.35` | `106.87` | `+2.52 SNPs` | `107.89 SNPs` |
| | `2. PipeB (v2_Curr)` | 85% | `1.23x` | `0-cycle:17, 3-cycle:2, CHIMERA:1` | `97.61` | `102.05` | `+4.44 SNPs` | `103.07 SNPs` |
| | `3. Upgraded v3` | 90% | `1.25x` | `0-cycle:14, 3-cycle:5, CHIMERA:1` | `101.66` | `101.66` | **`+0.00 SNPs`** | `101.76 SNPs` |
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
| **$K=10$** ($M_{10}=32554$) | `1. PipeA (sol4_v7)` | 30% | `5.65x` | `5+3:2, 7:5, 9:12, CHIMERA:1` | `99.28` | `105.97` | `+6.68 SNPs` | `105.97 SNPs` |
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

## 5. Quick-Start Reproduction Commands

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
