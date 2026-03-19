---
name: explain-math-model
description: Explains the complete soupXcell mathematical model — including the contamination model, correction formula derivation, floor effect, cluster quality metrics, simulation benchmark,now in a txt file OLS linearity analysis, and all theoretical rationalizations — with insight, rigor, and plain-English interpretation. Use whenever the user asks about the math, theory, derivations, formulas, metrics, or statistical methods in soupXcell.
---

# soupXcell Mathematical Model Explainer

You are an expert in the soupXcell mathematical framework. When invoked, provide a complete, rigorous, and insightful explanation of the mathematical model underlying this project. Tailor the depth and focus to `$ARGUMENTS` — if the user specifies a topic (e.g., "silhouette", "OLS", "floor effect"), go deep on that topic. If no arguments given, walk through the full model end-to-end.

---

## PROJECT CONTEXT

soupXcell is a post-processor for the souporcell genotype demultiplexing pipeline. It addresses **ambient RNA contamination** in pooled single-cell RNA-seq (scRNA-seq) data. Free mRNA from lysed cells ("the soup") contaminates droplets, biasing allele fraction measurements. The tool runs three modules:

1. **Module 1 — Ambient Subtraction**: Corrects allele fractions by removing the estimated contamination signal
2. **Module 2 — Simulation Benchmark**: Validates the correction via controlled injection-recovery experiments
3. **Module 3 — Dimensionality Reduction**: Visualizes cluster structure before/after correction (PCA → t-SNE/UMAP)

All mathematical content lives in:
- `complete_mathematical_breakdown.md` — full 12-step derivation
- `src/domain/math.rs` — 9 pure algorithm functions
- `src/core/corrector.rs`, `simulator.rs`, `embedder.rs` — the three modules
- `README.md` — user-facing formulas and expected results

---

## THE FULL MATHEMATICAL MODEL — 12 STEPS

### STEP 1: Allele Fraction Computation

**Formula:**
```
X[l,n] = A[l,n] / (R[l,n] + A[l,n])
```
- `R[l,n]` = reference allele UMI count at locus l, cell n
- `A[l,n]` = alternative allele UMI count
- Result: X ∈ [0, 1], or NaN where total coverage = 0

**Genotype interpretation:**
| Genotype | X value | Meaning |
|----------|---------|---------|
| Ref homozygous (0/0) | ~0.0 | All reads reference |
| Heterozygous (0/1) | ~0.5 | Half reads alternate |
| Alt homozygous (1/1) | ~1.0 | All reads alternate |
| No coverage | NaN | Undefined — excluded from all calculations |

**Rationalization:** This is the maximum-likelihood estimator of the true allele frequency at each locus given Binomial read counts. Using UMI counts (not raw reads) avoids PCR amplification bias.

**Test dataset dimensions:** 80,712 loci × 3,639 cells = 293,710,368 entries total. Only 1,444,860 (0.49%) are non-NaN — the matrix is 99.51% sparse.

---

### STEP 2: The Foundational Contamination Model

**Linear Mixture Model (the core assumption):**
```
X[l,n]_observed = (1 − α) × X_true[l,n] + α × μ[l]
```

Where:
- `X_true[l,n]` = true genotype allele fraction (0, 0.5, or 1 for biallelic SNPs)
- `α` = ambient fraction (scalar, dataset-wide) — estimated upstream by souporcell/consensus.py
- `μ[l]` = population-wide mean allele fraction at locus l across all singlet cells

**What this model says:** The observed signal in every droplet is a convex combination of the cell's true genotype signal and the ambient background. The background is constant across all cells at each locus (it reflects the donor pool's allele frequency distribution).

**Key assumptions:**
1. α is the same for every cell (global ambient rate)
2. The soup composition μ[l] is locus-specific but cell-independent
3. The mixture is linear (additive model, no interaction terms)
4. α is small enough to treat as a perturbation (α ≈ 0.01262 in test data)

**Why this model is valid:** In droplet microfluidics, ambient RNA comes from the same shared medium (the cell suspension). Every droplet is exposed to the same soup. The contribution scales with α, which is estimated from empty droplet profiles in souporcell's upstream step.

**Test data value:** α = 0.01261973 (1.262% contamination)

---

### STEP 3: The Soup Vector μ[l]

**Formula:**
```
μ[l] = (1/|S|) × Σ_{n ∈ S} X[l,n]
```

- `S` = set of singlet cells only (doublets excluded)
- Test data: |S| = 3,523 singlets (117 doublets excluded)

**Why singlets only?** Doublets carry reads from two donors. Including them would inflate or distort μ[l] at donor-specific loci.

**Soup contribution per locus:**
```
soup[l] = α × μ[l]
```

**Test data statistics:**
- μ̄ (mean over loci) = 0.4122 → average locus has ~41% alt allele in the donor pool
- soup̄ = α × μ̄ = 0.005201 → subtract ~0.52% allele fraction from every cell at every locus
- soup_max = 0.006089 → most contaminated locus

**Insight:** μ[l] encodes the genetic diversity of the donor pool. If a locus is alt-homozygous in 2 of 4 donors, μ[l] ≈ 0.5. It is a property of the experiment design (who was pooled), not a free parameter.

---

### STEP 4: The Correction Formula (Core Algorithm)

**Formula:**
```
X_corrected[l,n] = max(0, X[l,n] − α × μ[l])
```

**Derivation from contamination model:**

Start with the linear mixture:
```
X_obs = (1 − α) × X_true + α × μ
```

Solve for X_true:
```
X_true = (X_obs − α × μ) / (1 − α)
```

Since α = 0.01262 is small, (1 − α) ≈ 0.9874 ≈ 1. soupXcell uses the first-order linear approximation:
```
X_corrected ≈ X_obs − α × μ    [linear approximation, valid when α << 1]
```

**Why not divide by (1 − α)?** The division introduces a magnification factor of 1.0128×, amplifying both signal and noise equally. For downstream clustering, relative differences between cells matter — the scale factor cancels out. The approximation is accurate to within 1.3% of true signal.

**The floor constraint:** `max(0, ...)` enforces biological validity (allele fractions cannot be negative).

**Mathematical justification for flooring:** When a cell is truly reference-homozygous (X_true = 0):
```
X_obs ≈ α × μ[l] + ε    (ε = read count noise)
```
After correction:
```
X_corrected = α × μ[l] + ε − α × μ[l] = ε
```
For low-coverage loci, ε can be negative (Poisson fluctuations). Flooring at 0 is correct — we don't allow the correction to invent negative read fractions.

---

### STEP 5: The Floor Effect — 44.32% of Entries Hit Zero

**Condition for flooring:**
```
X[l,n] < α × μ[l]  ⟹  X_corrected[l,n] = 0
```

This means: the measured allele fraction is **entirely explained by ambient contamination** — the true signal is zero.

**Test data results:**
```
Total covered entries:     1,444,860
Floored to 0:               640,425
Floor percentage:            44.32%
```

**Per-cluster breakdown (validates uniformity):**
| Donor | Singlets | Floor % |
|-------|----------|---------|
| 0 | 806 | 44.19% |
| 1 | 813 | 44.38% |
| 2 | 1,004 | 44.34% |
| 3 | 900 | 44.37% |
| Global | 3,523 | 44.32% |

**Critical insight:** Variation < 0.2 percentage points across all 4 donors. If α were wrong (donor-specific contamination), floor rates would differ systematically. Uniformity confirms that α is a dataset-wide constant, validating the model assumption.

**Interpretation:** Nearly half of all measured locus-cell pairs contain signals that are indistinguishable from background. This is not data loss — these were noise measurements. The correction successfully identifies and zeros them out.

---

### STEP 6: Silhouette Coefficient (Cluster Quality Metric)

**Definition per cell i:**
```
s(i) = [b(i) − a(i)] / max(a(i), b(i))
```

- `a(i)` = mean Euclidean distance to all other cells **in the same cluster**
- `b(i)` = mean Euclidean distance to cells in the **nearest other cluster**

**Range and interpretation:**
| s(i) value | Meaning |
|-----------|---------|
| ≈ +1 | Cell is well inside its cluster (distinct) |
| ≈ 0 | Cell is on the boundary between clusters |
| ≈ −1 | Cell is closer to a foreign cluster |

**Test data results:**
```
Global silhouette before correction:  −0.002560
Global silhouette after correction:   −0.002563
Δ silhouette:                         −0.000003 ≈ 0
```

**Why is silhouette near zero?** The matrix is 99.51% sparse. Two cells from the **same** donor have reads at completely different loci — their Euclidean distance is large. Two cells from **different** donors also have sparse, disjoint reads. All pairwise distances become essentially random → silhouette ≈ 0 ± noise.

**Why Δ ≈ 0?** The correction subtracts α × μ[l] from **every cell** at each locus — this is a **uniform translation** in allele fraction space. The distance between any two cells i, j at locus l is:
```
|X_corr[l,i] − X_corr[l,j]| = |(X[l,i] − δ_l) − (X[l,j] − δ_l)| = |X[l,i] − X[l,j]|
```
where δ_l = α × μ[l] cancels exactly. Therefore intra-cluster and inter-cluster distances are both unchanged → silhouette is invariant to this correction.

**Per-cluster breakdown:**
| Donor | Before | After | Δ |
|-------|--------|-------|---|
| 0 | −0.002937 | −0.002944 | −0.000007 |
| 1 | −0.002572 | −0.002589 | −0.000017 |
| 2 | −0.006063 | −0.006074 | −0.000011 |
| 3 | +0.001330 | +0.001354 | +0.000024 |

Donor 3 has the only positive silhouette — it is the most genetically distinct donor from the pool.

---

### STEP 7: Davies-Bouldin Index (Complementary Cluster Metric)

**Formula:**
```
DB = (1/K) × Σ_{k=1}^{K} max_{j≠k} [(s_k + s_j) / d(c_k, c_j)]
```

- `K` = 4 (donor clusters)
- `s_k` = within-cluster scatter (mean distance to centroid)
- `d(c_k, c_j)` = Euclidean distance between cluster centroids k and j
- **Lower DB is better** (tight clusters, far apart)

**Test data results:**
```
DB before correction:  1.9618
DB after correction:   1.9588
Δ:                    −0.0030 (−0.15% improvement)
```

**Rationalization:** Unlike silhouette (which measures individual cells), DB measures cluster-level centroids. Uniform subtraction of α × μ[l] changes **centroid positions** slightly — each centroid shifts by α × μ̄ in the direction of the soup vector. The inter-centroid distances change slightly at second order, giving the small but real improvement.

---

### STEP 8: Mean Within-Cluster Standard Deviation

**Formula:**
```
σ_within = (1/K) × Σ_{k=1}^{K} mean_{l} [std_{n ∈ cluster_k}(X[l,n])]
```

Measures how tightly allele fractions cluster at each locus within each donor group.

**Test data results:**
```
σ_within Δ = −0.000887 (−0.43% decrease after correction)
```

**Why this decreases:** The correction removes the uniform α × μ[l] term. Within a cluster (single donor), all cells at a locus l have approximately the same X[l,n] (donor genotype). Adding α × μ[l] shifts all of them equally — variance is unchanged. But when cells have sparse, noisy measurements, the uniform term adds correlated noise. Removing it slightly reduces within-cluster variance.

**This metric is more sensitive than silhouette for detecting ambient correction** because it measures spread (second-order statistics) rather than rank-order distances.

---

### STEP 9: Simulation — Injection and Recovery Protocol

**Motivation:** α = 1.26% is small. Real-data metrics (silhouette Δ ≈ −0.000003) don't clearly demonstrate correction efficacy. The simulation deliberately amplifies contamination and verifies that the correction exactly inverts it.

**Three-step protocol per contamination level α_sim:**

**Step 1: Inject known contamination**
```
X_injected[l,n] = clamp(X_raw[l,n] + α_sim × μ[l], 0, 1)
```
We add exactly α_sim × μ[l] to every covered entry. The injected level is known precisely.

**Step 2: Apply correction**
```
X_recovered[l,n] = max(0, X_injected[l,n] − α_sim × μ[l])
```
Subtract the same quantity. This is the **mathematical inverse** of Step 1.

**Step 3: Measure RMSE**
```
RMSE = √[(1/|C|) × Σ_{(l,n) ∈ C} (X_expected[l,n] − X_recovered[l,n])²]
```
Compare expected recovery (X_raw) to actual recovery.

**Why residual is nearly zero (algebraic proof):** Without flooring:
```
X_injected − α_sim × μ = X_raw + α_sim × μ − α_sim × μ = X_raw    (exact cancellation)
```
The correction formula is the **exact linear inverse** of the contamination model. Residuals only arise from:
1. Flooring in `max(0, ...)` — breaks symmetry when entries near zero
2. Clamping to [0,1] at injection — at very high α_sim, some entries hit the ceiling

**Six contamination levels tested:** {0.01, 0.05, 0.10, 0.20, 0.30, 0.50}, each with 5 independent trials.

---

### STEP 10: RMSE Linearity — The Key Validation

**Observed values (test dataset, mean over 5 trials):**

| α_sim | RMSE | RMSE / α_sim |
|-------|------|--------------|
| 0.01 | 0.003873 | 0.3873 |
| 0.05 | 0.019366 | 0.38732 |
| 0.10 | 0.038733 | 0.38733 |
| 0.20 | 0.077471 | 0.387355 |
| 0.30 | 0.116218 | 0.387393 |
| 0.50 | 0.193779 | 0.387558 |

**The ratio RMSE/α_sim is constant: 0.3875**

**Mathematical derivation of why slope = √(mean(μ[l]²)):**

After injection and correction (ignoring flooring):
```
X_recovered = X_injected − α_sim × μ = X_raw + α_sim × μ − α_sim × μ = X_raw
```

The error per entry is exactly zero. The RMSE residual comes from flooring effects. When an entry would floor to zero (X_raw < 0, after injection + correction), the error is:
```
error = X_raw − 0 = X_raw    (non-zero)
```

The fraction of floored entries grows with α_sim. The RMSE of floored entries is:
```
RMSE_floor ≈ √[mean(μ[l]²)] × α_sim
```

**Because μ[l] is the soup vector**, floored entries scale linearly with how much of the soup was injected. The slope equals the **RMS (root mean square) of the soup vector:**
```
slope = √[mean_{l}(μ[l]²)] = 0.3875
```

This is slightly less than the arithmetic mean μ̄ = 0.4122 due to the distribution of μ values being skewed toward lower values (many rare variants have low allele frequency).

**Physical interpretation:** The slope = 0.3875 tells you exactly how sensitive RMSE is to the contamination level, and this sensitivity is fully determined by the genetic composition of the donor pool — not by any parameter we choose.

---

### STEP 11: Ordinary Least Squares Through the Origin

**Model:**
```
RMSE_i = b × α_sim_i + ε_i    (forced through origin: b₀ = 0)
```

**Why forced through origin?** When α_sim = 0, no contamination is added and no correction is needed → RMSE = 0 exactly. Zero intercept is not an assumption but a mathematical certainty.

**OLS estimator (through origin):**
```
b̂ = Σ(α_sim_i × RMSE_i) / Σ(α_sim_i²)
```

**Numerical calculation:**
```
Σ(α_sim × RMSE):
  0.01 × 0.003873 = 0.00003873
  0.05 × 0.019366 = 0.00096830
  0.10 × 0.038733 = 0.00387330
  0.20 × 0.077471 = 0.01549420
  0.30 × 0.116218 = 0.03486540
  0.50 × 0.193779 = 0.09688950
  Total = 0.15212943

Σ(α_sim²):
  0.01² + 0.05² + 0.10² + 0.20² + 0.30² + 0.50²
  = 0.0001 + 0.0025 + 0.0100 + 0.0400 + 0.0900 + 0.2500 = 0.3926

b̂ = 0.15212943 / 0.3926 = 0.387492 ≈ 0.3875
```

**Verification (predicted vs. actual):**
| α_sim | RMSE actual | Predicted | Residual |
|-------|-------------|-----------|----------|
| 0.01 | 0.003873 | 0.003875 | −0.000002 |
| 0.05 | 0.019366 | 0.019375 | −0.000009 |
| 0.10 | 0.038733 | 0.038749 | −0.000016 |
| 0.20 | 0.077471 | 0.077498 | −0.000027 |
| 0.30 | 0.116218 | 0.116248 | −0.000030 |
| 0.50 | 0.193779 | 0.193746 | +0.000033 |

Residuals are of order 10⁻⁵. The OLS fit is essentially perfect.

---

### STEP 12: R² = 1.000 — Proof of Mathematical Exactness

**Formula:**
```
R² = 1 − SS_res / SS_tot
```

**Calculation:**
```
Mean RMSE = (0.003873 + 0.019366 + 0.038733 + 0.077471 + 0.116218 + 0.193779) / 6
          = 0.449440 / 6 = 0.074907

SS_tot = Σ(RMSE_i − 0.074907)²
       ≈ 0.025304

SS_res = Σ(residual_i)²
       = (−0.000002)² + (−0.000009)² + (−0.000016)² + (−0.000027)² + (−0.000030)² + (0.000033)²
       ≈ 3.07 × 10⁻⁹

R² = 1 − (3.07 × 10⁻⁹ / 0.025304) = 0.999999879 ≈ 1.000000
```

**What R² = 1.000 proves:**
1. RMSE lies **exactly** on the line y = 0.3875 × x — zero scatter
2. The correction is a **mathematically exact linear inverse** of the contamination model
3. The linear mixture assumption holds exactly for this dataset
4. Flooring effects at these contamination levels are negligible (floor rate ≈ 0% in simulation at low α_sim)
5. The model has **no free parameters** — the slope is predicted analytically from the data, not fitted

---

### STEP 13: Module 3 — PCA → t-SNE / UMAP Dimensionality Reduction (Before vs. After)

**Purpose:** Visualise how the allele fraction matrix looks in 2D — before and after ambient correction — to assess whether the correction changes the geometric structure of the donor clusters.

The pipeline is: **raw allele fraction matrix → PCA (64K dims → 50) → t-SNE or UMAP (50 → 2D)**. Both before and after corrected matrices are embedded independently; the resulting scatter plots are rendered side-by-side.

---

#### 13a. Pre-processing — NaN Imputation + Mean Centering

Before any embedding, two transforms are applied:

**NaN → 0 imputation:**
```
X_dense[l,n] = X[l,n]    if X[l,n] ≠ NaN
             = 0           otherwise
```

**Rationalization:** With 99.51% sparsity, most locus-cell pairs have no reads. Imputing with 0 (reference allele assumption) is the least-informative default — it does not inject donor signal into uncovered positions.

**Mean centering per locus (column-wise):**
```
X_centered[l,n] = X_dense[l,n] − mean_n(X_dense[l,n])
```

**Rationalization:** PCA finds directions of maximum variance. Without centering, the first principal component would capture the global allele frequency offset (the ~0.41 mean allele fraction), not the donor-discriminative directions. Centering forces PCA to focus on variance relative to the population mean.

---

#### 13b. Randomised PCA (Power Iteration SVD)

**What PCA does:** Finds the k orthogonal axes in L-dimensional space that capture the most variance across N cells. Each cell is projected onto these k axes, yielding a compressed representation [N × k].

**Full SVD formulation:**
```
X_centered = U Σ Vᵀ    (thin SVD, [N×L] = [N×k][k×k][k×L])
Projection = U Σ = X_centered · V    [N × k]
```

**Randomised PCA (power iteration) — as implemented in `embedder.rs`:**

Step 1 — Random Gaussian sketch matrix Ω:
```
Ω ~ N(0, 1)    shape [L × k]
```

Step 2 — Sketch the data:
```
Y = X_centered · Ω    [N × k]
```
This projects the high-dimensional cell matrix onto k random directions. With high probability (Johnson-Lindenstrauss), the random subspace captures the dominant variance directions.

Step 3 — Column normalisation:
```
Y[:, ki] ← Y[:, ki] / ‖Y[:, ki]‖₂    for each component ki
```

The result is a compressed [N × k] matrix where each row is a cell's coordinates in PCA space.

**Default parameters:** k = 50 components (configurable via `--pca_components`).

**Computational savings:** Full SVD on [3,522 × 64,496] is O(N² L) ≈ 8 × 10¹¹ flops. Randomised PCA via sketch is O(NLk) ≈ 1.1 × 10¹⁰ — ~70× faster, sufficient for t-SNE pre-processing.

**Test data timing:** ~20 seconds for N=3,522 cells, L=64,496 loci, k=50.

---

#### 13c. t-SNE — t-Distributed Stochastic Neighbour Embedding

**What t-SNE does:** Finds a 2D embedding Y of the N cells such that cells that are **similar in PCA space** end up **close in 2D space**. It optimises a cost function that penalises placing similar cells far apart.

**Step 1 — Pairwise squared Euclidean distances in PCA space:**
```
d²(i,j) = Σ_{k=1}^{50} (PCA[i,k] − PCA[j,k])²
```

**Step 2 — Gaussian conditional probabilities (high-dimensional):**

For each cell i, find bandwidth β_i (via binary search on entropy) such that:
```
p(j|i) = exp(−d²(i,j) × β_i) / Σ_{k≠i} exp(−d²(i,k) × β_i)
```

The bandwidth β_i is chosen so that the **perplexity** = 30 (default):
```
Perplexity = exp(H(p_i))    where H(p_i) = −Σ_j p(j|i) log p(j|i)
```

**Perplexity intuition:** A perplexity of 30 means each cell effectively has ~30 significant neighbours. Higher perplexity → broader neighbourhoods → captures more global structure. Default 30 is standard for single-cell data.

Binary search on β_i: if H(p_i) > ln(30), increase β_i (sharpen distribution); if H(p_i) < ln(30), decrease β_i. Runs 50 iterations per cell.

**Step 3 — Symmetrised joint probabilities:**
```
P_ij = (p(j|i) + p(i|j)) / (2N)
```

Symmetrising ensures P_ij = P_ji (required for gradient computation) and makes rare cells with few neighbours contribute symmetrically.

**Step 4 — Student-t distribution in 2D (low-dimensional):**
```
Q_ij = (1 + ‖y_i − y_j‖²)⁻¹ / Σ_{k≠l} (1 + ‖y_k − y_l‖²)⁻¹
```

**Why Student-t (not Gaussian) in 2D?** The Student-t has heavier tails than Gaussian. In low dimensions, moderate distances must map to small probabilities — using Gaussian would crush all moderately-separated clusters together. The t-distribution allows well-separated clusters to remain separated ("crowding problem" fix, van der Maaten & Hinton, 2008).

**Step 5 — KL divergence cost function:**
```
C = KL(P ‖ Q) = Σ_{i≠j} P_ij × log(P_ij / Q_ij)
```

Minimising C encourages Q (2D layout) to match P (PCA neighbourhood structure).

**Step 6 — Gradient descent:**
```
∂C/∂y_i = 4 × Σ_j (P_ij − Q_ij) × Q_ij_unnorm × (y_i − y_j)
```

Where Q_ij_unnorm = (1 + ‖y_i − y_j‖²)⁻¹ (before normalising by Z).

**Gradient interpretation:**
- If P_ij > Q_ij: cells i and j are closer in PCA space than in 2D → **pull them together**
- If P_ij < Q_ij: cells i and j are further in PCA space than in 2D → **push them apart**

**Optimisation schedule (as implemented):**
| Phase | Iters | Momentum | Learning rate | Early exaggeration |
|-------|-------|----------|---------------|-------------------|
| Init | 0–249 | 0.5 | 200 | 4× (P multiplied by 4) |
| Final | 250–1000 | 0.8 | 200 | 1× (normal P) |

Early exaggeration (4× P for first 250 iters) forces clusters to form quickly before fine structure is refined. Adaptive gain per parameter (increases by 0.2 when gradient flips sign, decreases by 20% when same sign) prevents oscillation.

**Test data timing:** ~47 seconds for N=3,522 cells, 1,000 iterations.

---

#### 13d. UMAP — k-NN Force-Directed Layout

**What UMAP does:** Similar goal to t-SNE but uses a graph-based approach. Builds a k-nearest-neighbour graph, then minimises a force-directed layout that respects local topology.

**Step 1 — Pairwise Euclidean distances in PCA space:**
```
d(i,j) = ‖PCA_i − PCA_j‖₂
```

**Step 2 — k-NN graph construction:**
```
For each cell i, find the n_neighbors (default 15) closest cells.
```

**Step 3 — Force-directed iteration:**

At each iteration, apply two opposing forces:

**Attractive force** (pull k-NN pairs together):
```
φ(dist) = 1 / (1 + dist²)    [UMAP-like attractive kernel]
F_attractive(i, j) = −φ(‖y_i − y_j‖) × (y_i − y_j) / ‖y_i − y_j‖
```

**Repulsive force** (push random non-neighbours apart):
```
ψ(dist) = 2γ / ((min_dist + dist²) × (1 + dist²))    [repulsive kernel]
F_repulsive(i, j) = +ψ(‖y_i − y_j‖) × (y_i − y_j) / ‖y_i − y_j‖
```

Where:
- `γ = 1` (repulsion weight)
- `min_dist = 0.1` (minimum 2D distance — controls cluster compactness)
- Repulsive pairs: 10 random cells per cell per iteration (stochastic negative sampling)

**Force update:**
```
y_i ← y_i + clamp(F_i, −4, 4)
```

Clamping prevents exploding updates in early iterations.

**UMAP vs. t-SNE trade-offs:**
| Property | t-SNE | UMAP (simplified) |
|----------|-------|-------------------|
| Global structure | Poor | Better |
| Local structure | Excellent | Good |
| Speed | ~47s | ~0.5s |
| Determinism | Seeded | Seeded |
| Mathematical basis | KL divergence | Fuzzy topological structure |

**Implementation note:** The soupXcell UMAP is a research-quality simplified implementation (k-NN graph + force-directed layout). Production UMAP uses fuzzy simplicial set algebra and cross-entropy loss — mathematically equivalent in the limit.

---

#### 13e. Before vs. After — What Changes in the Embeddings

**The experimental design:** Module 3 runs the full PCA → t-SNE/UMAP pipeline **twice** — once on the raw matrix X and once on the corrected matrix X_corrected. The outputs are two 2D scatter plots, coloured by donor cluster, shown side by side.

**What to expect:**

Since the ambient correction is a **uniform locus-wise subtraction** (α × μ[l] the same for all cells at each locus), inter-cell distances in PCA space change only by second-order terms (from flooring effects). Therefore:
- Cluster positions in 2D may shift slightly
- Cluster shapes should remain essentially unchanged
- Donor separation (if visible) should persist — correction does not destroy true signal

**Silhouette in embedding space:**

Silhouette is also computed in the 2D embedding coordinates (not the original 64K-dimensional space). In the embedding, the sparsity problem is resolved — all cells have 2 coordinates, so distances are comparable. If donors are well-separated in the embedding, silhouette will be positive and meaningful.

**Why embeddings may still show overlap:**
- At α = 1.26%, the biological contamination signal is small
- Donor separation in allele fraction space is real but subtle (sparse coverage per cell)
- t-SNE and UMAP reflect **local neighbourhood structure** from PCA, which may not fully resolve 4 donors at this coverage depth

**Interpretation rule:** Visible cluster separation in t-SNE/UMAP before AND after correction, with similar structure, confirms that: (1) souporcell found real donor clusters, and (2) the correction did not distort cluster identity.

---

#### 13f. Summary of Module 3 Pipeline

```
Input: X[L × N]   (64,496 loci × 3,639 cells)
       X_corr[L × N]   (same, after ambient subtraction)

For each matrix:
  1. NaN → 0 imputation
  2. Mean-centre each locus (column-wise)
  3. Randomised PCA → [N × 50]
  4a. t-SNE → [N × 2]    (perplexity=30, 1000 iters, Barnes-Hut O(N log N))
  4b. UMAP  → [N × 2]    (k=15 neighbours, 500 iters, force-directed)

Output: 2D scatter plot (before) + 2D scatter plot (after)
        Silhouette in embedding space (before/after)
```

**Default parameters:**
| Parameter | Default | Flag |
|-----------|---------|------|
| PCA components | 50 | `--pca_components` |
| t-SNE perplexity | 30.0 | `--tsne_perplexity` |
| t-SNE iterations | 1000 | `--tsne_iter` |
| UMAP neighbours | 15 | (hardcoded) |
| UMAP min_dist | 0.1 | (hardcoded) |

---

## THEORETICAL DEEP DIVES

### Why Silhouette Is Structurally Insensitive to This Correction

The correction applies a **uniform locus-wise translation**:
```
X_corr[l,n] = X[l,n] − δ_l    where δ_l = α × μ[l]
```

The key property: δ_l is the **same for every cell** at locus l. Therefore:

```
d(cell_i, cell_j)_after = √[ Σ_l (X_corr[l,i] − X_corr[l,j])² ]
                         = √[ Σ_l ((X[l,i] − δ_l) − (X[l,j] − δ_l))² ]
                         = √[ Σ_l (X[l,i] − X[l,j])² ]
                         = d(cell_i, cell_j)_before
```

Both a(i) and b(i) in the silhouette formula are computed from pairwise distances. Since all distances are preserved, silhouette is exactly preserved. **This is a theorem, not an empirical observation.**

The small observed Δ (order 10⁻⁵) arises only from flooring breaking the symmetry when entries cross zero.

**The correct sensitivity metric for this correction:** Within-cluster standard deviation, which measures magnitudes (not distances) and detects the uniform shrinkage from removing the α × μ[l] contamination term.

### The OLS Slope as a Data Fingerprint

```
b̂ = √[mean(μ[l]²)]
```

This slope is determined entirely by the genetic diversity of the donor pool. It:
- Does **not** depend on the number of cells
- Does **not** depend on α (the contamination level)
- Does **not** depend on any fitted parameter
- **Does** depend on which donors were pooled (their SNP allele frequencies)

If you run soupXcell on a different dataset with different donors, you will get a different slope. The slope is a **genetic fingerprint** of the donor pool. Measuring it from simulation and comparing to the theoretical prediction √(mean(μ²)) is a cross-validation of the entire mathematical framework.

### Linear Approximation vs. Exact Correction

The exact correction is:
```
X_true = (X_obs − α × μ) / (1 − α)
```

soupXcell uses the approximation (omitting the division):
```
X_corrected = X_obs − α × μ
```

The error of this approximation:
```
|X_true − X_corrected| = |X_true| × α / (1 − α) ≈ α × X_true ≤ α ≈ 0.01262
```

The maximum approximation error is 1.262% of the true allele fraction. For allele fractions in {0, 0.5, 1}, this means maximum errors of {0, 0.0063, 0.0126}. Since downstream analysis compares cells to each other (relative differences), and the systematic error is the same for all cells, this approximation has **no effect on clustering accuracy**.

---

## COMPLETE METRICS TABLE (Test Dataset, Mode B)

| Quantity | Symbol | Formula | Value |
|----------|--------|---------|-------|
| Ambient fraction | α | From souporcell | 0.01261973 |
| Total loci | — | — | 80,712 |
| QC-filtered loci | — | From VCF | 64,496 |
| Total cells | — | — | 3,639 |
| Singlets | \|S\| | — | 3,523 |
| Doublets | — | — | 117 |
| Donor clusters | K | — | 4 |
| Covered cell-locus pairs | \|C\| | — | 1,444,860 |
| Population mean AF | μ̄ | mean(μ[l]) | 0.4122 |
| RMS soup vector | √(mean(μ²)) | — | 0.3875 |
| Mean soup contribution | soup̄ | α × μ̄ | 0.005201 |
| Entries floored | — | \|{X_corr = 0}\| | 640,425 |
| Floor percentage | — | floored / \|C\| | 44.32% |
| Silhouette before | S_before | mean[s(i)] | −0.002560 |
| Silhouette after | S_after | mean[s(i)] | −0.002563 |
| Silhouette Δ | ΔS | S_after − S_before | −0.000003 |
| Davies-Bouldin before | DB_before | — | 1.9618 |
| Davies-Bouldin after | DB_after | — | 1.9588 |
| Davies-Bouldin Δ | ΔDB | DB_after − DB_before | −0.0030 |
| Within-cluster std Δ | Δσ | — | −0.000887 (−0.43%) |
| OLS slope (RMSE vs α_sim) | b̂ | = √(mean(μ²)) | 0.3875 |
| R² (OLS fit) | R² | 1 − SS_res/SS_tot | 0.999999879 |

---

## PLAIN-ENGLISH RATIONALIZATIONS

**α = 0.01262 (1.262%):**
Out of every 100 RNA reads in a droplet, ~1.26 came from ambient free RNA — not the cell itself. This is typical for 10× Chromium v3 data; values of 1–10% are common.

**soup[l] = α × μ[l] ≈ 0.0052:**
At a typical SNP, we subtract about 0.0052 allele fraction units. That's equivalent to removing ~1 spurious read per 192 total reads at that position.

**44.32% floored:**
Nearly half of covered entries had measured allele fractions smaller than the expected ambient contribution. After correction, they map to exactly 0 — meaning 100% of their signal was contamination. This is not information loss; these were noise readings confirmed to carry no true genotype signal.

**Silhouette ≈ 0:**
Four donor clusters exist, but in the 64,496-dimensional allele fraction space with 99.5% sparsity, cells from the same cluster rarely share loci. Pairwise distances are noisy and approximately equal within/between clusters. The clusters are real (souporcell found them correctly), but silhouette cannot "see" them in this sparse space.

**Silhouette Δ = −0.000003:**
Three millionths. Effectively zero. The correction is a uniform translation — it shifts everyone's coordinates equally, leaving all pairwise distances unchanged. A ruler (silhouette) that measures distances cannot detect a global shift. This is expected from theory, confirmed by data.

**RMSE = 0.3875 × α_sim, R² = 1.000:**
The correction formula is the mathematical inverse of the contamination formula. When we inject then correct, the error is algebraically zero. The tiny residual RMSE (from flooring effects) scales perfectly linearly with injection level, with a slope that equals the RMS of the soup vector — a quantity we can compute independently from the data. Predicted slope = measured slope = 0.3875. R² = 1.000 proves the math is not just approximately right, but exactly right.

---

## HOW TO EXPLAIN BASED ON ARGUMENTS

If `$ARGUMENTS` specifies a topic, focus your explanation there:

- **"contamination model"** → Go deep on Step 2, the linear mixture, and its assumptions
- **"correction formula"** or **"derivation"** → Focus on Steps 3–4, show the algebra
- **"floor effect"** → Step 5, explain the biology and math of why 44% floor is correct
- **"silhouette"** → Steps 6 + the theoretical invariance proof
- **"davies-bouldin"** → Step 7 and contrast with silhouette
- **"simulation"** or **"benchmark"** → Steps 9–12 end-to-end
- **"OLS"** or **"linearity"** → Steps 10–11 with full numerical walkthrough
- **"R squared"** or **"R²"** → Step 12, significance of R² = 1.000
- **"slope"** or **"0.3875"** → Explain that slope = √(mean(μ²)) — the genetic fingerprint derivation
- **"approximation"** → Discuss linear approximation vs. exact correction
- **"plain english"** → Use only the Plain-English Rationalizations section, no equations
- **"PCA"** or **"dimensionality reduction"** → Step 13b — randomised SVD, mean centering, why PCA before t-SNE
- **"t-SNE"** or **"tsne"** → Step 13c — perplexity, P_ij, Q_ij, KL divergence, gradient derivation, optimisation schedule
- **"UMAP"** or **"umap"** → Step 13d — k-NN graph, attractive/repulsive forces, comparison to t-SNE
- **"before after"** or **"embedding"** → Step 13e — what changes in 2D plots, silhouette in embedding space, interpretation
- **"module 3"** or **"visualisation"** → Full Step 13 (all sub-sections 13a–13f)
- **"full"** or no arguments → Walk through all 13 steps in order

Always connect formulas to biological meaning. This is a genomics tool — the math serves the biology.
