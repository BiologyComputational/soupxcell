# soupXcell v1.0 — Complete Step-by-Step Execution Guide

> **Purpose:** Every command you need to run, in exact order, from a fresh machine
> to a fully analysed output. Copy-paste ready. Nothing is assumed.

---

## Module map — what each module does and which flags control it

```
┌─────────────────────────────────────────────────────────────────────────────┐
│  MODULE 1 — Ambient Subtraction                          ALWAYS runs        │
│  Source: src/core/corrector.rs                                              │
│  What: X_corrected[l,n] = max(0, X[l,n] − α×μ[l])                         │
│  Flags: --ref --alt --clusters --ambient --alpha --singlets_only            │
│  Flags: --souporcell_vcf --freebayes_vcf --min_ref --min_alt               │
│  Outputs: X_corrected.mtx, X_raw.mtx, soup_vector.csv,                     │
│           correction_summary.csv, metrics.txt                              │
│  Figures: silhouette_comparison.svg  correction_distribution.svg           │
│           soup_vector_profile.svg    per_cluster_correction.svg            │
│           floor_per_locus.svg                                               │
├─────────────────────────────────────────────────────────────────────────────┤
│  MODULE 2 — Simulation Benchmark                         --simulate flag    │
│  Source: src/core/simulator.rs                                              │
│  What: inject α_sim×μ[l] → recover → measure RMSE + silhouette Δ          │
│  Flags: --simulate --sim_levels --sim_trials                                │
│  Outputs: simulation_benchmark.csv                                          │
│  Figures: simulation_curves.svg    rmse_linearity.svg                      │
├─────────────────────────────────────────────────────────────────────────────┤
│  MODULE 3 — Dimensionality Reduction                     --embed flag       │
│  Source: src/core/embedder.rs                                               │
│  What: PCA → t-SNE/UMAP before + after correction                          │
│  Flags: --embed --pca_components --tsne_perplexity --tsne_iter             │
│  Outputs: (embedded in report)                                              │
│  Figures: tsne_before_after.svg    umap_before_after.svg                   │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## PART 0 — One-time setup (skip if already done)

### 0.1 Install Rust

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
# When prompted: choose option 1 (default install)
# Reload your shell
source $HOME/.cargo/env
```

### 0.2 Verify Rust version ≥ 1.75

```bash
rustc --version
# Must print: rustc 1.75.0 or higher
# If older: rustup update stable
```

### 0.3 macOS only — install Xcode Command Line Tools (linker)

```bash
xcode-select --install
# Click "Install" in the dialog — only command-line tools are needed
# Verify the linker is available
clang --version
```

---

## PART 1 — Get the source and build

### 1.1 Unzip the source archive

```bash
cd ~/Setups/ambientRNA_research/souporcell_pipeline_components
unzip soupxcell_v1.0.zip
cd soupxcell
```

### 1.2 Confirm the source tree is complete

```bash
ls
# Must show: Cargo.toml  .soupxcell.env  profiles/  src/
find src -name '*.rs' | wc -l
# Must print ≥ 20
ls src/analysis/report/
# Must show 11 .rs files: mod.rs  cards.rs  css.rs  experiment.rs ...
wc -l src/analysis/plots.rs
# Must print ≥ 1174  (contains all 9 plot functions)
wc -l src/core/corrector.rs src/core/simulator.rs src/core/embedder.rs
# Module 1, 2, 3 source files — each must be > 0 lines
```

### 1.3 Build the release binary

```bash
cargo build --release 2>&1 | tee build.log
# First build: 3–8 min (downloads ~106 crates)
# Subsequent builds: < 30s (incremental)
```

### 1.4 Confirm the build succeeded

```bash
tail -3 build.log
# Must contain: Finished release [optimized] target(s) in ...
# If you see errors: check PART 0 and re-run
```

### 1.5 Verify the binary runs

```bash
./target/release/soupxcell --version
# Must print: soupxcell 1.0.0
./target/release/soupxcell --help 2>&1 | head -5
# Must print the help banner
```

### 1.6 Set the BIN shorthand

```bash
BIN=./target/release/soupxcell
# Every run command below uses $BIN
```

---

## PART 2 — Set input file paths

### 2.1 Set all shell variables

```bash
REF=/Volumes/extra_space/demux_test/ref.mtx        # Module 1 input
ALT=/Volumes/extra_space/demux_test/alt.mtx        # Module 1 input
CLU=/Volumes/extra_space/demux_test/clusters.tsv   # Module 1 input
AMB=/Volumes/extra_space/demux_test/ambient_rna.txt  # Module 1 input (α)
VCF=/Volumes/extra_space/demux_test/cluster_genotypes.vcf  # Module 1 locus list
OUT=~/Setups/ambientRNA_research/souporcell_pipeline_components/playground/output_soupXCell
```

### 2.2 Create output directory

```bash
mkdir -p $OUT
```

### 2.3 Verify every input file exists

```bash
for f in $REF $ALT $CLU $AMB $VCF; do
  [ -f "$f" ] && echo "OK       $f" || echo "MISSING  $f"
done
# All five lines must print OK
# If MISSING: volume may not be mounted → ls /Volumes/extra_space/demux_test/
```

### 2.4 Verify file contents are readable

```bash
# Module 1 — check ambient fraction (α)
cat $AMB
# Must print a float: 0.01261973 or "ambient RNA estimated as 1.26%"
# Module 1 — check clusters.tsv has singlet/doublet status column
head -3 $CLU
# Must have columns: barcode, status, cluster_id
# status values must include "singlet" and "doublet"
# Module 1 — check matrices are non-empty
wc -l $REF $ALT
# Both should have > 80000 lines
# Module 1 — check VCF locus list is complete
wc -l $VCF
# Should print > 8000 lines
```

---

## PART 3 — Configure the environment file (one-time)

The `.soupxcell.env` file sets defaults for all three modules. Edit it once.

### 3.1 Open and edit

```bash
nano .soupxcell.env
```

### 3.2 Set these keys (covers all three modules)

```ini
# ── Module 1 inputs ───────────────────────────────────────────────────────────
SOUPX_REF=/Volumes/extra_space/demux_test/ref.mtx
SOUPX_ALT=/Volumes/extra_space/demux_test/alt.mtx
SOUPX_CLUSTERS=/Volumes/extra_space/demux_test/clusters.tsv
SOUPX_AMBIENT=/Volumes/extra_space/demux_test/ambient_rna.txt

# ── Module 1 locus selection ──────────────────────────────────────────────────
SOUPX_SOUPORCELL_VCF=/Volumes/extra_space/demux_test/cluster_genotypes.vcf

# ── Module 1 correction settings ─────────────────────────────────────────────
SOUPX_THREADS=8
SOUPX_SINGLETS_ONLY=true
SOUPX_SEED=42
# SOUPX_ALPHA=0.05    # uncomment only to override ambient_rna.txt

# ── Module 3 embedding settings ───────────────────────────────────────────────
SOUPX_EMBED=tsne
SOUPX_PCA_COMPONENTS=50
SOUPX_TSNE_PERPLEXITY=30.0
SOUPX_TSNE_ITER=1000

# ── Module 2 simulation settings ──────────────────────────────────────────────
SOUPX_SIMULATE=false
SOUPX_SIM_LEVELS=0.01,0.05,0.10,0.20,0.30,0.50
SOUPX_SIM_TRIALS=3

# ── Output ────────────────────────────────────────────────────────────────────
SOUPX_OUTPUT=~/Setups/ambientRNA_research/souporcell_pipeline_components/playground/output_soupXCell
```

> Save: `Ctrl+O`, `Enter`, `Ctrl+X`

---

## PART 4 — Dry run (preview all three modules before executing)

**Always do this before the first real run.**
Nothing is read or written until you approve.

```bash
$BIN --dry_run \
  --ref $REF --alt $ALT --clusters $CLU --ambient $AMB \
  --output $OUT \
  --souporcell_vcf $VCF \
  --embed tsne \
  --simulate \
  --threads 8
```

> The pre-flight output shows:
> - **Module 1:** input paths, α value, locus source (Mode A/B/C), correction formula
> - **Module 2:** simulation levels, trials per level
> - **Module 3:** embedding method, PCA components, t-SNE parameters
> - All expected output files and estimated runtime
>
> Type **`y`** to proceed, **`n`** to abort.

---

## PART 5 — Smoke test (verify all three modules work end-to-end, ~15s)

This runs a reduced version of all three modules to confirm the build is healthy.

```bash
$BIN \
  --ref $REF --alt $ALT --clusters $CLU --ambient $AMB \
  --output $OUT/smoke_test \
  --souporcell_vcf $VCF \
  --embed tsne \
  --pca_components 10 \
  --tsne_iter 100 \
  --threads 4 \
  --seed 42
```

> Note: `--simulate` is not included — Module 2 is slow and not needed for a smoke test.

### 5.1 Verify all three modules ran

```bash
ls $OUT/smoke_test/
# Must contain: X_corrected.mtx  metrics.txt  soupxcell_report.html  figures/
# Module 1 — check correction result
cat $OUT/smoke_test/metrics.txt
# Must show: alpha_used = 0.01261973  entries_floored = ...
# Module 1 — check the 5 always-generated figures exist
ls $OUT/smoke_test/figures/
# Must contain all five Module 1 figures:
#   silhouette_comparison.svg       ← Module 1: cluster quality metric
#   correction_distribution.svg     ← Module 1: per-cell correction histogram
#   soup_vector_profile.svg         ← Module 1: locus contamination profile
#   per_cluster_correction.svg      ← Module 1: per-cluster summary
#   floor_per_locus.svg             ← Module 1: floor fraction per locus
#   tsne_before_after.svg           ← Module 3: before/after scatter
```

> If `metrics.txt` exists and `tsne_before_after.svg` exists, all modules
> ran successfully. Proceed to Part 6.

---

## PART 6 — Run modes by module combination

---

### RUN A — MODULE 1 ONLY: correction without embedding (~2s)

**Modules active: Module 1 only**
`--embed none` disables Module 3. `--simulate` is absent, so Module 2 is also off.
Use this for a quick correction pass, downstream analysis, or CI pipelines.

```bash
$BIN \
  --ref $REF --alt $ALT --clusters $CLU --ambient $AMB \
  --output $OUT \
  --souporcell_vcf $VCF \
  --singlets_only \
  --threads 8 \
  --embed none
```

```bash
$BIN \
  --ref $REF --alt $ALT --clusters $CLU --ambient $AMB \
  --output $OUT \
  --souporcell_vcf $VCF \
  --singlets_only \
  --threads 8 \
  --embed none
```

**Verify Module 1 output:**

```bash
# Core correction metrics
cat $OUT/metrics.txt
grep 'loci_from_vcf' $OUT/metrics.txt   # expect: 8219
grep 'entries_floored' $OUT/metrics.txt # expect: 640425
grep 'mean_correction' $OUT/metrics.txt # expect: 0.00484321
# Module 1 figures (all 5 should exist)
ls $OUT/figures/
open $OUT/figures/soup_vector_profile.svg      # which loci are most contaminated?
open $OUT/figures/per_cluster_correction.svg   # is correction uniform across clusters?
open $OUT/figures/floor_per_locus.svg          # which loci floor the most cells?
open $OUT/figures/silhouette_comparison.svg    # cluster quality before/after
open $OUT/figures/correction_distribution.svg  # per-cell correction histogram
```

**Module 1 output files:**
```
X_corrected.mtx                ← corrected allele fraction matrix [L × N]
X_raw.mtx                      ← original allele fraction matrix  [L × N]
soup_vector.csv                 ← per-locus: μ[l], soup[l]=α×μ[l], α
correction_summary.csv          ← per-cell: barcode, cluster, mean_correction
metrics.txt                     ← α, loci, cells, floor%, Δsilhouette

figures/silhouette_comparison.svg     ← MODULE 1 figure
figures/correction_distribution.svg  ← MODULE 1 figure
figures/soup_vector_profile.svg       ← MODULE 1 figure
figures/per_cluster_correction.svg    ← MODULE 1 figure
figures/floor_per_locus.svg           ← MODULE 1 figure
```

---

### RUN B — MODULES 1 + 3: correction + t-SNE before/after (~100s)

**Modules active: Module 1 + Module 3 (t-SNE)**
`--embed tsne` enables Module 3. Module 2 is still off.
This is the professor's core experiment: *"do the t-SNE before and after
subtraction to see if it tightens up the clusters."*

```bash
$BIN \
  --ref $REF --alt $ALT --clusters $CLU --ambient $AMB \
  --output $OUT \
  --souporcell_vcf $VCF \
  --singlets_only \
  --threads 8 \
  --embed tsne \
  --pca_components 50 \
  --tsne_perplexity 30 \
  --tsne_iter 1000 \
  --seed 42
```

> **Clean version:**

```bash
$BIN \
  --ref $REF --alt $ALT --clusters $CLU --ambient $AMB \
  --output $OUT \
  --souporcell_vcf $VCF \
  --singlets_only \
  --threads 8 \
  --embed tsne \
  --pca_components 50 \
  --tsne_perplexity 30 \
  --tsne_iter 1000 \
  --seed 42
```

**Verify Module 1 + Module 3 output:**

```bash
# Module 1 verification (same as Run A)
grep 'entries_floored\|silhouette_delta' $OUT/metrics.txt
# Module 3 verification — the key figure
open $OUT/figures/tsne_before_after.svg
# Left panel  = before correction (raw allele fractions)
# Right panel = after correction
# Cells coloured by cluster: cyan=0, violet=1, green=2, amber=3...
# Each panel title shows its silhouette score
# Footer: "Silhouette Δ = +0.0000" — this is EXPECTED at α=1.26%
open $OUT/soupxcell_report.html
```

**Output files added by Module 3 (on top of Run A):**
```
figures/tsne_before_after.svg         ← MODULE 3 figure (t-SNE before/after)
```

---

### RUN C — MODULES 1 + 2: correction + simulation benchmark (~300s)

**Modules active: Module 1 + Module 2**
`--simulate` enables Module 2. `--embed none` keeps Module 3 off.
This directly addresses the supervisor's comment: *"1.26% is small — we might
need to simulate more to see a difference."*

```bash
$BIN \
  --ref $REF --alt $ALT --clusters $CLU --ambient $AMB \
  --output $OUT \
  --souporcell_vcf $VCF \
  --singlets_only \
  --threads 8 \
  --simulate \
  --sim_levels 0.01 0.05 0.10 0.20 0.30 0.50 \
  --sim_trials 3 \
  --embed none
```

> **Clean version:**

```bash
$BIN \
  --ref $REF --alt $ALT --clusters $CLU --ambient $AMB \
  --output $OUT \
  --souporcell_vcf $VCF \
  --singlets_only \
  --threads 8 \
  --simulate \
  --sim_levels 0.01 0.05 0.10 0.20 0.30 0.50 \
  --sim_trials 3 \
  --embed none
```

**Verify Module 2 output:**

```bash
# Module 2 — simulation results table
column -t -s, $OUT/simulation_benchmark.csv
# Expected (key insight: RMSE = 0.3668 × α_sim, R² = 1.000):
#   alpha_sim  rmse      sil_before  sil_after  delta    floor_pct
#   0.01       0.003668  -0.0195     -0.0195    +0.0000  0.000
#   0.05       0.018340  -0.0193     -0.0195    -0.0002  0.000
#   0.10       0.036680  -0.0191     -0.0195    -0.0004  0.000
#   0.20       0.073360  -0.0188     -0.0197    -0.0009  0.000
#   0.30       0.110040  -0.0185     -0.0198    -0.0013  0.000
#   0.50       0.183400  -0.0175     -0.0201    -0.0026  0.000
# Module 2 figures
open $OUT/figures/simulation_curves.svg
# Shows: RMSE (cyan, left axis) and silhouette before/after (amber/green, right axis)
# vs α_sim. Vertical amber dashed line marks the true run α = 0.0126.
open $OUT/figures/rmse_linearity.svg
# Shows: scatter of (α_sim, RMSE) + OLS fit line + annotation:
#   slope = 0.3668  (RMSE / α)
#   R²    = 1.000000  ← perfect linear recovery
#   ✓ Perfect linearity confirmed
```

**Output files added by Module 2 (on top of Run A):**
```
simulation_benchmark.csv              ← MODULE 2 output
figures/simulation_curves.svg         ← MODULE 2 figure
figures/rmse_linearity.svg            ← MODULE 2 figure
```

---

### RUN D — ALL THREE MODULES: full experiment (~10 min)

**Modules active: Module 1 + Module 2 + Module 3**
Complete run — correction + simulation + t-SNE + UMAP.
All 9 figures are generated.

```bash
$BIN \
  --ref $REF --alt $ALT --clusters $CLU --ambient $AMB \
  --output $OUT \
  --souporcell_vcf $VCF \
  --singlets_only \
  --threads 8 \
  --simulate \
  --sim_levels 0.01 0.05 0.10 0.20 0.30 0.50 \
  --sim_trials 5 \
  --embed both \
  --pca_components 50 \
  --tsne_perplexity 30 \
  --tsne_iter 1000 \
  --seed 42
```

> **Clean version:**

```bash
$BIN \
  --ref $REF --alt $ALT --clusters $CLU --ambient $AMB \
  --output $OUT \
  --souporcell_vcf $VCF \
  --singlets_only \
  --threads 8 \
  --simulate \
  --sim_levels 0.01 0.05 0.10 0.20 0.30 0.50 \
  --sim_trials 5 \
  --embed both \
  --pca_components 50 \
  --tsne_perplexity 30 \
  --tsne_iter 1000 \
  --seed 42
```

**Using the bundled profile instead (equivalent):**

```bash
$BIN --config profiles/simulate_high.json \
  --ref $REF --alt $ALT --clusters $CLU --ambient $AMB \
  --output $OUT --souporcell_vcf $VCF
```

**Verify all three modules ran:**

```bash
# All 9 figures must exist
ls $OUT/figures/
# Module 1 figures (always generated):
#   silhouette_comparison.svg
#   correction_distribution.svg
#   soup_vector_profile.svg
#   per_cluster_correction.svg
#   floor_per_locus.svg
#
# Module 2 figures (generated because --simulate was set):
#   simulation_curves.svg
#   rmse_linearity.svg
#
# Module 3 figures (generated because --embed both was set):
#   tsne_before_after.svg
#   umap_before_after.svg
open $OUT/soupxcell_report.html
# The report embeds all 9 figures, 12 metric cards, 9 interpretation cards,
# 4 data tables, and a complete methods section with all three modules documented
```

**Complete output file list for Run D (all three modules):**
```
X_corrected.mtx                       ← MODULE 1 output
X_raw.mtx                             ← MODULE 1 output
soup_vector.csv                        ← MODULE 1 output
correction_summary.csv                 ← MODULE 1 output
metrics.txt                            ← MODULE 1 output
simulation_benchmark.csv              ← MODULE 2 output
soupxcell_report.html                  ← all modules combined

figures/silhouette_comparison.svg      ← MODULE 1 figure
figures/correction_distribution.svg   ← MODULE 1 figure
figures/soup_vector_profile.svg        ← MODULE 1 figure
figures/per_cluster_correction.svg     ← MODULE 1 figure
figures/floor_per_locus.svg            ← MODULE 1 figure
figures/simulation_curves.svg          ← MODULE 2 figure
figures/rmse_linearity.svg             ← MODULE 2 figure
figures/tsne_before_after.svg          ← MODULE 3 figure
figures/umap_before_after.svg          ← MODULE 3 figure
```

---

### RUN E — MODULES 1 + 3: UMAP variant

**Modules active: Module 1 + Module 3 (UMAP only)**

```bash
$BIN \
  --ref $REF --alt $ALT --clusters $CLU --ambient $AMB \
  --output $OUT \
  --souporcell_vcf $VCF \
  --singlets_only \
  --threads 8 \
  --embed umap \
  --pca_components 50 \
  --seed 42
```

**Output:** `figures/umap_before_after.svg` (Module 3) + all 5 Module 1 figures.

---

### RUN F — MODULE 1: alpha override (compare contamination levels)

**Modules active: Module 1 only (at a user-specified α)**
Use to compare the correction at the measured α=1.26% vs a hypothetical level.

```bash
$BIN \
  --ref $REF --alt $ALT --clusters $CLU --ambient $AMB \
  --output ${OUT}_alpha5pct \
  --souporcell_vcf $VCF \
  --alpha 0.05 \
  --singlets_only \
  --threads 8 \
  --embed tsne \
  --seed 42
open ${OUT}_alpha5pct/figures/tsne_before_after.svg
# At α=5%, the left/right t-SNE panels should look more different than at α=1.26%
```

---

### RUN G — MODULE 1 + VCF Mode A (maximum locus reproducibility)

**Modules active: Module 1 (with exact CHROM:POS locus matching)**
Requires both `--souporcell_vcf` (souporcell VCF) and `--freebayes_vcf`
(freebayes VCF). Handles chr/no-chr naming mismatches automatically.

```bash
FB_VCF=/Volumes/extra_space/demux_test/souporcell_merged_sorted_vcf.vcf.gz
$BIN \
  --ref $REF --alt $ALT --clusters $CLU --ambient $AMB \
  --output $OUT \
  --souporcell_vcf $VCF \
  --freebayes_vcf  $FB_VCF \
  --singlets_only \
  --threads 8 \
  --embed tsne \
  --seed 42
# Log will show: mode = A — exact CHROM:POS cross-reference
```

---

### RUN H — Web GUI (all modules via browser)

```bash
$BIN --gui
# Open: http://localhost:7878
```

The GUI exposes all three modules:
- **Configure phase** — Module 1/2/3 params on separate cards, live Mode A/B/C locus indicator
- **Pre-flight phase** — module-by-module plan overview
- **Execute phase** — live log with per-module progress, Module 2 simulation table
- **Outputs phase** — figures grouped by module, download links

---

## PART 7 — Post-run verification by module

### 7.1 Module 1 verification

```bash
# Core correction worked
cat $OUT/metrics.txt
# Expected values
grep 'alpha_used' $OUT/metrics.txt      # 0.01261973
grep 'loci_from_vcf\|loci_qc' $OUT/metrics.txt  # 8219
grep 'entries_floored' $OUT/metrics.txt # 640425
grep 'pct_floored' $OUT/metrics.txt     # 44.30%
grep 'mean_correction' $OUT/metrics.txt # 0.00484321
grep 'silhouette_delta' $OUT/metrics.txt # +0.00000 (expected at α=1.26%)
grep 'within_std_delta' $OUT/metrics.txt # -0.00235 (the real signal)
# Soup vector — which loci are most contaminated?
sort -t, -k3 -rn $OUT/soup_vector.csv | head -6
# Top loci by soup[l]=α×μ[l]
# Per-cell corrections — which cells were most corrected?
sort -t, -k4 -rn $OUT/correction_summary.csv | head -6
```

### 7.2 Module 2 verification

```bash
# Simulation results table
column -t -s, $OUT/simulation_benchmark.csv
# Key: RMSE should be proportional to alpha_sim (R²=1.000)
# Key: silhouette delta ≈ 0 at low α, becomes visible at α ≥ 0.10
# Key figures
open $OUT/figures/rmse_linearity.svg
# Confirm: slope ≈ 0.3668, R² = 1.000, residuals ≈ 0
open $OUT/figures/simulation_curves.svg
# Confirm: RMSE line is straight, silhouette lines are nearly flat at low α
```

### 7.3 Module 3 verification

```bash
# t-SNE/UMAP figures
open $OUT/figures/tsne_before_after.svg
# Expected: left and right panels look near-identical at α=1.26%
# This is CORRECT — the effect is small at real α
open $OUT/figures/umap_before_after.svg
# Same expectation — near-identical at true α
# Silhouette bar chart (Module 1 computes these, Module 3 plots them per-cluster)
open $OUT/figures/silhouette_comparison.svg
# Expected: bars nearly equal height before and after
```

### 7.4 Open all 9 figures at once

```bash
# Module 1 figures
open $OUT/figures/silhouette_comparison.svg
open $OUT/figures/correction_distribution.svg
open $OUT/figures/soup_vector_profile.svg
open $OUT/figures/per_cluster_correction.svg
open $OUT/figures/floor_per_locus.svg
# Module 2 figures
open $OUT/figures/simulation_curves.svg
open $OUT/figures/rmse_linearity.svg
# Module 3 figures
open $OUT/figures/tsne_before_after.svg
open $OUT/figures/umap_before_after.svg
```

### 7.5 Open the HTML report (all modules combined)

```bash
open $OUT/soupxcell_report.html
# §2  Metric cards — 12 cards covering Modules 1 + 3 results
# §3  Experiment design — locus source (Module 1), correction formula
# §4  Correction formula detail — Module 1 step-by-step
# §5  Tables — per-cluster (Module 1), metrics before/after (Module 3),
#              simulation (Module 2), locus QC (Module 1)
# §6  Inline SVGs — all 9 figures from Modules 1, 2, 3
# §8  Interpretation cards — 9 cards interpreting Module 1, 2, 3 results
# §9  Parameters — 27 resolved params with module labels
# §10 Methods — Module 1/2/3 algorithms and references
```

---

## PART 8 — Flag-to-module quick reference

| Flag | Module | Effect |
|---|---|---|
| `--ref` `--alt` | Module 1 | Raw allele count matrix inputs |
| `--clusters` | Module 1 | Cell barcode → cluster → singlet/doublet |
| `--ambient` | Module 1 | Source of α (ambient fraction) |
| `--alpha` | Module 1 | Override α from ambient_rna.txt |
| `--singlets_only` | Module 1 | Use only singlet cells to estimate μ[l] |
| `--souporcell_vcf` | Module 1 | Locus list (Mode B) |
| `--freebayes_vcf` | Module 1 | Upgrades locus matching to Mode A |
| `--min_ref` `--min_alt` | Module 1 | Locus QC fallback (Mode C) |
| `--simulate` | Module 2 | Enable simulation benchmark |
| `--sim_levels` | Module 2 | α values to inject and recover |
| `--sim_trials` | Module 2 | Trials per level |
| `--embed` | Module 3 | `tsne` \| `umap` \| `both` \| `none` |
| `--pca_components` | Module 3 | PCA dims before t-SNE/UMAP |
| `--tsne_perplexity` | Module 3 | t-SNE neighbourhood size |
| `--tsne_iter` | Module 3 | t-SNE optimisation iterations |
| `--seed` | Module 3 | RNG seed for PCA + t-SNE + UMAP |
| `--threads` | Modules 1+2+3 | Parallel threads (Rayon) |
| `--output` | All | Output directory |

---

## PART 9 — Troubleshooting

| Symptom | Module | Command to diagnose | Fix |
|---|---|---|---|
| `cargo build` linker error | — | `clang --version` | `xcode-select --install` |
| `rustc` too old | — | `rustc --version` | `rustup update stable` |
| `File not found: ref.mtx` | Module 1 | `ls /Volumes/extra_space/demux_test/` | Mount the volume |
| `0 loci from VCF` | Module 1 | `wc -l $VCF` | VCF must have >8000 lines |
| `n_cells_singlet = 0` | Module 1 | `head $CLU` | clusters.tsv needs `status` column |
| `alpha = 0.0000` | Module 1 | `cat $AMB` | Must contain float `0.01261973` |
| Silhouette Δ ≈ 0 | Module 1/3 | Expected | Run Module 2 to demonstrate at higher α |
| Floor% > 80% | Module 1 | `grep alpha $OUT/metrics.txt` | α may be from wrong run |
| Simulation RMSE not linear | Module 2 | `column -t -s, $OUT/simulation_benchmark.csv` | Need ≥3 levels spanning 10× range |
| t-SNE looks noisy | Module 3 | Re-run | `--tsne_iter 3000 --tsne_perplexity 50` |
| UMAP not generated | Module 3 | `ls $OUT/figures/` | Must use `--embed umap` or `--embed both` |
| Out of memory | Module 3 | — | `--embed none` first (saves ~300 MB) |
| Figures missing in report | All | `ls $OUT/figures/` | Open HTML from `$OUT/` directory directly |

---

## PART 10 — Expected numeric results (all three modules)

All values confirmed on `/Volumes/extra_space/demux_test/` with souporcell v3.0.

**Module 1 results:**

| Metric | Expected value |
|---|---|
| Matrix dimensions | 80,712 loci × 3,639 cells |
| Loci from VCF (Mode B) | 8,219 |
| α | 0.01261973 (1.26%) |
| Entries floored at 0 | 640,425 (44.3% of covered) |
| Mean \|correction\| | 0.00484321 |
| soup_max | 0.01261973 (= α × max μ[l]) |
| soup_mean | 0.00484321 (= α × mean μ[l]) |
| Silhouette before | −0.01950 |
| Silhouette after | −0.01950 (Δ ≈ 0 — expected) |
| Within-std Δ | −0.00235 (the real signal from Module 1) |

**Module 2 results:**

| α_sim | RMSE | Silhouette Δ |
|---|---|---|
| 0.01 | 0.003668 | ≈ 0 |
| 0.05 | 0.018340 | −0.0002 |
| 0.10 | 0.036680 | −0.0004 |
| 0.20 | 0.073360 | −0.0009 |
| 0.30 | 0.110040 | −0.0013 |
| 0.50 | 0.183400 | −0.0026 |
| **Pattern** | **RMSE = 0.3668 × α_sim** | **R² = 1.000** |

**Module 3 results:**

| Metric | Expected value |
|---|---|
| Silhouette before (t-SNE) | −0.01950 |
| Silhouette after (t-SNE) | −0.01950 (Δ ≈ 0) |
| t-SNE runtime | ~42s per embedding (1000 iters, 3522 cells) |
| PCA runtime | ~0.2s |

**Runtimes:**

| Run | Modules active | Expected time |
|---|---|---|
| Smoke test | 1 + 3 (fast) | ~15s |
| Run A | 1 only | ~2s |
| Run B | 1 + 3 (t-SNE) | ~100s |
| Run C | 1 + 2 | ~300s |
| Run D | 1 + 2 + 3 (both embeds) | ~10 min |

---

*soupXcell v1.0 · Jahidul Arafat*
