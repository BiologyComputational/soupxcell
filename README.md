# soupXcell v1.0 — Ambient RNA Correction for Genotype Demultiplexing

<p align="center">
  <strong>Post-processor for the souporcell pipeline — correct, visualise, and benchmark ambient RNA contamination</strong>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/version-1.0.0-8B5CF6" alt="version"/>
  <img src="https://img.shields.io/badge/rust-≥1.75-orange" alt="rust"/>
  <img src="https://img.shields.io/badge/edition-2021-blue" alt="edition"/>
  <img src="https://img.shields.io/badge/figures-9%20SVG%20plots-10B981" alt="figures"/>
</p>

---

> **souporcell pipeline:** Heaton et al., *Nature Methods* 2020 · [doi:10.1038/s41592-020-0820-1](https://doi.org/10.1038/s41592-020-0820-1)
>
> **soupXcell:** Jahidul Arafat — ambient RNA correction, simulation benchmark, before/after t-SNE/UMAP visualisation

---

## Table of Contents

1. [What soupXcell does](#1-what-soupxcell-does)
2. [Where it fits in the pipeline](#2-where-it-fits-in-the-pipeline)
3. [The biological problem and correction formula](#3-the-biological-problem-and-correction-formula)
4. [Repository layout](#4-repository-layout)
5. [Prerequisites and building](#5-prerequisites-and-building)
6. [Quick start](#6-quick-start)
7. [All run modes](#7-all-run-modes)
8. [Complete parameter reference](#8-complete-parameter-reference)
9. [Module 1 — Ambient Subtraction](#9-module-1--ambient-subtraction)
10. [Module 2 — Simulation Benchmark](#10-module-2--simulation-benchmark)
11. [Module 3 — Dimensionality Reduction](#11-module-3--dimensionality-reduction)
12. [Locus selection — VCF modes A / B / C](#12-locus-selection--vcf-modes-a--b--c)
13. [Step-by-step run commands](#13-step-by-step-run-commands)
14. [All 9 output figures](#14-all-9-output-figures)
15. [Complete output file reference](#15-complete-output-file-reference)
16. [Config profiles (.json)](#16-config-profiles-json)
17. [Environment file (.soupxcell.env)](#17-environment-file-soupxcellenv)
18. [Reading the logs](#18-reading-the-logs)
19. [Web GUI](#19-web-gui)
20. [Troubleshooting](#20-troubleshooting)

---

## 1. What soupXcell does

soupXcell is a **post-processor** for the souporcell pipeline. It takes the
allele count matrices and cluster assignments produced by the pipeline and:

1. **Corrects** each cell's allele fraction measurements by subtracting the
   estimated ambient RNA contribution (Module 1)
2. **Benchmarks** the correction at multiple contamination levels using synthetic
   injection — directly addressing the supervisor's concern that 1.26% is too
   small to see a visible effect (Module 2)
3. **Visualises** the effect via t-SNE and/or UMAP embeddings with **9 diagnostic
   SVG plots** and a complete 10-section HTML report (Module 3)

Runs as a **CLI tool** (with `--dry_run` pre-flight), as a **web GUI** (`--gui`),
and produces a self-contained HTML report plus 9 SVG figures.

---

## 2. Where it fits in the pipeline

```
souporcell_pipeline_v2.py
  Stage 1   samtools merge      Merge BAM files
  Stage 2   minimap2            Remap reads
  Stage 3   retag.py            Restore CB/UMI tags
  Stage 4   freebayes           Call SNVs
  Stage 5   vartrix             Count alleles → ref.mtx, alt.mtx
  Stage 6   souporcell          Cluster cells → clusters.tsv
  Stage 7   troublet            Flag doublets → clusters.tsv (status column)
  Stage 8   consensus.py        → ambient_rna.txt
                                → cluster_genotypes.vcf  ← definitive locus list
                                                                    │
                                             ┌──────────────────────┘
                                             ▼
                                       soupXcell v1.0
                                  (this tool — takes all outputs)
                                             │
                                             ▼
                              X_corrected.mtx, soup_vector.csv
                              9 SVG figures, metrics.txt
                              soupxcell_report.html
```

All four required inputs already exist after a full pipeline run.

---

## 3. The biological problem and correction formula

### Why ambient RNA contaminates allele fractions

When cells are captured in droplets, free ambient mRNA (the "soup") from
lysed cells is also captured. In a multiplexed experiment, this ambient RNA
carries alleles from **all donors in the pool**, pulling every cell's allele
fractions toward the population-wide mean.

`consensus.py` estimates the ambient fraction α (test dataset: α = 0.012619,
i.e. 1.26%). soupXcell uses that estimate to correct the allele fraction matrices.

### The correction formula (supervisor's exact specification)

```
ref.mtx    →  R[l,n]   ref allele UMI count at locus l, cell n
alt.mtx    →  A[l,n]   alt allele UMI count at locus l, cell n
α          =  0.01261973  (from ambient_rna.txt)

Step 1:  X[l,n] = A[l,n] / (R[l,n] + A[l,n])    NaN where total = 0
Step 2:  μ[l]   = mean_n( X[l,n] )  over singlet cells
         soup[l] = α × μ[l]
Step 3:  X_corrected[l,n] = max( 0,  X[l,n] − soup[l] )
```

The `max(0, ...)` floor prevents negative allele fractions (supervisor's
explicit requirement).

---

## 4. Repository layout

```
soupxcell/
├── Cargo.toml                        — 18 dependencies
├── .soupxcell.env                    — copy to CWD, edit paths
├── profiles/
│   ├── default.json                  — standard run (t-SNE only)
│   ├── fast_test.json                — smoke test (100 iterations)
│   └── simulate_high.json           — full benchmark α=[0.01…1.0], both embeds
└── src/
    ├── main.rs                731L   — composition root
    ├── domain/
    │   ├── types.rs           200L   — all structs (CorrectionResult, RunSummary, etc.)
    │   └── math.rs            337L   — allele_fraction, soup_vector, subtract_ambient,
    │                                   silhouette, davies_bouldin, inject_contamination
    ├── config/
    │   ├── params.rs          350L   — Params struct, 3-layer resolve, expand_tilde
    │   ├── params.yml         155L   — all CLI flags with help text
    │   ├── env_loader.rs       59L   — SOUPX_* reader
    │   └── config_loader.rs    39L   — JSON profile reader
    ├── core/
    │   ├── corrector.rs        86L   — MODULE 1: run_correction()
    │   ├── simulator.rs       115L   — MODULE 2: run_simulation()
    │   └── embedder.rs        340L   — MODULE 3: run_pca, run_tsne, run_umap
    ├── infra/
    │   ├── io.rs              380L   — load_raw_matrices, load_matrices_qc,
    │   │                               load_clusters, parse_ambient_fraction,
    │   │                               parse_vcf_locus_indices,
    │   │                               write_corrected_mtx, write_soup_vector
    │   ├── logger.rs           88L   — structured stderr, banner
    │   ├── preflight.rs       122L   — --dry_run plan + [Y/n] gate
    │   ├── gui_server.rs      465L   — --gui HTTP server
    │   └── soupxcell_gui.html 1321L  — web GUI (4-phase mission control)
    └── analysis/
        ├── plots.rs          1174L   — 9 SVG plot functions
        └── report/           1805L   — 11-module HTML report system
            ├── mod.rs                  entry point, ReportData struct
            ├── css.rs, utils.rs        stylesheet, colour tokens
            ├── cards.rs                12 metric cards
            ├── experiment.rs           experiment design section
            ├── tables.rs               4 data tables
            ├── svgs.rs, plot_embeds.rs figure embeds
            ├── interpretation.rs       9 PhD-level interpretation cards
            ├── params.rs               27-parameter table
            └── methods.rs              methods + references
```

**Dependency rule:** `domain ← config ← core ← infra ← analysis ← main.rs`.
No inner layer imports from an outer layer. Algorithms never write to disk.

---

## 5. Prerequisites and building

| Tool | Version | Notes |
|---|---|---|
| Rust + Cargo | ≥ 1.75 (edition 2021) | [rustup.rs](https://rustup.rs) |
| C linker | system default | macOS: `xcode-select --install` |

```bash
# Install Rust (one-time — skip if rustc --version ≥ 1.75 already)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env
rustc --version    # must be ≥ 1.75.0

# macOS Apple Silicon — linker prerequisite
xcode-select --install

# Unzip and build
cd ~/Setups/ambientRNA_research/souporcell_pipeline_components
unzip soupxcell_v1.0.zip
cd soupxcell

cargo build --release 2>&1 | tee build.log
# First build: 3–8 min. Look for: Finished release [optimized]

# Verify
./target/release/soupxcell --version    # expect: soupxcell 1.0.0
```

---

## 6. Quick start

```bash
BIN=./target/release/soupxcell
REF=/Volumes/extra_space/demux_test/ref.mtx
ALT=/Volumes/extra_space/demux_test/alt.mtx
CLU=/Volumes/extra_space/demux_test/clusters.tsv
AMB=/Volumes/extra_space/demux_test/ambient_rna.txt
VCF=/Volumes/extra_space/demux_test/cluster_genotypes.vcf
OUT=~/Setups/ambientRNA_research/souporcell_pipeline_components/playground/output_soupXCell
mkdir -p $OUT

# Smoke test (~15s)
$BIN --ref $REF --alt $ALT --clusters $CLU --ambient $AMB \
     --output $OUT/smoke \
     --embed tsne --pca_components 10 --tsne_iter 100 --threads 4

# Standard run (~100s)
$BIN --ref $REF --alt $ALT --clusters $CLU --ambient $AMB \
     --output $OUT --souporcell_vcf $VCF \
     --embed tsne --pca_components 50 --tsne_perplexity 30 --tsne_iter 1000 \
     --threads 8 --seed 42

open $OUT/soupxcell_report.html
open $OUT/figures/tsne_before_after.svg
```

---

## 7. All run modes

**Priority chain:** `CLI flags > --config JSON > .soupxcell.env > built-in defaults`

```bash
# Mode 1: Full CLI
./target/release/soupxcell \
  --ref $REF --alt $ALT --clusters $CLU --ambient $AMB \
  --output $OUT --souporcell_vcf $VCF --threads 8 --seed 42

# Mode 2: JSON profile + CLI overrides
./target/release/soupxcell --config profiles/simulate_high.json \
  --ref $REF --alt $ALT --clusters $CLU --ambient $AMB \
  --output $OUT --souporcell_vcf $VCF

# Mode 3: .env file (auto-discovered in CWD)
cp .soupxcell.env /your/working/directory/
cd /your/working/directory/
./target/release/soupxcell

# Mode 4: Dry run pre-flight
./target/release/soupxcell --dry_run \
  --ref $REF --alt $ALT --clusters $CLU --ambient $AMB \
  --output $OUT --embed tsne --simulate

# Mode 5: Web GUI
./target/release/soupxcell --gui              # http://localhost:7878
./target/release/soupxcell --gui --gui_port 8080
```

---

## 8. Complete parameter reference

### Required inputs

| Flag | Env key | Description |
|---|---|---|
| `--ref <path>` | `SOUPX_REF` | ref.mtx from vartrix (Stage 5) |
| `--alt <path>` | `SOUPX_ALT` | alt.mtx from vartrix (Stage 5) |
| `--clusters <path>` | `SOUPX_CLUSTERS` | clusters.tsv from troublet (Stage 7) |
| `--ambient <path>` | `SOUPX_AMBIENT` | ambient_rna.txt from consensus.py (Stage 8) |

### Locus selection (VCF modes)

| Flag | Env key | Default | Description |
|---|---|---|---|
| `--souporcell_vcf <path>` | `SOUPX_SOUPORCELL_VCF` | — | cluster_genotypes.vcf — enables Mode B positional matching |
| `--freebayes_vcf <path>` | `SOUPX_FREEBAYES_VCF` | — | freebayes VCF — combined with above enables Mode A exact CHROM:POS |
| `--min_ref <int>` | `SOUPX_MIN_REF` | `4` | Mode C fallback: min cells with ≥1 ref read |
| `--min_alt <int>` | `SOUPX_MIN_ALT` | `4` | Mode C fallback: min cells with ≥1 alt read |

### Output and config

| Flag | Env key | Default | Description |
|---|---|---|---|
| `--output <path>` | `SOUPX_OUTPUT` | `./soupxcell_output` | Output directory |
| `--alpha <float>` | `SOUPX_ALPHA` | from file | Override α from ambient_rna.txt |
| `--config <path>` | — | — | JSON parameter profile |
| `--env <path>` | — | auto | Explicit .soupxcell.env path |

### Module 1 — Correction

| Flag | Env key | Default | Description |
|---|---|---|---|
| `--threads <int>` | `SOUPX_THREADS` | `1` | Parallel threads (Rayon) |
| `--seed <u64>` | `SOUPX_SEED` | `42` | RNG seed for reproducibility |
| `--singlets_only` | `SOUPX_SINGLETS_ONLY` | `true` | Use only singlet cells for μ estimation |

### Module 3 — Embedding

| Flag | Env key | Default | Description |
|---|---|---|---|
| `--embed <str>` | `SOUPX_EMBED` | `tsne` | `tsne` \| `umap` \| `both` \| `none` |
| `--pca_components <int>` | `SOUPX_PCA_COMPONENTS` | `50` | PCA dims before t-SNE/UMAP |
| `--tsne_perplexity <float>` | `SOUPX_TSNE_PERPLEXITY` | `30.0` | t-SNE neighbourhood size |
| `--tsne_iter <int>` | `SOUPX_TSNE_ITER` | `1000` | t-SNE iterations |

### Module 2 — Simulation

| Flag | Env key | Default | Description |
|---|---|---|---|
| `--simulate` | `SOUPX_SIMULATE` | off | Enable simulation benchmark |
| `--sim_levels <floats>` | `SOUPX_SIM_LEVELS` | `0.01 0.05 0.10 0.20 0.30 0.50` | Space-separated α levels |
| `--sim_trials <int>` | `SOUPX_SIM_TRIALS` | `3` | Trials per level |

### Metrics tuning

| Flag | Env key | Default | Description |
|---|---|---|---|
| `--metric_loci <int>` | `SOUPX_METRIC_LOCI` | `0` | Max loci for silhouette (0=all) |
| `--metric_cells <int>` | `SOUPX_METRIC_CELLS` | `150` | Max cells/cluster for silhouette |

### Pre-flight and GUI

| Flag | Env key | Default | Description |
|---|---|---|---|
| `--dry_run` | `SOUPX_DRY_RUN` | off | Print plan, prompt [Y/n] |
| `--dry_run_yes` | `SOUPX_DRY_RUN_YES` | off | Auto-approve (requires `--dry_run`) |
| `--gui` | — | — | Launch web GUI at http://localhost:7878 |
| `--gui_port <int>` | — | `7878` | GUI server port |

---

## 9. Module 1 — Ambient Subtraction

**Source:** `src/core/corrector.rs`

```
Step 1:  X[l,n] = A[l,n] / (R[l,n] + A[l,n])     NaN where coverage=0
Step 2:  μ[l]   = mean of X[l,n] over singlet cells
         soup[l] = α × μ[l]
Step 3:  X_corrected[l,n] = max(0, X[l,n] − soup[l])
```

**Test dataset results:**
- 80,712 loci × 3,639 cells → 0.49% covered
- α = 0.01261973, entries floored = 640,425 (44.3%), mean correction = 0.00484
- Silhouette Δ ≈ 0 at α=1.26% — expected; within-std Δ = −0.0025 is the real signal

**Outputs:** `X_corrected.mtx`, `X_raw.mtx`, `soup_vector.csv`, `correction_summary.csv`

---

## 10. Module 2 — Simulation Benchmark

**Source:** `src/core/simulator.rs`

Supervisor: *"That is a pretty small amount of ambient RNA. We might need to simulate more to see a difference."*

For each `α_sim` level: inject `α_sim × μ[l]` → recover → measure RMSE and silhouette Δ.

**Test dataset results:** RMSE = 0.3668 × α_sim (R² = 1.000 — perfect linear recovery).
Silhouette effect visible at α ≥ 0.10.

**Outputs:** `simulation_benchmark.csv`, `figures/simulation_curves.svg`, `figures/rmse_linearity.svg`

---

## 11. Module 3 — Dimensionality Reduction

**Source:** `src/core/embedder.rs`

```
X_raw and X_corrected [L_qc × N_singlets]
  → PCA (randomised SVD): → [N × 50]         ~0.2s
  → t-SNE (Barnes-Hut):  → [N × 2]           ~42s for N=3522, 1000 iters
  → UMAP (k-NN force):   → [N × 2]
  → Cluster metrics: silhouette, Davies-Bouldin, mean within-cluster std
     — computed before AND after correction
```

**Outputs:** `figures/tsne_before_after.svg`, `figures/umap_before_after.svg`,
`figures/silhouette_comparison.svg`, `figures/correction_distribution.svg`

---

## 12. Locus selection — VCF modes A / B / C

```
Mode A  --souporcell_vcf + --freebayes_vcf  →  exact CHROM:POS cross-reference
Mode B  --souporcell_vcf alone              →  positional (VCF row N = matrix row N)
Mode C  --min_ref / --min_alt fallback      →  reproduce souporcell QC filter
```

`cluster_genotypes.vcf` from Stage 8 is the definitive locus list — use Mode B
or Mode A whenever possible.

```bash
# Mode B (recommended)
$BIN --souporcell_vcf $VCF ...

# Mode A (maximum reproducibility)
$BIN --souporcell_vcf $VCF \
     --freebayes_vcf /Volumes/extra_space/demux_test/souporcell_merged_sorted_vcf.vcf.gz ...

# Mode C (no VCF available)
$BIN --min_ref 4 --min_alt 4 ...
```

---

## 13. Step-by-step run commands

### Variables (set once)

```bash
BIN=./target/release/soupxcell
REF=/Volumes/extra_space/demux_test/ref.mtx
ALT=/Volumes/extra_space/demux_test/alt.mtx
CLU=/Volumes/extra_space/demux_test/clusters.tsv
AMB=/Volumes/extra_space/demux_test/ambient_rna.txt
VCF=/Volumes/extra_space/demux_test/cluster_genotypes.vcf
OUT=~/Setups/ambientRNA_research/souporcell_pipeline_components/playground/output_soupXCell
mkdir -p $OUT

# Verify all inputs exist
for f in $REF $ALT $CLU $AMB $VCF; do
  [ -f "$f" ] && echo "OK  $f" || echo "MISSING  $f"
done
```

### Step 0 — Dry run

```bash
$BIN --dry_run \
  --ref $REF --alt $ALT --clusters $CLU --ambient $AMB \
  --output $OUT --souporcell_vcf $VCF \
  --embed tsne --simulate --threads 8
```

### Step 1 — Module 1 only (~2s)

```bash
$BIN \
  --ref $REF --alt $ALT --clusters $CLU --ambient $AMB \
  --output $OUT --souporcell_vcf $VCF \
  --embed none --threads 8 --singlets_only
```

### Step 2 — Standard run: Module 1 + t-SNE (~100s)

```bash
$BIN \
  --ref $REF --alt $ALT --clusters $CLU --ambient $AMB \
  --output $OUT --souporcell_vcf $VCF \
  --embed tsne --pca_components 50 --tsne_perplexity 30 --tsne_iter 1000 \
  --threads 8 --seed 42 --singlets_only
```

### Step 3 — Simulation benchmark (~300s)

```bash
$BIN \
  --ref $REF --alt $ALT --clusters $CLU --ambient $AMB \
  --output $OUT --souporcell_vcf $VCF \
  --simulate --sim_levels 0.01 0.05 0.10 0.20 0.30 0.50 --sim_trials 3 \
  --embed none --threads 8
```

### Step 4 — Full experiment (~10 min)

```bash
$BIN \
  --ref $REF --alt $ALT --clusters $CLU --ambient $AMB \
  --output $OUT --souporcell_vcf $VCF \
  --simulate --sim_levels 0.01 0.05 0.10 0.20 0.30 0.50 --sim_trials 5 \
  --embed both --pca_components 50 --tsne_perplexity 30 --tsne_iter 1000 \
  --threads 8 --seed 42 --singlets_only
```

### Step 5 — Alpha override

```bash
$BIN \
  --ref $REF --alt $ALT --clusters $CLU --ambient $AMB \
  --output ${OUT}_alpha5pct --souporcell_vcf $VCF \
  --alpha 0.05 --embed tsne --threads 8 --seed 42
```

### Step 6 — Post-run checks

```bash
cat $OUT/metrics.txt
grep 'silhouette_delta' $OUT/metrics.txt
column -t -s, $OUT/simulation_benchmark.csv
open $OUT/soupxcell_report.html
```

### Step 7 — Web GUI

```bash
$BIN --gui    # open http://localhost:7878
```

---

## 14. All 9 output figures

All written to `<output>/figures/`. Self-contained SVG files.

| # | Filename | Canvas | When generated | What it shows |
|---|---|---|---|---|
| 1 | `tsne_before_after.svg` | 1200×580 | `--embed tsne/both` | Side-by-side t-SNE scatter before/after correction. Cells coloured by cluster (16-colour palette). Silhouette score on each panel. Δ annotation at bottom. |
| 2 | `umap_before_after.svg` | 1200×580 | `--embed umap/both` | Identical layout to #1 using UMAP 2D coordinates. |
| 3 | `silhouette_comparison.svg` | 800×480 | Always | Grouped bar chart: per-cluster silhouette before (semi-transparent) and after (solid). Dashed global reference lines. Per-cluster Δ labels in green/red. |
| 4 | `correction_distribution.svg` | 720×400 | Always | 40-bin histogram of per-cell mean correction magnitude. Stats box: n_cells, mean correction, floor%, α. |
| 5 | `soup_vector_profile.svg` | 900×700 | Always | Top 30 loci by `soup[l]=α×μ[l]` (sorted desc, coloured cyan→violet by intensity) + histogram of all L loci's soup values with median line. |
| 6 | `per_cluster_correction.svg` | 900×h | Always | One row per cluster with 4 metric bars: cell count, mean correction (violet), correction proxy (amber), silhouette Δ (green/red dot). |
| 7 | `floor_per_locus.svg` | 900×h | Always | Top-25 loci by floor count (red >60%, amber >30%, violet <30%) + distribution of floor counts across all QC loci. |
| 8 | `simulation_curves.svg` | 900×520 | `--simulate` | RMSE (cyan, left axis) and silhouette before/after (amber/green, right axis) vs α_sim. True α marked with dashed amber line. |
| 9 | `rmse_linearity.svg` | 800×580 | `--simulate` | Scatter of (α_sim, RMSE) + OLS fit line + R² + slope annotation + residual panel. Confirms RMSE = 0.367×α (R²=1.000). |

---

## 15. Complete output file reference

```
<output>/
├── X_corrected.mtx           corrected allele fraction matrix [L_raw × N]  ~15 MB
├── X_raw.mtx                 original allele fraction matrix  [L_raw × N]  ~15 MB
├── soup_vector.csv           locus_idx, locus_mean, soup_contribution, α    ~3 MB
├── correction_summary.csv    barcode, cluster, status, mean_correction      ~151 KB
├── metrics.txt               α, loci, cells, floor%, Δsilhouette, ΔDB      ~470 B
├── simulation_benchmark.csv  α_sim, rmse, sil_before, sil_after, floor%    ~0.5 KB
├── soupxcell_report.html     self-contained HTML, 10 sections               ~530 KB
└── figures/
    ├── tsne_before_after.svg          plot 1  (if --embed tsne/both)
    ├── umap_before_after.svg          plot 2  (if --embed umap/both)
    ├── silhouette_comparison.svg      plot 3  always
    ├── correction_distribution.svg    plot 4  always
    ├── soup_vector_profile.svg        plot 5  always
    ├── per_cluster_correction.svg     plot 6  always
    ├── floor_per_locus.svg            plot 7  always
    ├── simulation_curves.svg          plot 8  if --simulate
    └── rmse_linearity.svg             plot 9  if --simulate
```

---

## 16. Config profiles (.json)

```json
{
  "profile": "my_profile",
  "description": "Human-readable description",
  "params": {
    "singlets_only": true,
    "embed": "tsne",
    "pca_components": 50,
    "tsne_perplexity": 30.0,
    "tsne_iter": 1000,
    "seed": 42,
    "simulate": false,
    "sim_levels": "0.01,0.05,0.10,0.20,0.30,0.50",
    "sim_trials": 3,
    "threads": 8
  }
}
```

| Profile | Purpose | Key settings |
|---|---|---|
| `profiles/default.json` | Standard run | t-SNE, no simulation |
| `profiles/fast_test.json` | Smoke test | PCA=10, perplexity=10, iter=100 |
| `profiles/simulate_high.json` | Full benchmark | Both embeds, α=[0.01…1.0], 5 trials |

---

## 17. Environment file (.soupxcell.env)

Auto-discovered in CWD then `~/.soupxcell.env`. All keys must start with `SOUPX_`.

```ini
SOUPX_REF=/Volumes/extra_space/demux_test/ref.mtx
SOUPX_ALT=/Volumes/extra_space/demux_test/alt.mtx
SOUPX_CLUSTERS=/Volumes/extra_space/demux_test/clusters.tsv
SOUPX_AMBIENT=/Volumes/extra_space/demux_test/ambient_rna.txt
SOUPX_OUTPUT=~/Setups/ambientRNA_research/souporcell_pipeline_components/playground/output_soupXCell
SOUPX_THREADS=8
SOUPX_SEED=42
SOUPX_SINGLETS_ONLY=true
SOUPX_EMBED=tsne
SOUPX_PCA_COMPONENTS=50
SOUPX_TSNE_PERPLEXITY=30.0
SOUPX_TSNE_ITER=1000
SOUPX_SIMULATE=false
SOUPX_SIM_LEVELS=0.01,0.05,0.10,0.20,0.30,0.50
SOUPX_SIM_TRIALS=3
# SOUPX_SOUPORCELL_VCF=/Volumes/extra_space/demux_test/cluster_genotypes.vcf
# SOUPX_ALPHA=0.05
```

---

## 18. Reading the logs

| Prefix | Meaning |
|---|---|
| `[SOUPX] SECTION` | Major pipeline stage |
| `[SOUPX] METRIC` | Numeric result: `key = value  unit` |
| `[SOUPX] OK` | Stage complete with elapsed time |
| `[SOUPX] WRITTEN` | File written |
| `[SOUPX] WARN / ERROR` | Non-fatal / fatal issue |
| `SOUPX_SIM\t...` | Tab-delimited simulation row |

```bash
# Capture all output
$BIN ... 2> run.log

# Extract metrics
grep '\[SOUPX\] METRIC' run.log

# Simulation rows as TSV
grep '^SOUPX_SIM' run.log

# Timing
grep '\[SOUPX\] OK' run.log
```

---

## 19. Web GUI

```bash
./target/release/soupxcell --gui            # http://localhost:7878
```

Four phases: **Configure** (4 presets, all params, live CLI preview, live Mode A/B/C indicator)
→ **Pre-flight** (full plan preview) → **Execute** (live log, progress, simulation table)
→ **Outputs** (8 metric cards, interpretation panel, file download links).

---

## 20. Troubleshooting

| Problem | Fix |
|---|---|
| `cargo build` linker error (macOS) | `xcode-select --install` then retry |
| `rustc --version` < 1.75 | `rustup update stable` |
| `File not found: ref.mtx` | Check volume: `ls /Volumes/extra_space/demux_test/` |
| `0 loci from VCF` | Check path: `wc -l $VCF` should print >8000 |
| `n_cells_singlet = 0` | clusters.tsv must have `status` column from troublet |
| `alpha = 0.0000` | `cat $AMB` should print a float like `0.01261973` |
| Silhouette Δ ≈ 0 | Expected at α=1.26%. Use `--simulate` to demonstrate larger effects |
| Floor% > 50% | α may be overestimated — verify `$AMB` matches this souporcell run |
| t-SNE looks noisy | `--tsne_iter 3000 --tsne_perplexity 50` or try `--embed umap` |
| Out of memory | `--embed none` first (saves ~300 MB). Peak ~500 MB with t-SNE |
| Report figures missing | Open `soupxcell_report.html` from the output directory — SVGs use relative paths |

---

## Expected results (test dataset)

| Metric | Value |
|---|---|
| Matrix | 80,712 loci × 3,639 cells (0.49% covered) |
| Loci after QC (VCF Mode B) | 8,219 |
| α | 0.01261973 (1.26%) |
| Entries floored | 640,425 (44.3% of covered) |
| Mean \|correction\| | 0.00484 |
| Silhouette Δ | ≈ 0 (expected at 1.26%) |
| Within-std Δ | −0.0025 (primary signal) |
| RMSE linearity | RMSE = 0.3668 × α (R² = 1.000) |
| Standard run time | ~100s (t-SNE, 1000 iters) |
| Full experiment | ~10 min (both embeds + simulation) |

---

*soupXcell v1.0 · Jahidul Arafat · souporcell: Heaton et al., Nature Methods 2020*
