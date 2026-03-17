# soupXcell v1.0 — Ambient RNA Correction for Genotype Demultiplexing

<p align="center">
  <strong>Post-processor for the souporcell pipeline — correct, visualise, and benchmark ambient RNA contamination in pooled single-cell RNA-seq</strong>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/version-1.0.0-8B5CF6" alt="version"/>
  <img src="https://img.shields.io/badge/rust-≥1.75-orange" alt="rust"/>
  <img src="https://img.shields.io/badge/edition-2021-blue" alt="edition"/>
  <img src="https://img.shields.io/badge/figures-9%20HTML%20+%20JSON-10B981" alt="figures"/>
  <img src="https://img.shields.io/badge/report-12%20sections-6366F1" alt="report"/>
</p>

---

> **souporcell pipeline:** Heaton et al., *Nature Methods* 2020 · [doi:10.1038/s41592-020-0820-1](https://doi.org/10.1038/s41592-020-0820-1)
>
> **soupXcell:** Jahidul Arafat — ambient RNA correction, simulation benchmark, interactive HTML figures, storytelling HTML report

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
14. [Output figures — HTML + JSON architecture](#14-output-figures--html--json-architecture)
15. [Complete output file reference](#15-complete-output-file-reference)
16. [The HTML report — 12 sections](#16-the-html-report--12-sections)
17. [Config profiles (.json)](#17-config-profiles-json)
18. [Environment file (.soupxcell.env)](#18-environment-file-soupxcellenv)
19. [Reading the logs](#19-reading-the-logs)
20. [Web GUI](#20-web-gui)
21. [Troubleshooting](#21-troubleshooting)

---

## 1. What soupXcell does

soupXcell is a **post-processor** for the souporcell pipeline. It takes the
allele count matrices and cluster assignments produced by the pipeline and:

1. **Corrects** each cell's allele fraction measurements by subtracting the
   estimated ambient RNA contribution (Module 1)
2. **Benchmarks** the correction at multiple contamination levels using synthetic
   injection — directly addressing the supervisor's concern that 1.26% is too
   small to see a visible effect (Module 2)
3. **Visualises** the effect via t-SNE and/or UMAP embeddings with **9 interactive
   HTML figures** backed by **JSON data files** and a complete **12-section HTML
   report** with storytelling narrative (Module 3)

Runs as a **CLI tool** (with `--dry_run` pre-flight), as a **web GUI** (`--gui`),
and produces a self-contained HTML report, 9 interactive HTML figures, and all
underlying data as JSON files in `figures/image_data/`.

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
                              9 interactive HTML figures
                              figures/image_data/*.json (9 data files)
                              figures/cluster_data.json
                              soupxcell_report.html (12-section report)
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
    ├── main.rs                870L   — composition root; writes JSON data files
    │                                   + copies HTML templates to output figures/
    ├── domain/
    │   ├── types.rs           200L   — all structs (CorrectionResult, RunSummary,
    │   │                               BenchmarkResult, ClusterMetrics, CellMeta…)
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
    │   ├── soupxcell_gui.html 1321L  — web GUI (4-phase mission control)
    │   └── figures_html/             — HTML figure templates (11 files)
    │       ├── tsne_before_after.html
    │       ├── umap_before_after.html
    │       ├── silhouette_comparison.html
    │       ├── correction_distribution.html
    │       ├── simulation_curves.html
    │       ├── soup_vector_profile.html
    │       ├── per_cluster_correction.html
    │       ├── floor_per_locus.html
    │       ├── rmse_linearity.html
    │       ├── cluster_composition.html
    │       └── contamination_progression.html
    └── analysis/
        └── report/            ~900L  — 11-module HTML report system
            ├── mod.rs                  entry point, ReportData struct,
            │                           4 new builder fns (research, findings,
            │                           story, sim_bars)
            ├── css.rs                  full-width professional CSS
            │                           (max-width:1240px, scroll-reveal,
            │                            finding-boxes, timeline, sim-bars)
            ├── utils.rs                colour constants, chrono_now()
            ├── cards.rs                12 metric cards
            ├── experiment.rs           §3 experiment design section
            ├── tables.rs               §5 results tables
            ├── svgs.rs                 §6 inline SVG diagnostics
            ├── plot_embeds.rs          §7 HTML iframe embeds
            ├── interpretation.rs       §8 per-metric commentary
            ├── params.rs               §9 27-parameter table
            └── methods.rs              §10 methods + references
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

# IMPORTANT: always export $VCF before running so the variable is set
export VCF
echo $VCF    # must print a path, not a blank line

# Smoke test (~15s)
$BIN --ref $REF --alt $ALT --clusters $CLU --ambient $AMB \
     --output $OUT/smoke \
     --souporcell_vcf $VCF \
     --embed tsne --pca_components 10 --tsne_iter 100 --threads 4

# Standard run (~25 min with simulation)
$BIN --ref $REF --alt $ALT --clusters $CLU --ambient $AMB \
     --output $OUT --souporcell_vcf $VCF \
     --singlets_only \
     --simulate --sim_levels 0.01 0.05 0.10 0.20 0.30 0.50 --sim_trials 5 \
     --embed both --pca_components 50 --tsne_perplexity 30 --tsne_iter 1000 \
     --threads 8 --seed 42

# Open report and figures
open $OUT/soupxcell_report.html
open $OUT/figures/tsne_before_after.html
open $OUT/figures/contamination_progression.html
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
| `--souporcell_vcf <path>` | `SOUPX_SOUPORCELL_VCF` | — | cluster_genotypes.vcf — enables Mode B positional matching. **Always set `export VCF` before running.** |
| `--freebayes_vcf <path>` | `SOUPX_FREEBAYES_VCF` | — | freebayes VCF — combined with above enables Mode A exact CHROM:POS |
| `--min_ref <int>` | `SOUPX_MIN_REF` | `4` | Mode C fallback: min cells with ≥1 ref read |
| `--min_alt <int>` | `SOUPX_MIN_ALT` | `4` | Mode C fallback: min cells with ≥1 alt read |

> **Common mistake:** Running with `--souporcell_vcf $VCF` when `$VCF` is not
> exported causes the binary to receive an empty string, falling back to Mode C.
> The log will show `souporcell_vcf resolved: <path>` so you can verify.

### Output and config

| Flag | Env key | Default | Description |
|---|---|---|---|
| `--output <path>` | `SOUPX_OUTPUT` | `./soupxcell_output` | Output directory (created if missing) |
| `--alpha <float>` | `SOUPX_ALPHA` | from file | Override α from ambient_rna.txt |
| `--config <path>` | — | — | JSON parameter profile |
| `--env <path>` | — | auto | Explicit .soupxcell.env path |

### Module 1 — Correction

| Flag | Env key | Default | Description |
|---|---|---|---|
| `--threads <int>` | `SOUPX_THREADS` | `1` | Parallel threads (Rayon) |
| `--seed <u64>` | `SOUPX_SEED` | `42` | RNG seed for reproducibility |
| `--singlets_only` | `SOUPX_SINGLETS_ONLY` | `true` | Use only singlet cells for μ estimation |

### Module 2 — Simulation

| Flag | Env key | Default | Description |
|---|---|---|---|
| `--simulate` | `SOUPX_SIMULATE` | off | Enable Module 2 simulation benchmark |
| `--sim_levels <floats>` | `SOUPX_SIM_LEVELS` | `0.01 0.05 0.10 0.20 0.30 0.50` | Space-separated α levels to test |
| `--sim_trials <int>` | `SOUPX_SIM_TRIALS` | `3` | Trials averaged per level |

### Module 3 — Embedding

| Flag | Env key | Default | Description |
|---|---|---|---|
| `--embed <str>` | `SOUPX_EMBED` | `tsne` | `tsne` \| `umap` \| `both` \| `none` |
| `--pca_components <int>` | `SOUPX_PCA_COMPONENTS` | `50` | PCA dims before t-SNE/UMAP |
| `--tsne_perplexity <float>` | `SOUPX_TSNE_PERPLEXITY` | `30.0` | t-SNE neighbourhood size |
| `--tsne_iter <int>` | `SOUPX_TSNE_ITER` | `1000` | t-SNE iterations |

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

**Test dataset results (Mode B, VCF-selected loci):**
- 80,712 loci × 3,639 cells → 0.49% covered
- QC loci (from VCF): 64,496
- α = 0.01261973, entries floored = 640,425 (44.3% of covered)
- Mean |correction| = 0.005201, Silhouette Δ ≈ 0 (expected at 1.26%)
- Within-cluster std Δ = −0.000887 (primary signal of correction working)

**Outputs:** `X_corrected.mtx`, `X_raw.mtx`, `soup_vector.csv`, `correction_summary.csv`

---

## 10. Module 2 — Simulation Benchmark

**Source:** `src/core/simulator.rs`

Supervisor: *"That is a pretty small amount of ambient RNA. We might need to simulate more to see a difference."*

For each `α_sim` level: inject `α_sim × μ[l]` → correct → measure RMSE and silhouette Δ.
Each level averaged over `--sim_trials` independent trials.

**Test dataset results (Mode B, 5 trials):**

| α_sim | RMSE | Sil Δ | Interpretation |
|---|---|---|---|
| 1% | 0.003873 | ≈0 | Below detection threshold |
| 5% | 0.019366 | ≈0 | Below detection threshold |
| 10% | 0.038733 | ≈0 | Below detection threshold |
| 20% | 0.077471 | ≈0 | Donor identity encoded in allele pattern, not magnitude |
| 30% | 0.116218 | ≈0 | Same — correction is still exact |
| 50% | 0.193779 | −0.0001 | Barely detectable |

**Key finding:** RMSE = 0.3875 × α_sim (R² = 1.000) — perfect linear recovery.
Silhouette Δ ≈ 0 even at 50% because donors are distinguished by *which loci*
show alt/ref, not by allele fraction magnitude. This is scientifically correct.

**Outputs:**
- `simulation_benchmark.csv`
- `figures/image_data/simulation_curves.json` ← data for HTML figures
- `figures/image_data/rmse_linearity.json`
- `figures/simulation_curves.html`
- `figures/rmse_linearity.html`
- `figures/contamination_progression.html` ← interactive level selector

---

## 11. Module 3 — Dimensionality Reduction

**Source:** `src/core/embedder.rs`

```
X_raw and X_corrected [L_qc × N_singlets]
  → PCA (randomised SVD): → [N × 50]         ~20s for N=3522
  → t-SNE (Barnes-Hut):  → [N × 2]           ~47s for N=3522, 1000 iters
  → UMAP (k-NN force):   → [N × 2]           ~0.5s (simplified layout)
  → Cluster metrics: silhouette, Davies-Bouldin, mean within-cluster std
     — computed before AND after correction on QC loci
```

Embedding coordinates are written as JSON to `figures/image_data/tsne_before_after.json`
and `figures/image_data/umap_before_after.json` containing all cell (x, y) coordinates
and cluster labels. The HTML figures read these files directly.

**Outputs:**
- `figures/image_data/tsne_before_after.json` (240 KB for 3,522 cells)
- `figures/image_data/umap_before_after.json` (225 KB)
- `figures/tsne_before_after.html`
- `figures/umap_before_after.html`

---

## 12. Locus selection — VCF modes A / B / C

```
Mode A  --souporcell_vcf + --freebayes_vcf  →  exact CHROM:POS cross-reference
Mode B  --souporcell_vcf alone              →  positional (VCF row N = matrix row N)
Mode C  --min_ref / --min_alt fallback      →  reproduce souporcell QC filter
```

The log always prints `souporcell_vcf resolved: <path>` so you can verify
the path was received. A blank line here means `$VCF` was not exported.

`cluster_genotypes.vcf` from Stage 8 is the definitive locus list — use Mode B
or Mode A whenever possible. Mode B gives 64,496 loci on the test dataset vs
8,219 for Mode C.

```bash
# Mode B (recommended)
export VCF=/Volumes/extra_space/demux_test/cluster_genotypes.vcf
$BIN --souporcell_vcf $VCF ...

# Mode A (maximum reproducibility)
$BIN --souporcell_vcf $VCF \
     --freebayes_vcf /Volumes/extra_space/demux_test/souporcell_merged_sorted_vcf.vcf.gz ...

# Mode C (no VCF available)
$BIN --min_ref 4 --min_alt 4 ...
```

---

## 13. Step-by-step run commands

### Variables (set once, always export)

```bash
export BIN=./target/release/soupxcell
export REF=/Volumes/extra_space/demux_test/ref.mtx
export ALT=/Volumes/extra_space/demux_test/alt.mtx
export CLU=/Volumes/extra_space/demux_test/clusters.tsv
export AMB=/Volumes/extra_space/demux_test/ambient_rna.txt
export VCF=/Volumes/extra_space/demux_test/cluster_genotypes.vcf
export OUT=~/Setups/ambientRNA_research/souporcell_pipeline_components/playground/output_soupXCell
mkdir -p $OUT

# Verify all inputs and that $VCF is set
for f in $REF $ALT $CLU $AMB $VCF; do
  [ -f "$f" ] && echo "OK  $f" || echo "MISSING  $f"
done
echo "VCF variable: $VCF"    # must not be blank
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

### Step 2 — Module 1 + t-SNE (~2 min)

```bash
$BIN \
  --ref $REF --alt $ALT --clusters $CLU --ambient $AMB \
  --output $OUT --souporcell_vcf $VCF \
  --embed tsne --pca_components 50 --tsne_perplexity 30 --tsne_iter 1000 \
  --threads 8 --seed 42 --singlets_only
```

### Step 3 — Simulation benchmark only (~22 min for 6 levels × 5 trials)

```bash
$BIN \
  --ref $REF --alt $ALT --clusters $CLU --ambient $AMB \
  --output $OUT --souporcell_vcf $VCF \
  --simulate --sim_levels 0.01 0.05 0.10 0.20 0.30 0.50 --sim_trials 5 \
  --embed none --threads 8
```

### Step 4 — Full experiment (~25 min)

```bash
$BIN \
  --ref $REF --alt $ALT --clusters $CLU --ambient $AMB \
  --output $OUT --souporcell_vcf $VCF \
  --singlets_only \
  --simulate --sim_levels 0.01 0.05 0.10 0.20 0.30 0.50 --sim_trials 5 \
  --embed both --pca_components 50 --tsne_perplexity 30 --tsne_iter 1000 \
  --threads 8 --seed 42
```

### Step 5 — Post-run checks

```bash
cat $OUT/metrics.txt
grep 'silhouette_delta' $OUT/metrics.txt
column -t -s, $OUT/simulation_benchmark.csv

open $OUT/soupxcell_report.html
open $OUT/figures/tsne_before_after.html
open $OUT/figures/contamination_progression.html
open $OUT/figures/cluster_composition.html
```

---

## 14. Output figures — HTML + JSON architecture

All figures follow a **two-file pattern**:

```
figures/
  image_data/<name>.json      ← exact data written by Rust at runtime
  <name>.html                 ← interactive figure that reads the JSON
```

The HTML files are templates copied from `src/infra/figures_html/` to the output
directory at runtime. They `fetch()` their corresponding JSON file and render
using Canvas 2D. **No data is hardcoded** — every number comes from the JSON.

If a JSON file is missing (e.g., `--simulate` not run), the HTML displays a
clear error message rather than showing fake data.

### All 11 figures

| # | HTML file | JSON data file | When written | What it shows |
|---|---|---|---|---|
| 1 | `tsne_before_after.html` | `image_data/tsne_before_after.json` | `--embed tsne/both` | Side-by-side scatter of all cells in t-SNE space, before and after correction. Coloured by donor cluster. Hover any cell for coordinates, cluster, silhouette. |
| 2 | `umap_before_after.html` | `image_data/umap_before_after.json` | `--embed umap/both` | Same layout as #1 using UMAP coordinates. |
| 3 | `silhouette_comparison.html` | `image_data/silhouette_comparison.json` | Always | Grouped bars: per-cluster silhouette before (faded) and after (solid). Dashed global reference lines. Δ per cluster with colour coding. |
| 4 | `correction_distribution.html` | `image_data/correction_distribution.json` | Always | Histogram of per-cell mean correction magnitude. Mean and median lines. Hover for bin range and count. |
| 5 | `soup_vector_profile.html` | `image_data/soup_vector_profile.json` | Always | Top-30 loci by μ[l] with horizontal bars + log-scale histogram of all loci's soup values. |
| 6 | `per_cluster_correction.html` | `image_data/per_cluster_correction.json` | Always | Table-style row per cluster: cell count, mean correction, floor%, silhouette Δ. |
| 7 | `floor_per_locus.html` | `image_data/floor_per_locus.json` | Always | Top-25 loci by floor count + log-scale histogram of all floor counts. |
| 8 | `simulation_curves.html` | `image_data/simulation_curves.json` | `--simulate` | RMSE (left axis) and silhouette before/after (right axis) vs α_sim. OLS fit line, R², true α marker. |
| 9 | `rmse_linearity.html` | `image_data/rmse_linearity.json` | `--simulate` | Scatter (α_sim, RMSE) + OLS line + residuals panel. |
| 10 | `cluster_composition.html` | `cluster_data.json` | Always | Real per-cluster breakdown from clusters.tsv: singlet/doublet counts, silhouette, correction. |
| 11 | `contamination_progression.html` | `image_data/simulation_curves.json` | `--simulate` | **Interactive level selector** — click any α_sim level button to animate the scatter to that contamination state. PhD impact table below. |

All figures have **PhD-level tooltips**: hover any data point for exact values
with 6–8 decimal places, formula derivations, and biological interpretation.
All figures have **⬇ Download PNG** buttons.

---

## 15. Complete output file reference

```
<output>/
├── X_corrected.mtx            corrected allele fraction matrix [L_raw × N]  ~15 MB
├── X_raw.mtx                  original allele fraction matrix  [L_raw × N]  ~15 MB
├── soup_vector.csv            locus_idx, locus_mean, soup_contribution, α    ~3 MB
├── correction_summary.csv     barcode, cluster, status, mean_correction      ~151 KB
├── metrics.txt                α, loci, cells, floor%, Δsilhouette, ΔDB, …   ~470 B
├── simulation_benchmark.csv   α_sim, rmse, sil_before, sil_after, floor%    ~0.5 KB
├── soupxcell_report.html      12-section HTML report, full-width layout      ~65 KB
└── figures/
    ├── image_data/                 ← all figure data as JSON (no estimates)
    │   ├── tsne_before_after.json      cell coordinates + cluster labels     ~240 KB
    │   ├── umap_before_after.json      UMAP coordinates + cluster labels     ~225 KB
    │   ├── silhouette_comparison.json  per-cluster sil before/after          ~0.5 KB
    │   ├── correction_distribution.json  all cell correction values          ~31 KB
    │   ├── simulation_curves.json      all sim level results                 ~0.8 KB
    │   ├── rmse_linearity.json         same as above (different view)        ~0.8 KB
    │   ├── soup_vector_profile.json    top loci + all soup values            ~711 KB
    │   ├── per_cluster_correction.json per-cluster exact metrics             ~0.6 KB
    │   └── floor_per_locus.json        top-25 loci + all floor counts        ~139 KB
    ├── cluster_data.json           exact per-cluster breakdown from clusters.tsv
    ├── tsne_before_after.html      interactive figure #1
    ├── umap_before_after.html      interactive figure #2
    ├── silhouette_comparison.html  interactive figure #3
    ├── correction_distribution.html  interactive figure #4
    ├── simulation_curves.html      interactive figure #8  (if --simulate)
    ├── soup_vector_profile.html    interactive figure #5
    ├── per_cluster_correction.html interactive figure #6
    ├── floor_per_locus.html        interactive figure #7
    ├── rmse_linearity.html         interactive figure #9  (if --simulate)
    ├── cluster_composition.html    cluster viewer (reads cluster_data.json)
    └── contamination_progression.html  interactive progression (if --simulate)
```

---

## 16. The HTML report — 12 sections

`soupxcell_report.html` is a single self-contained file (no external dependencies
except Google Fonts). Full-width layout (`max-width: 1240px`). All sections animate
in on scroll. Content is always visible — scroll animation uses `transform` only.

| Section | ID | Contents |
|---|---|---|
| Cover | — | Title, α badge, run ID, quick stats, generated timestamp |
| Table of Contents | `#toc` | Two-column TOC with anchor links to all 12 sections |
| Research Context | `#research` | Background explanation, research objectives, research questions, 3 hypotheses, timeline of what happened |
| Key Findings | `#findings` | 6 colour-coded finding boxes with plain-English interpretation |
| Metrics Dashboard | `#metrics` | 12 metric cards (value + Δ + colour-coded bar) |
| Story Narrative | `#story` | 4-Act narrative of the run: loading → correction → measurement → simulation proof |
| Experiment Design | `#experiment` | Locus source, correction narrative, singlets/doublets |
| Formula Detail | `#formula` | Step-by-step formula with plain-English notes + run-specific values table |
| Simulation Benchmark | `#simulation` | Animated RMSE bars per level + level-by-level breakdown table |
| Results Tables | `#tables` | 4 detailed tables: correction summary, cluster metrics, simulation results, locus stats |
| Interactive Figures | `#figures` | All 11 HTML figures embedded as iframes + open links |
| Scientific Interpretation | `#interpretation` | 9 interpretation cards with technical + biological commentary |
| Parameters | `#params` | Complete 27-parameter table for reproducibility |
| Methods | `#methods` | Algorithm description + Heaton et al. 2020 reference |

### Opening the report

```bash
# macOS
open $OUT/soupxcell_report.html

# The report embeds figure iframes using relative paths.
# Always open from the output directory, not by dragging to browser.
# The figures/ subdirectory must be present alongside the HTML file.
```

---

## 17. Config profiles (.json)

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

## 18. Environment file (.soupxcell.env)

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
SOUPX_EMBED=both
SOUPX_PCA_COMPONENTS=50
SOUPX_TSNE_PERPLEXITY=30.0
SOUPX_TSNE_ITER=1000
SOUPX_SIMULATE=true
SOUPX_SIM_LEVELS=0.01,0.05,0.10,0.20,0.30,0.50
SOUPX_SIM_TRIALS=5
SOUPX_SOUPORCELL_VCF=/Volumes/extra_space/demux_test/cluster_genotypes.vcf
# SOUPX_ALPHA=0.05   # uncomment to override α
```

---

## 19. Reading the logs

| Prefix | Meaning |
|---|---|
| `[SOUPX] SECTION` | Major pipeline stage |
| `[SOUPX] METRIC` | Numeric result: `key = value  unit` |
| `[SOUPX] OK` | Stage complete with elapsed time |
| `[SOUPX] WRITTEN` | File written (path + size) |
| `[SOUPX] WARN / ERROR` | Non-fatal / fatal issue |
| `SOUPX_SIM\t...` | Tab-delimited simulation row (one per α_sim level) |

```bash
# Capture all output
$BIN ... 2> run.log

# Check VCF was received correctly
grep 'souporcell_vcf resolved' run.log

# Extract all metrics
grep '\[SOUPX\] METRIC' run.log

# Simulation rows as TSV
grep '^SOUPX_SIM' run.log

# Timing summary
grep '\[SOUPX\] OK' run.log

# Check all JSON files written
grep '\[SOUPX\] WRITTEN' run.log | grep '\.json'
```

---

## 20. Web GUI

```bash
./target/release/soupxcell --gui            # http://localhost:7878
```

Four phases: **Configure** (4 presets, all params, live CLI preview, live Mode A/B/C indicator)
→ **Pre-flight** (full plan preview) → **Execute** (live log, progress, simulation table)
→ **Outputs** (8 metric cards, interpretation panel, file download links).

---

## 21. Troubleshooting

| Problem | Fix |
|---|---|
| `cargo build` linker error (macOS) | `xcode-select --install` then retry |
| `rustc --version` < 1.75 | `rustup update stable` |
| `File not found: ref.mtx` | Check volume: `ls /Volumes/extra_space/demux_test/` |
| `0 loci from VCF` | Check path: `wc -l $VCF` should print >8000 |
| Log shows `souporcell_vcf not set` | `$VCF` was empty. Run `export VCF=<path>` first |
| `n_cells_singlet = 0` | clusters.tsv must have `status` column from troublet |
| `alpha = 0.0000` | `cat $AMB` should print a float like `0.01261973` |
| Silhouette Δ ≈ 0 | Expected at α=1.26%. Use `--simulate` to demonstrate larger effects |
| Floor% > 50% | α may be overestimated — verify `$AMB` matches this souporcell run |
| t-SNE looks noisy / no cluster structure | `--tsne_iter 3000 --tsne_perplexity 50`. At α=1.26% clusters may genuinely overlap in allele-fraction space — this is biologically correct |
| HTML figures show "JSON not found" | Open figures from the `output/figures/` directory, not from `src/infra/figures_html/` |
| Report sections invisible on load | Fixed in current version — `opacity:.85` on load, not `opacity:0` |
| Out of memory | `--embed none` first (saves ~300 MB). Peak ~500 MB with t-SNE on 3,522 cells |
| `contamination_progression.html` shows old data | Rebuild and rerun — this file is now copied *after* simulation completes |

---

## Expected results (test dataset — Mode B, full run)

| Metric | Value |
|---|---|
| Matrix | 80,712 loci × 3,639 cells (0.49% covered) |
| Loci after QC (VCF Mode B) | 64,496 |
| α | 0.01261973 (1.26%) |
| Entries floored | 640,425 (44.3% of covered) |
| Mean \|correction\| | 0.005201 |
| Silhouette Δ | −0.000003 (≈0, expected at 1.26%) |
| Within-cluster std Δ | −0.000887 (primary correction signal) |
| Davies-Bouldin Δ | −0.002964 (slight improvement) |
| RMSE linearity slope | 0.3875 (RMSE = 0.3875 × α_sim, R² = 1.000) |
| t-SNE time | ~47s per embedding (1000 iters, 3522 cells) |
| Full experiment time | ~25 min (both embeds + 6 levels × 5 trials) |

---

*soupXcell v1.0 · Jahidul Arafat · souporcell: Heaton et al., Nature Methods 2020*
