// ============================================================================
// infra/preflight.rs — --dry_run plan printer and [Y/n] approval gate
// ============================================================================

use std::io::{self, BufRead, Write};
use crate::config::params::Params;
use crate::infra::logger::Logger;

pub enum ApprovalResult { Approved, Rejected }

pub fn show_preflight(params: &Params, _log: &Logger) -> ApprovalResult {
    let sep = "─".repeat(78);

    eprintln!("\n\x1b[1m\x1b[36m╔══════════════════════════════════════════════════════════════════╗");
    eprintln!("║   soupxcell v1.0  ·  DRY-RUN PLAN                               ║");
    eprintln!("╚══════════════════════════════════════════════════════════════════╝\x1b[0m\n");

    // Sources
    eprintln!("\x1b[33m§ Parameter Sources\x1b[0m");
    eprintln!("{}", sep);
    if let Some(s) = params.param_sources.first() {
        eprintln!("  {}", s.1);
    }

    // Inputs
    eprintln!("\n\x1b[33m§ Input Files\x1b[0m");
    eprintln!("{}", sep);
    eprintln!("  ref matrix         : {}", params.ref_matrix);
    eprintln!("  alt matrix         : {}", params.alt_matrix);
    eprintln!("  clusters.tsv       : {}", params.clusters);
    eprintln!("  ambient_rna.txt    : {}", params.ambient);
    eprintln!("  α from file        : will be parsed at runtime");
    if let Some(a) = params.alpha_override {
        eprintln!("  α override         : {:.8}  (CLI override — file value ignored)", a);
    }

    // Core settings
    eprintln!("\n\x1b[33m§ Correction Settings\x1b[0m");
    eprintln!("{}", sep);
    eprintln!("  singlets only      : {}", params.singlets_only);
    eprintln!("  threads            : {}", params.threads);
    eprintln!("  seed               : {}", params.seed);

    // Module 1
    eprintln!("\n\x1b[33m§ Module 1 — Ambient Subtraction\x1b[0m");
    eprintln!("{}", sep);
    eprintln!("  STEP 1  Compute allele fraction matrix");
    eprintln!("          X[l,n] = alt[l,n] / (ref[l,n] + alt[l,n])");
    eprintln!("          NaN where coverage = 0");
    eprintln!("  STEP 2  Compute per-locus soup contribution");
    eprintln!("          μ[l] = mean(X[l,:]) over all cells with coverage");
    eprintln!("          soup[l] = α × μ[l]");
    eprintln!("  STEP 3  Subtract and floor at 0");
    eprintln!("          X_corrected[l,n] = max(0,  X[l,n] − soup[l])");

    // Embedding
    eprintln!("\n\x1b[33m§ Module 3 — Embedding & Visualisation\x1b[0m");
    eprintln!("{}", sep);
    eprintln!("  embed              : {:?}", params.embed);
    eprintln!("  PCA components     : {}", params.pca_components);
    eprintln!("  t-SNE perplexity   : {}", params.tsne_perplexity);
    eprintln!("  t-SNE iterations   : {}", params.tsne_iter);

    // Module 2
    if params.simulate {
        eprintln!("\n\x1b[33m§ Module 2 — Simulation Benchmark\x1b[0m");
        eprintln!("{}", sep);
        eprintln!("  α simulation levels: {:?}", params.sim_levels);
        eprintln!("  trials per level   : {}", params.sim_trials);
        eprintln!("  For each α_sim:");
        eprintln!("    X_contaminated = X_raw + α_sim × μ  (inject)");
        eprintln!("    X_recovered    = max(0, X_contaminated − α_sim × μ)  (recover)");
        eprintln!("    Metrics: RMSE, silhouette before/after, floor fraction");
    }

    // Outputs
    eprintln!("\n\x1b[33m§ Expected Outputs\x1b[0m");
    eprintln!("{}", sep);
    eprintln!("  {}/", params.output);
    eprintln!("  ├── X_corrected.mtx           corrected allele fraction matrix");
    eprintln!("  ├── X_raw.mtx                 original allele fraction matrix");
    eprintln!("  ├── soup_vector.csv            per-locus soup contribution");
    eprintln!("  ├── correction_summary.csv     per-cell correction stats");
    eprintln!("  ├── metrics.txt                key quality metrics");
    eprintln!("  ├── soupxcell_report.html      self-contained HTML report");
    eprintln!("  └── figures/");
    if params.embed.run_tsne() {
        eprintln!("      ├── tsne_before_after.svg");
    }
    if params.embed.run_umap() {
        eprintln!("      ├── umap_before_after.svg");
    }
    eprintln!("      ├── silhouette_comparison.svg");
    eprintln!("      └── correction_distribution.svg");
    if params.simulate {
        eprintln!("      └── simulation_curves.svg");
    }

    // Prompt
    if params.dry_run_yes {
        eprintln!("\n\x1b[32m  [dry_run_yes] Auto-approved — proceeding.\x1b[0m\n");
        return ApprovalResult::Approved;
    }

    eprint!("\n  Proceed with this run? [Y/n]: ");
    io::stderr().flush().ok();
    io::stdout().flush().ok();

    let stdin = io::stdin();
    let line = stdin.lock().lines().next()
        .and_then(|l| l.ok())
        .unwrap_or_default();
    let answer = line.trim().to_lowercase();

    if answer == "n" || answer == "no" {
        eprintln!("\n\x1b[33m  Aborted — nothing was read or written.\x1b[0m\n");
        ApprovalResult::Rejected
    } else {
        eprintln!("\x1b[32m  Approved — starting run.\x1b[0m\n");
        ApprovalResult::Approved
    }
}
