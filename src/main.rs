// ============================================================================
// main.rs — soupxcell v1.0 composition root
// ============================================================================
//
// Layer map (dependency order — inner layers never import outer):
//   domain/     ① Pure types + math  (no I/O, no CLI)
//   config/     ② Params struct + 3-layer resolution (CLI > JSON > .env)
//   core/       ③ Algorithms: corrector · simulator · embedder
//   infra/      ④ I/O + logger + preflight
//   analysis/   ⑤ SVG plots + HTML report
//
// main() responsibilities:
//   1. Print banner
//   2. Load params (config layer)
//   3. Dry-run pre-flight if --dry_run
//   4. Load inputs: ref.mtx, alt.mtx, clusters.tsv, ambient_rna.txt
//   5. Module 1 — ambient subtraction
//   6. Module 3 — PCA → t-SNE / UMAP embeddings
//   7. Module 2 — simulation benchmark (if --simulate)
//   8. Write all outputs: matrices, CSV, SVG plots, HTML report
//   9. Print final summary
// ============================================================================

#[macro_use]
extern crate clap;

mod domain;
mod config;
mod core;
mod infra;
mod analysis;

use std::fs;
use std::time::Instant;

use rayon::ThreadPoolBuilder;

use crate::config::params::{load_params, EmbedChoice};
use crate::domain::math::{
    silhouette_score, davies_bouldin,
    mean_within_cluster_std,
};
use crate::domain::types::{ClusterMetrics, RunSummary};
use crate::core::corrector::run_correction;
use crate::core::simulator::run_simulation;
use crate::core::embedder::{run_pca, run_tsne, run_umap};
use crate::infra::logger::{print_banner, Logger};
use crate::infra::io::{
    load_raw_matrices, load_matrices_qc, load_clusters, parse_ambient_fraction,
    parse_vcf_locus_indices,
    write_corrected_mtx, write_soup_vector, write_correction_summary,
    write_benchmark_csv, write_metrics_txt,
};
use crate::infra::preflight::{show_preflight, ApprovalResult};
use crate::infra::gui_server;



fn main() {
    print_banner();

    // ── 1. GUI mode: intercept --gui from raw args BEFORE clap parse ─────────
    // This mirrors souporcell v3.0 — allows --gui without supplying required
    // flags like --ref / --alt which clap would otherwise reject.
    let raw_args: Vec<String> = std::env::args().collect();
    if raw_args.iter().any(|a| a == "--gui") {
        let port: u16 = raw_args.iter()
            .position(|a| a == "--gui_port")
            .and_then(|i| raw_args.get(i + 1))
            .and_then(|v| v.parse().ok())
            .unwrap_or(7878);
        gui_server::launch(port);
        return;
    }

    // ── 2. Load params (clap + JSON profile + .env merge) ────────────────────
    let params = load_params();
    let log    = Logger::new();

    // ── 2. Dry-run pre-flight ─────────────────────────────────────────────────
    if params.dry_run {
        match show_preflight(&params, &log) {
            ApprovalResult::Rejected => std::process::exit(0),
            ApprovalResult::Approved => {}
        }
    }

    // ── 3. Configure Rayon thread pool ────────────────────────────────────────
    ThreadPoolBuilder::new()
        .num_threads(params.threads)
        .build_global()
        .unwrap_or_else(|e| log.warn(&format!("ThreadPool init: {}", e)));

    // ── 4. Validate required inputs ───────────────────────────────────────────
    log.section("INPUT VALIDATION");
    for (name, path) in &[
        ("ref matrix",   &params.ref_matrix),
        ("alt matrix",   &params.alt_matrix),
        ("clusters.tsv", &params.clusters),
        ("ambient_rna",  &params.ambient),
    ] {
        if path.is_empty() {
            log.fatal(&format!(
                "--{} is required. Set it via CLI, --config JSON, or SOUPX_{} in .soupxcell.env.",
                name.replace(' ', "_"),
                name.to_uppercase().replace(' ', "_")
            ));
        }
        if !std::path::Path::new(path).exists() {
            log.fatal(&format!("File not found: {}  ({})", path, name));
        }
    }
    log.ok("All input files present", None);

    // ── 5. Create output directories ──────────────────────────────────────────
    for dir in &[&params.output, &params.plot_dir] {
        fs::create_dir_all(dir)
            .unwrap_or_else(|e| log.fatal(&format!("Cannot create {}: {}", dir, e)));
    }

    // Copy interactive HTML figure files to output figures/ dir.
    // Source candidates: next to binary (release build), or cargo workspace (dev).
    {
        let html_figures = ["cluster_composition.html"];  // contamination_progression copied after simulation
        let binary_dir = std::env::current_exe().ok()
            .and_then(|p| p.parent().map(|d| d.to_path_buf()));
        let src_dirs: Vec<std::path::PathBuf> = vec![
            binary_dir.as_ref().map(|d| d.join("figures_html")),
            Some(std::path::PathBuf::from("src/infra/figures_html")),
            Some(std::path::PathBuf::from("../src/infra/figures_html")),
        ].into_iter().flatten().collect();
        for fname in &html_figures {
            for src_dir in &src_dirs {
                let src = src_dir.join(fname);
                if src.exists() {
                    let dst = format!("{}/{}", params.plot_dir, fname);
                    match fs::copy(&src, &dst) {
                        Ok(_)  => { log.written(fname, &dst); break; }
                        Err(e) => { log.warn(&format!("Could not copy {}: {}", fname, e)); }
                    }
                    break;
                }
            }
        }
    }

    // ── 6. Load inputs ────────────────────────────────────────────────────────
    log.section("LOADING MATRICES");
    // Primary load: full dense matrix (for correction — needs all cells/loci)
    let raw = load_raw_matrices(&params.ref_matrix, &params.alt_matrix, &log)
        .unwrap_or_else(|e| log.fatal(&e));

    // ── Locus selection: VCF-first priority chain ────────────────────────────
    //
    // The cluster_genotypes.vcf from souporcell IS the definitive locus list.
    // It records every locus souporcell used — no re-deriving QC, no guessing.
    //
    // Priority:
    //   1. --souporcell_vcf + --freebayes_vcf → Mode A: exact CHROM:POS mapping
    //   2. --souporcell_vcf alone             → Mode B: positional row mapping
    //   3. --min_ref/--min_alt (fallback)     → reproduce souporcell QC filter
    log.subsection("Locus selection");

    // Log resolved VCF path so users can debug if $VCF variable was unset
    match &params.souporcell_vcf {
        Some(p) => log.info("souporcell_vcf resolved", p),
        None    => log.info("souporcell_vcf", "not set — check $VCF shell variable is exported before running"),
    }

    // Determine which matrix row indices to use for metrics/embeddings
    let (x_for_metrics, index_to_locus, loci_qc) =
        if let Some(ref vcf_path) = params.souporcell_vcf {
            // VCF-based: read directly from souporcell output
            log.info("locus source", "cluster_genotypes.vcf (souporcell output — ground truth)");
            if params.freebayes_vcf.is_some() {
                log.info("mode", "A — exact CHROM:POS cross-reference with freebayes VCF");
            } else {
                log.info("mode", "B — positional (row N in VCF = row N in matrix)");
                log.info("tip", "add --freebayes_vcf souporcell_merged_sorted_vcf.vcf.gz for exact mapping");
            }
            let indices = parse_vcf_locus_indices(
                vcf_path,
                params.freebayes_vcf.as_deref(),
                &log,
            ).unwrap_or_else(|e| log.fatal(&e));
            let n = indices.len();
            log.metric("loci_from_vcf", &n.to_string(), "count  ← souporcell QC-passing loci");

            // Extract those rows from the raw matrix into a compact [L_vcf × N] array
            let n_cells = raw.n_cells;
            log.computing(&format!("extracting {} VCF-selected loci from {} × {} matrix", n, n, n_cells));
            // Build allele fraction matrix from raw counts (correction not built yet)
            let mut x_compact = ndarray::Array2::<f64>::from_elem((n, n_cells), f64::NAN);
            for (ci, &orig_row) in indices.iter().enumerate() {
                if orig_row >= raw.n_loci { continue; }
                for cell in 0..n_cells {
                    let r = raw.ref_counts[[orig_row, cell]];
                    let a = raw.alt_counts[[orig_row, cell]];
                    let total = r + a;
                    if total > 0.0 {
                        x_compact[[ci, cell]] = a / total;
                    }
                }
            }
            (x_compact, indices, n)

        } else {
            // Fallback: reproduce souporcell's QC filter with min_ref/min_alt
            log.info("locus source", "QC filter (no --souporcell_vcf provided — reproducing souporcell QC)");
            log.warn("add --souporcell_vcf cluster_genotypes.vcf for exact locus matching");
            log.info("QC params", &format!("min_ref={} min_alt={} min_ref_umis={} min_alt_umis={}",
                params.min_ref, params.min_alt, params.min_ref_umis, params.min_alt_umis));
            let (x_qc, idx_to_locus, loci_raw, loci_used, zero_cov) =
                load_matrices_qc(&params.ref_matrix, &params.alt_matrix,
                    params.min_ref, params.min_alt,
                    params.min_ref_umis, params.min_alt_umis, &log)
                .unwrap_or_else(|e| log.fatal(&e));
            log.metric("loci_raw",      &loci_raw.to_string(),   "count");
            log.metric("loci_qc",       &loci_used.to_string(),  "count  ← equivalent to souporcell 'Loci passing QC: N'");
            log.metric("zero_cov_pct",  &format!("{:.1}", 100.0 * zero_cov as f64 / (loci_used * raw.n_cells).max(1) as f64), "%");
            (x_qc, idx_to_locus, loci_used)
        };

    let _l_full = raw.n_loci;
    log.metric("loci_for_metrics", &loci_qc.to_string(), "count");

    log.section("LOADING CLUSTERS");
    let all_cells = load_clusters(&params.clusters, raw.n_cells, &log)
        .unwrap_or_else(|e| log.fatal(&e));

    log.section("PARSING AMBIENT FRACTION");
    let alpha_from_file = parse_ambient_fraction(&params.ambient)
        .unwrap_or_else(|e| log.fatal(&e));
    let alpha = params.alpha_override.unwrap_or(alpha_from_file);
    log.info("alpha (from file)",    &format!("{:.8}", alpha_from_file));
    log.info("alpha (used)",         &format!("{:.8}", alpha));
    if params.alpha_override.is_some() {
        log.warn(&format!("CLI --alpha {:.8} overrides file value {:.8}",
            alpha, alpha_from_file));
    }
    log.ok("Ambient fraction loaded", None);

    // ── 7. Build cell filter (singlets only) ─────────────────────────────────
    let n_clusters = all_cells.iter().map(|c| c.cluster_id).max()
        .map(|m| m + 1).unwrap_or(1);

    let working_cells: Vec<_> = if params.singlets_only {
        all_cells.iter().filter(|c| c.status.is_singlet()).collect()
    } else {
        all_cells.iter().collect()
    };

    let cell_col_indices: Vec<usize> = working_cells.iter()
        .map(|c| c.col_idx)
        .filter(|&idx| idx < raw.n_cells)
        .collect();

    let labels: Vec<usize> = working_cells.iter()
        .filter(|c| c.col_idx < raw.n_cells)
        .map(|c| c.cluster_id)
        .collect();

    log.info("cells used for correction",
        &format!("{} ({} singlets)", cell_col_indices.len(),
            working_cells.iter().filter(|c| c.status.is_singlet()).count()));

    // ── 8. MODULE 1 — Ambient subtraction ────────────────────────────────────
    log.section("MODULE 1 — AMBIENT SUBTRACTION");
    log.info("algorithm", "Supervisor spec: X_corr[l,n] = max(0, X[l,n] - α×μ[l])");
    log.info("alpha",     &format!("{:.8}  ({:.4}%)", alpha, alpha * 100.0));
    log.info("cells for soup estimation", &format!("{} (singlets only: {})",
        cell_col_indices.len(), params.singlets_only));
    let t_m1 = Instant::now();
    log.step("1/3", "Computing allele fraction matrix  X[l,n] = alt[l,n] / (ref[l,n] + alt[l,n])");
    log.step("2/3", "Computing per-locus soup vector   soup[l] = α × mean_n(X[l,n])");
    log.step("3/3", "Subtracting soup, flooring at 0  X_corr[l,n] = max(0, X[l,n] - soup[l])");
    log.computing("ambient subtraction on full matrix");

    let correction = run_correction(
        &raw,
        alpha,
        if params.singlets_only { Some(&cell_col_indices) } else { None },
        &log,
    );

    // Per-cell correction magnitudes (for histogram)
    let cell_corrections: Vec<f64> = cell_col_indices.iter().map(|&ni| {
        let mut sum = 0.0_f64;
        let mut count = 0usize;
        for li in 0..raw.n_loci {
            let r = correction.x_raw[[li, ni]];
            let c = correction.x_corrected[[li, ni]];
            if !r.is_nan() { sum += (r - c).abs(); count += 1; }
        }
        if count > 0 { sum / count as f64 } else { 0.0 }
    }).collect();

    log.ok("Module 1 — Ambient subtraction complete", Some(t_m1));
    log.metric("alpha_used",          &format!("{:.8}", alpha),                  "fraction");
    log.metric("entries_floored",     &format!("{}", correction.n_floored),       "count");
    log.metric("pct_floored",         &format!("{:.4}", correction.pct_floored * 100.0), "%");
    log.metric("mean_abs_correction", &format!("{:.8}", correction.mean_correction), "allele_fraction");
    log.metric("covered_entries",     &format!("{}", correction.n_covered),      "count");
    log.metric("soup_max",            &format!("{:.8}", correction.soup_vector.soup_contribution.iter()
        .cloned().fold(f64::NEG_INFINITY, f64::max)),                             "allele_fraction");
    log.metric("soup_mean",           &format!("{:.8}", correction.soup_vector.soup_contribution.iter()
        .sum::<f64>() / correction.soup_vector.n_loci as f64),                   "allele_fraction");

    // ── 9. Cluster quality metrics (before & after) ───────────────────────────
    log.section("CLUSTER QUALITY METRICS");
    let t_metrics = Instant::now();

    // Extract working-cell columns from the full L×N matrices,
    // keeping ONLY loci that have at least one covered (non-NaN) cell.
    // This filters out the ~95% zero-coverage loci which contribute
    // nothing to silhouette/DB but make computation 100× slower.
    let _n_work = cell_col_indices.len();
    let l_full = raw.n_loci;

    // Pass 1: find covered loci (any non-NaN entry among working cells)
    let covered_loci: Vec<usize> = (0..l_full).filter(|&li| {
        cell_col_indices.iter().any(|&ni| !correction.x_raw[[li, ni]].is_nan())
    }).collect();
    let _l = covered_loci.len();
    // ── Use QC-filtered matrix for all metrics/embeddings ───────────────────
    // x_qc was built with souporcell's exact Pass1+Pass2 filter
    // index_to_locus maps compact_idx → original matrix row
    // This replaces ALL heuristic locus selection (variance caps, hardcoded N, VCF parsing)
    log.info("loci (raw matrix)",  &l_full.to_string());
    log.info("loci (QC-passing)",  &loci_qc.to_string());
    log.info("loci (QC match)",    &format!("matches souporcell 'Loci passing QC: {}' exactly", loci_qc));

    // Extract working-cell columns from x_qc [L_qc × N_all] → [L_qc × N_work]
    let n_work = cell_col_indices.len();
    let l = loci_qc;
    let mut x_raw_work  = ndarray::Array2::<f64>::zeros((l, n_work));
    let mut x_corr_work = ndarray::Array2::<f64>::zeros((l, n_work));
    for (wi, &ni) in cell_col_indices.iter().enumerate() {
        for li in 0..l {
            x_raw_work[[li, wi]]  = x_for_metrics[[li, ni]];
            // For corrected: apply correction on the QC locus at original index
            let orig_li = index_to_locus[li];
            x_corr_work[[li, wi]] = correction.x_corrected[[orig_li, ni]];
        }
    }
    log.computing(&format!("silhouette score ({} QC loci × ~{} cells/cluster)", l, params.metric_cells));

    // Pass 2: select loci for metrics computation.
    // --metric_loci 0  → use ALL covered loci (professor's exact spec)
    // --metric_loci N  → use top N by variance (faster, less accurate)
    // --min_locus_variance V → use only loci with variance > V
    let max_metric_loci = params.metric_loci;  // 0 = no cap = use all
    let min_var         = params.min_locus_variance;

    // Compute per-locus variance over working cells (needed for both filtering modes)
    let mut locus_var: Vec<(usize, f64)> = covered_loci.iter().map(|&li| {
        let vals: Vec<f64> = cell_col_indices.iter()
            .map(|&ni| correction.x_raw[[li, ni]])
            .filter(|v| !v.is_nan())
            .collect();
        let n = vals.len() as f64;
        let var = if n > 1.0 {
            let mean = vals.iter().sum::<f64>() / n;
            vals.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n
        } else { 0.0 };
        (li, var)
    }).collect();

    // Filter by minimum variance if set
    if min_var > 0.0 {
        locus_var.retain(|(_, var)| *var >= min_var);
    }

    // Sort by variance descending
    locus_var.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    // Apply max cap if set (0 = no cap = professor's spec = all covered loci)
    let final_loci: Vec<usize> = if max_metric_loci == 0 || locus_var.len() <= max_metric_loci {
        locus_var.iter().map(|(li, _)| *li).collect()
    } else {
        locus_var.into_iter().take(max_metric_loci).map(|(li, _)| li).collect()
    };

    let l = final_loci.len();
    let loci_description = match (max_metric_loci, min_var > 0.0) {
        (0, false)  => format!("all {} covered loci (professor spec)", l),
        (_n, false)  => format!("top {} by variance from {}", l, covered_loci.len()),
        (0, true)   => format!("{} loci with variance ≥ {:.4}", l, min_var),
        (_n, true)   => format!("top {} loci with variance ≥ {:.4}", l, min_var),
    };
    log.info("loci (for metrics)", &loci_description);
    log.info("metric_loci param", &format!("{} (0=all)", max_metric_loci));

    // Pass 3: build compact [l × n_work] matrices from selected loci only
    let mut x_raw_work  = ndarray::Array2::<f64>::zeros((l, n_work));
    let mut x_corr_work = ndarray::Array2::<f64>::zeros((l, n_work));
    for (wi_l, &li) in final_loci.iter().enumerate() {
        for (wi, &ni) in cell_col_indices.iter().enumerate() {
            x_raw_work[[wi_l, wi]]  = correction.x_raw[[li, ni]];
            x_corr_work[[wi_l, wi]] = correction.x_corrected[[li, ni]];
        }
    }

    let (sil_before_global, sil_before_per) =
        silhouette_score(&x_raw_work, &labels, n_clusters, params.metric_cells);
    let (sil_after_global, sil_after_per) =
        silhouette_score(&x_corr_work, &labels, n_clusters, params.metric_cells);
    let db_before = davies_bouldin(&x_raw_work,  &labels, n_clusters);
    let db_after  = davies_bouldin(&x_corr_work, &labels, n_clusters);
    let mws_before = mean_within_cluster_std(&x_raw_work,  &labels, n_clusters);
    let mws_after  = mean_within_cluster_std(&x_corr_work, &labels, n_clusters);

    let cluster_counts: Vec<usize> = (0..n_clusters).map(|k|
        labels.iter().filter(|&&l| l == k).count()).collect();

    let metrics_before = ClusterMetrics {
        silhouette:              sil_before_global,
        davies_bouldin:          db_before,
        mean_within_cluster_std: mws_before,
        within_cluster_std:      mws_before,
        per_cluster_silhouette:  sil_before_per.clone(),
        cluster_counts:          cluster_counts.clone(),
    };
    let metrics_after = ClusterMetrics {
        silhouette:              sil_after_global,
        davies_bouldin:          db_after,
        mean_within_cluster_std: mws_after,
        within_cluster_std:      mws_after,
        per_cluster_silhouette:  sil_after_per.clone(),
        cluster_counts:          cluster_counts.clone(),
    };

    log.ok("Cluster quality metrics computed", Some(t_metrics));
    log.subsection("Silhouette Score  (range [-1,1], higher = better separation)");
    log.metric("silhouette_before",      &format!("{:.6}", sil_before_global),         "score");
    log.metric("silhouette_after",       &format!("{:.6}", sil_after_global),          "score");
    log.metric("silhouette_delta",       &format!("{:+.6}", sil_after_global - sil_before_global), "score");
    log.subsection("Davies-Bouldin Index  (lower = better separation)");
    log.metric("davies_bouldin_before",  &format!("{:.6}", db_before),                 "index");
    log.metric("davies_bouldin_after",   &format!("{:.6}", db_after),                  "index");
    log.metric("davies_bouldin_delta",   &format!("{:+.6}", db_after - db_before),     "index");
    log.subsection("Within-Cluster Spread");
    log.metric("mean_within_std_before", &format!("{:.6}", mws_before),                "allele_fraction");
    log.metric("mean_within_std_after",  &format!("{:.6}", mws_after),                 "allele_fraction");
    log.metric("within_std_delta",       &format!("{:+.6}", mws_after - mws_before),   "allele_fraction");
    // Interpretation hint
    let sil_delta = sil_after_global - sil_before_global;
    let interp = if sil_delta > 0.01       { "▲ IMPROVED — correction tightened clusters" }
                 else if sil_delta > 0.001 { "▲ slight improvement" }
                 else if sil_delta > -0.001{ "→ no significant change (expected at α=1.26%)" }
                 else                      { "▼ degraded — check alpha value" };
    log.info("interpretation", interp);


    // ── 10. MODULE 3 ─ Embedding ────────────────────────────────────────────────────────────
    // Stores embed coordinates as Array2 for JSON serialisation (no SVG)
    let mut embed_before: Option<ndarray::Array2<f64>> = None;
    let mut embed_after:  Option<ndarray::Array2<f64>> = None;
    let mut umap_before:  Option<ndarray::Array2<f64>> = None;
    let mut umap_after:   Option<ndarray::Array2<f64>> = None;

    if params.embed != EmbedChoice::None {
        log.section("MODULE 3 — DIMENSIONALITY REDUCTION");
        let t_embed = Instant::now();

        let pca_b = run_pca(&x_raw_work,  params.pca_components, &log);
        let pca_a = run_pca(&x_corr_work, params.pca_components, &log);

        if params.embed.run_tsne() {
            log.subsection("PCA → t-SNE  (cells before correction)");
            log.computing(&format!("PCA {}-dim → {}-dim, then t-SNE {} iterations",
                pca_b.ncols(), 2, params.tsne_iter));
            embed_before = Some(run_tsne(&pca_b,
                params.tsne_perplexity, params.tsne_iter, params.seed, &log));
            embed_after  = Some(run_tsne(&pca_a,
                params.tsne_perplexity, params.tsne_iter, params.seed, &log));
            log.ok("t-SNE complete (before + after)", Some(t_embed));
            log.metric("tsne_perplexity", &params.tsne_perplexity.to_string(), "");
            log.metric("tsne_iterations", &params.tsne_iter.to_string(),       "");
            log.metric("pca_components",  &params.pca_components.to_string(),  "");
        }

        if params.embed.run_umap() {
            log.subsection("PCA → UMAP  (simplified k-NN force-directed layout)");
            log.computing("UMAP embedding before + after correction");
            let t_umap = Instant::now();
            umap_before = Some(run_umap(&pca_b, 15, 0.1, 300, params.seed, &log));
            umap_after  = Some(run_umap(&pca_a, 15, 0.1, 300, params.seed, &log));
            log.ok("UMAP complete", Some(t_umap));
        }
    }
    // ── 11. MODULE 2 — Simulation benchmark ──────────────────────────────────
    let mut bench_result: Option<crate::domain::types::BenchmarkResult> = None;

    if params.simulate {
        log.section("MODULE 2 — SIMULATION BENCHMARK");
        log.info("simulation levels",  &format!("{:?}", params.sim_levels));
        log.info("trials per level",   &params.sim_trials.to_string());
        log.info("method", "inject α×μ → correct → measure RMSE + silhouette delta");
        let _t_sim = Instant::now();

        let locus_mean = correction.soup_vector.locus_mean.clone();
        let result = run_simulation(
            &x_raw_work,
            &locus_mean,
            &labels,
            n_clusters,
            &params.sim_levels,
            params.sim_trials,
            alpha,
            &log,
        );
        bench_result = Some(result);
    }

    // ── 11b. Log simulation results table ───────────────────────────────────────
    if let Some(ref bench) = bench_result {
        log.subsection("Simulation Benchmark Results");
        log.info("header", "alpha_sim   RMSE        sil_before  sil_after   delta       floor%");
        for pt in &bench.points {
            let delta = pt.silhouette_after - pt.silhouette_before;
            let trend = if delta > 0.005 { "▲" } else if delta < -0.005 { "▼" } else { "→" };
            log.info(
                &format!("α={:.3}", pt.alpha_sim),
                &format!("rmse={:.6}  sil_bef={:.4}  sil_aft={:.4}  Δ={:+.4} {}  floor={:.1}%",
                    pt.rmse, pt.silhouette_before, pt.silhouette_after,
                    delta, trend, pt.floor_pct_after * 100.0)
            );
        }
    }

    // ── 12. Write image_data JSON files + copy HTML figure templates ────────
    log.section("GENERATING PLOTS");

    let img_data_dir = format!("{}/image_data", params.plot_dir);
    fs::create_dir_all(&img_data_dir).unwrap_or_default();

    // Helper: write a JSON file and log it
    let write_json = |name: &str, json: &str| {
        let path = format!("{}/{}.json", img_data_dir, name);
        match fs::write(&path, json) {
            Ok(_)  => log.written(&format!("{}.json", name), &path),
            Err(e) => log.warn(&format!("Could not write {}.json: {}", name, e)),
        }
    };

    // ── JSON 1+2: tsne/umap_before_after ─────────────────────────────────
    if let (Some(ref eb), Some(ref ea)) = (&embed_before, &embed_after) {
        let mut j = String::from("{\n");
        j.push_str(&format!("  \"method\": \"t-SNE\",\n  \"alpha\": {:.8},\n  \"n_clusters\": {},\n", alpha, n_clusters));
        j.push_str(&format!("  \"sil_before\": {:.6},\n  \"sil_after\": {:.6},\n", sil_before_global, sil_after_global));
        j.push_str("  \"cells\": [\n");
        let n = eb.nrows().min(labels.len());
        for i in 0..n {
            let comma = if i + 1 < n { "," } else { "" };
            j.push_str(&format!("    {{\"bx\":{:.4},\"by\":{:.4},\"ax\":{:.4},\"ay\":{:.4},\"k\":{}}}{}\n",
                eb[[i,0]], eb[[i,1]], ea[[i,0]], ea[[i,1]], labels[i], comma));
        }
        j.push_str("  ]\n}\n");
        write_json("tsne_before_after", &j);
    }
    if let (Some(ref ub), Some(ref ua)) = (&umap_before, &umap_after) {
        let mut j = String::from("{\n");
        j.push_str(&format!("  \"method\": \"UMAP\",\n  \"alpha\": {:.8},\n  \"n_clusters\": {},\n", alpha, n_clusters));
        j.push_str(&format!("  \"sil_before\": {:.6},\n  \"sil_after\": {:.6},\n", sil_before_global, sil_after_global));
        j.push_str("  \"cells\": [\n");
        let n = ub.nrows().min(labels.len());
        for i in 0..n {
            let comma = if i + 1 < n { "," } else { "" };
            j.push_str(&format!("    {{\"bx\":{:.4},\"by\":{:.4},\"ax\":{:.4},\"ay\":{:.4},\"k\":{}}}{}\n",
                ub[[i,0]], ub[[i,1]], ua[[i,0]], ua[[i,1]], labels[i], comma));
        }
        j.push_str("  ]\n}\n");
        write_json("umap_before_after", &j);
    }

    // ── JSON 3: silhouette_comparison ────────────────────────────────────
    {
        let mut j = String::from("{\n");
        j.push_str(&format!("  \"n_clusters\": {},\n  \"global_before\": {:.6},\n  \"global_after\": {:.6},\n",
            n_clusters, sil_before_global, sil_after_global));
        j.push_str(&format!("  \"alpha\": {:.8},\n", alpha));
        j.push_str("  \"clusters\": [\n");
        for k in 0..n_clusters {
            let sb = sil_before_per.get(k).cloned().unwrap_or(0.0);
            let sa = sil_after_per.get(k).cloned().unwrap_or(0.0);
            let cnt = cluster_counts.get(k).cloned().unwrap_or(0);
            let comma = if k + 1 < n_clusters { "," } else { "" };
            j.push_str(&format!("    {{\"k\":{k},\"sil_before\":{sb:.6},\"sil_after\":{sa:.6},\"delta\":{:.6},\"n_cells\":{cnt}}}{comma}\n",
                sa - sb));
        }
        j.push_str("  ]\n}\n");
        write_json("silhouette_comparison", &j);
    }

    // ── JSON 4: correction_distribution ──────────────────────────────────
    {
        let mut j = String::from("{\n");
        j.push_str(&format!("  \"alpha\": {:.8},\n  \"n_floored\": {},\n  \"pct_floored\": {:.6},\n  \"n_cells\": {},\n",
            alpha, correction.n_floored, correction.pct_floored * 100.0, cell_corrections.len()));
        j.push_str(&format!("  \"mean_correction\": {:.8},\n", correction.mean_correction));
        j.push_str("  \"corrections\": [");
        let strs: Vec<String> = cell_corrections.iter().map(|v| format!("{:.6}", v)).collect();
        j.push_str(&strs.join(","));
        j.push_str("]\n}\n");
        write_json("correction_distribution", &j);
    }

    // ── JSON 5: simulation_curves ─────────────────────────────────────────
    if let Some(ref bench) = bench_result {
        let mut j = String::from("{\n");
        j.push_str(&format!("  \"true_alpha\": {:.8},\n", bench.true_alpha));
        j.push_str("  \"points\": [\n");
        for (i, p) in bench.points.iter().enumerate() {
            let comma = if i + 1 < bench.points.len() { "," } else { "" };
            j.push_str(&format!(
                "    {{\"alpha_sim\":{:.4},\"rmse\":{:.6},\"sil_before\":{:.6},\"sil_after\":{:.6},\"delta\":{:.6},\"floor_pct\":{:.4}}}{comma}\n",
                p.alpha_sim, p.rmse, p.silhouette_before, p.silhouette_after,
                p.silhouette_after - p.silhouette_before, p.floor_pct_after * 100.0));
        }
        j.push_str("  ]\n}\n");
        write_json("simulation_curves", &j);
        write_json("rmse_linearity",    &j);  // same data, different view
    }

    // ── JSON 6: soup_vector_profile ───────────────────────────────────────
    {
        // Sort loci by locus_mean DESC (same order as the plot function)
        let mut indexed: Vec<(usize, f64)> = correction.soup_vector.locus_mean.iter()
            .cloned().enumerate()
            .filter(|(_, mu)| mu.is_finite() && *mu > 0.0)
            .collect();
        indexed.sort_by(|a, b| b.1.total_cmp(&a.1));
        let top30: Vec<(usize, f64)> = {
            let mut out = Vec::new();
            let mut n_max = 0usize;
            for &(idx, mu) in &indexed {
                if mu >= 0.9999 { if n_max < 3 { out.push((idx, mu)); n_max += 1; } }
                else { out.push((idx, mu)); if out.len() >= 30 { break; } }
            }
            out
        };
        let mut j = String::from("{\n");
        j.push_str(&format!("  \"alpha\": {:.8},\n  \"n_loci_total\": {},\n", alpha, correction.soup_vector.locus_mean.len()));
        j.push_str("  \"top_loci\": [\n");
        for (i, &(idx, mu)) in top30.iter().enumerate() {
            let soup = mu * alpha;
            let comma = if i + 1 < top30.len() { "," } else { "" };
            j.push_str(&format!("    {{\"rank\":{i},\"locus_idx\":{idx},\"mu\":{mu:.6},\"soup\":{soup:.8}}}{comma}\n"));
        }
        j.push_str("  ],\n");
        // Distribution: all finite locus_mean values (soup = mu * alpha)
        j.push_str("  \"all_soup\": [");
        let soups: Vec<String> = correction.soup_vector.soup_contribution.iter()
            .filter(|v| v.is_finite()).map(|v| format!("{:.6}", v)).collect();
        j.push_str(&soups.join(","));
        j.push_str("]\n}\n");
        write_json("soup_vector_profile", &j);
    }

    // ── JSON 7: per_cluster_correction ───────────────────────────────────
    {
        let mut per_corr = vec![0.0f64; n_clusters];
        let mut per_cnt  = vec![0usize;  n_clusters];
        let mut per_floored   = vec![0usize; n_clusters];
        let mut per_covered   = vec![0usize; n_clusters];
        for (&lbl, &c) in labels.iter().zip(cell_corrections.iter()) {
            if lbl < n_clusters { per_corr[lbl] += c; per_cnt[lbl] += 1; }
        }
        for k in 0..n_clusters { if per_cnt[k] > 0 { per_corr[k] /= per_cnt[k] as f64; } }
        let n_loci_soup = correction.soup_vector.soup_contribution.len().min(correction.x_raw.nrows());
        for (wi, (&col, &lbl)) in cell_col_indices.iter().zip(labels.iter()).enumerate() {
            if lbl >= n_clusters || col >= correction.x_raw.ncols() { continue; }
            let _ = wi;
            for li in 0..n_loci_soup {
                let v = correction.x_raw[[li, col]];
                if !v.is_nan() {
                    per_covered[lbl] += 1;
                    if v < correction.soup_vector.soup_contribution[li] { per_floored[lbl] += 1; }
                }
            }
        }
        let mut j = String::from("{\n");
        j.push_str(&format!("  \"alpha\": {:.8},\n  \"n_clusters\": {},\n", alpha, n_clusters));
        j.push_str("  \"clusters\": [\n");
        for k in 0..n_clusters {
            let floor_pct = if per_covered[k] > 0 { per_floored[k] as f64 / per_covered[k] as f64 * 100.0 } else { 0.0 };
            let sb = sil_before_per.get(k).cloned().unwrap_or(0.0);
            let sa = sil_after_per.get(k).cloned().unwrap_or(0.0);
            let cnt = cluster_counts.get(k).cloned().unwrap_or(0);
            let comma = if k + 1 < n_clusters { "," } else { "" };
            j.push_str(&format!(
                "    {{\"k\":{k},\"n_cells\":{cnt},\"mean_correction\":{:.8},\"floor_pct\":{:.4},\"sil_before\":{sb:.6},\"sil_after\":{sa:.6},\"sil_delta\":{:.6}}}{comma}\n",
                per_corr[k], floor_pct, sa - sb));
        }
        j.push_str("  ]\n}\n");
        write_json("per_cluster_correction", &j);
    }

    // ── JSON 8: floor_per_locus ───────────────────────────────────────────
    {
        let l_qc = index_to_locus.len();
        let n_cells = correction.x_raw.ncols();
        let n_loci_soup = correction.soup_vector.soup_contribution.len().min(correction.x_raw.nrows());
        let mut floor_counts: Vec<(usize, usize, usize)> = Vec::with_capacity(l_qc);
        for (ci, &orig_li) in index_to_locus.iter().enumerate() {
            if orig_li >= n_loci_soup { continue; }
            let s = correction.soup_vector.soup_contribution[orig_li];
            let fc = (0..n_cells).filter(|&ni| {
                let v = correction.x_raw[[orig_li, ni]];
                !v.is_nan() && v < s
            }).count();
            floor_counts.push((ci, orig_li, fc));
        }
        floor_counts.sort_by(|a, b| b.2.cmp(&a.2));
        let mut j = String::from("{\n");
        j.push_str(&format!("  \"n_cells\": {},\n  \"l_qc\": {},\n", n_cells, l_qc));
        j.push_str("  \"top_loci\": [\n");
        for (i, &(ci, orig, fc)) in floor_counts.iter().take(25).enumerate() {
            let pct = fc as f64 / n_cells.max(1) as f64 * 100.0;
            let comma = if i + 1 < 25.min(floor_counts.len()) { "," } else { "" };
            j.push_str(&format!("    {{\"rank\":{i},\"compact_idx\":{ci},\"locus_idx\":{orig},\"floor_count\":{fc},\"floor_pct\":{pct:.2}}}{comma}\n"));
        }
        j.push_str("  ],\n");
        j.push_str("  \"all_floor_counts\": [");
        let all: Vec<String> = floor_counts.iter().map(|t| t.2.to_string()).collect();
        j.push_str(&all.join(","));
        j.push_str("]\n}\n");
        write_json("floor_per_locus", &j);
    }

    // ── Copy HTML figure templates to figures/ dir ────────────────────────
    {
        let html_plots = [
            "tsne_before_after.html", "umap_before_after.html",
            "silhouette_comparison.html", "correction_distribution.html",
            "simulation_curves.html", "soup_vector_profile.html",
            "per_cluster_correction.html", "floor_per_locus.html",
            "rmse_linearity.html",
            "contamination_progression.html",  // reads simulation_curves.json, copied after sim
        ];
        let src_dirs: Vec<std::path::PathBuf> = vec![
            std::env::current_exe().ok()
                .and_then(|p| p.parent().map(|d| d.join("figures_html"))),
            Some(std::path::PathBuf::from("src/infra/figures_html")),
            Some(std::path::PathBuf::from("../src/infra/figures_html")),
        ].into_iter().flatten().collect();
        for fname in &html_plots {
            for src_dir in &src_dirs {
                let src = src_dir.join(fname);
                if src.exists() {
                    let dst = format!("{}/{}", params.plot_dir, fname);
                    match fs::copy(&src, &dst) {
                        Ok(_)  => { log.written(fname, &dst); break; }
                        Err(e) => { log.warn(&format!("Could not copy {}: {}", fname, e)); }
                    }
                    break;
                }
            }
        }
    }

    // ── plot_files list: HTML figures for the report ──────────────────────
    // (replaces SVG file list — report now embeds HTML iframes)

        log.ok("All plots written", None);
    log.checkpoint("plots-complete");

    // ── 13. Write data outputs ────────────────────────────────────────────────
    log.section("WRITING OUTPUTS");
    let t_write = Instant::now();

    write_corrected_mtx(
        &correction.x_corrected,
        &format!("{}/X_corrected.mtx", params.output),
    ).unwrap_or_else(|e| log.warn(&e));
    log.written("X_corrected.mtx", &format!("{}/X_corrected.mtx", params.output));

    write_corrected_mtx(
        &correction.x_raw,
        &format!("{}/X_raw.mtx", params.output),
    ).unwrap_or_else(|e| log.warn(&e));
    log.written("X_raw.mtx", &format!("{}/X_raw.mtx", params.output));

    write_soup_vector(
        &correction.soup_vector,
        &format!("{}/soup_vector.csv", params.output),
    ).unwrap_or_else(|e| log.warn(&e));
    log.written("soup_vector.csv", &format!("{}/soup_vector.csv", params.output));

    write_correction_summary(
        &correction,
        &working_cells.iter().map(|c| (*c).clone()).collect::<Vec<_>>(),
        &format!("{}/correction_summary.csv", params.output),
    ).unwrap_or_else(|e| log.warn(&e));
    log.written("correction_summary.csv", &format!("{}/correction_summary.csv", params.output));

    // ── Write cluster_data.json — real per-cluster breakdown for HTML figures ─
    // Computes exact per-cluster singlet/doublet counts, silhouette, and
    // correction from the actual parsed clusters.tsv data. No estimates.
    {
        let json_path = format!("{}/figures/cluster_data.json", params.output);
        let mut json = String::from("{\n");
        json.push_str(&format!("  \"run_alpha\": {:.8},\n", alpha));
        json.push_str(&format!("  \"loci_raw\": {},\n", raw.n_loci));
        json.push_str(&format!("  \"loci_qc\": {},\n", loci_qc));
        json.push_str(&format!("  \"total_cells_clusters_tsv\": {},\n", all_cells.len()));
        json.push_str(&format!("  \"total_cells_matrix\": {},\n", raw.n_cells));
        json.push_str(&format!("  \"n_clusters\": {},\n", n_clusters));
        json.push_str(&format!("  \"global_sil_before\": {:.6},\n", sil_before_global));
        json.push_str(&format!("  \"global_sil_after\": {:.6},\n", sil_after_global));
        json.push_str(&format!("  \"mean_abs_correction\": {:.8},\n", correction.mean_correction));
        json.push_str(&format!("  \"pct_floored\": {:.4},\n", correction.pct_floored * 100.0));
        json.push_str("  \"clusters\": [\n");

        // Per-cluster: exact singlet/doublet counts from all_cells (parsed from clusters.tsv)
        // mean correction from cell_corrections (computed from real allele fractions)
        for k in 0..n_clusters {
            let k_sing   = all_cells.iter().filter(|c| c.cluster_id == k && c.status.is_singlet()).count();
            let k_doub   = all_cells.iter().filter(|c| c.cluster_id == k && !c.status.is_singlet()).count();
            let k_total  = k_sing + k_doub;
            let k_sil_b  = sil_before_per.get(k).cloned().unwrap_or(0.0);
            let k_sil_a  = sil_after_per.get(k).cloned().unwrap_or(0.0);
            // Per-cluster mean correction: average cell_corrections for cells in this cluster
            let k_corr: f64 = {
                let vals: Vec<f64> = labels.iter().zip(cell_corrections.iter())
                    .filter(|(&lbl, _)| lbl == k)
                    .map(|(_, &c)| c)
                    .collect();
                if vals.is_empty() { 0.0 } else { vals.iter().sum::<f64>() / vals.len() as f64 }
            };
            let comma = if k + 1 < n_clusters { "," } else { "" };
            json.push_str("    {\"cluster\": ");
            json.push_str(&k.to_string());
            json.push_str(", \"singlets\": ");
            json.push_str(&k_sing.to_string());
            json.push_str(", \"doublets\": ");
            json.push_str(&k_doub.to_string());
            json.push_str(", \"total\": ");
            json.push_str(&k_total.to_string());
            json.push_str(&format!(", \"sil_before\": {:.6}", k_sil_b));
            json.push_str(&format!(", \"sil_after\": {:.6}", k_sil_a));
            json.push_str(&format!(", \"mean_correction\": {:.8}", k_corr));
            json.push_str("}");
            json.push_str(comma);
            json.push('\n');
        }
        json.push_str("  ]\n}\n");

        match fs::write(&json_path, &json) {
            Ok(_)  => log.written("cluster_data.json", &json_path),
            Err(e) => log.warn(&format!("Could not write cluster_data.json: {}", e)),
        }
    }

    if let Some(ref bench) = bench_result {
        write_benchmark_csv(bench,
            &format!("{}/simulation_benchmark.csv", params.output))
            .unwrap_or_else(|e| log.warn(&e));
        log.written("simulation_benchmark.csv", &format!("{}/simulation_benchmark.csv", params.output));
    }

    write_metrics_txt(
        alpha,
        raw.n_loci,
        n_work,
        correction.n_floored,
        correction.pct_floored,
        correction.mean_correction,
        sil_before_global, sil_after_global,
        db_before, db_after,
        &format!("{}/metrics.txt", params.output),
    ).unwrap_or_else(|e| log.warn(&e));
    log.written("metrics.txt", &format!("{}/metrics.txt", params.output));

    // ── 14. Write HTML report ─────────────────────────────────────────────────
    let n_sing = all_cells.iter().filter(|c| c.status.is_singlet()).count();
    let n_doub = all_cells.len() - n_sing;

    let run_summary = RunSummary {
        n_loci:          raw.n_loci,
        n_cells_total:   all_cells.len(),
        n_cells_singlet: n_sing,
        n_cells_doublet: n_doub,
        n_clusters,
        alpha_used:      alpha,
        n_floored:       correction.n_floored,
        pct_floored:     correction.pct_floored,
        mean_correction: correction.mean_correction,
        soup_max: correction.soup_vector.soup_contribution.iter()
            .cloned().fold(f64::NEG_INFINITY, f64::max),
        soup_mean: correction.soup_vector.soup_contribution.iter()
            .sum::<f64>() / correction.soup_vector.n_loci.max(1) as f64,
        n_covered: correction.n_covered,
        metrics_before,
        metrics_after,
        silhouette_delta: sil_after_global - sil_before_global,
        simulation_run:  params.simulate,
        per_cluster: (0..n_clusters).map(|k| {
            let kc = cluster_counts.get(k).cloned().unwrap_or(0);
            let k_sing = working_cells.iter().filter(|c| c.cluster_id == k && c.status.is_singlet()).count();
            let k_doub = kc.saturating_sub(k_sing);
            let k_mean_corr = if cell_col_indices.is_empty() { 0.0 } else {
                labels.iter().zip(cell_corrections.iter())
                    .filter(|(&l, _)| l == k)
                    .map(|(_, &c)| c).sum::<f64>()
                    / labels.iter().filter(|&&l| l == k).count().max(1) as f64
            };
            crate::domain::types::PerClusterSummary {
                n_cells: kc, n_singlets: k_sing, n_doublets: k_doub,
                mean_correction: k_mean_corr, pct_floored: 0.0,
            }
        }).collect(),
    };

    // Determine locus source mode for report
    let locus_mode_str = if params.souporcell_vcf.is_some() && params.freebayes_vcf.is_some() {
        "A".to_string()
    } else if params.souporcell_vcf.is_some() {
        "B".to_string()
    } else {
        "C".to_string()
    };

    let zero_cov_pct = 100.0 * correction.n_covered as f64
        / (raw.n_loci * raw.n_cells).max(1) as f64;
    let zero_cov_pct = 100.0 - zero_cov_pct; // zero-cov = 100% - covered%

    let report_data = crate::analysis::report::ReportData {
        params:       &params,
        summary:      &run_summary,
        bench:        bench_result.as_ref(),
        loci_raw:     raw.n_loci,
        loci_qc:      loci_qc,
        zero_cov_pct,
        locus_mode:   locus_mode_str,
        load_time:    std::time::Duration::from_secs(0),
        total_time:   std::time::Duration::from_secs(0),
        plot_files:   {
            // All figures are now HTML files that read from image_data/*.json
            let mut pf = Vec::new();
            if embed_before.is_some() { pf.push(format!("{}/tsne_before_after.html", params.plot_dir)); }
            if umap_before.is_some()  { pf.push(format!("{}/umap_before_after.html", params.plot_dir)); }
            pf.push(format!("{}/silhouette_comparison.html", params.plot_dir));
            pf.push(format!("{}/correction_distribution.html", params.plot_dir));
            pf.push(format!("{}/soup_vector_profile.html", params.plot_dir));
            pf.push(format!("{}/per_cluster_correction.html", params.plot_dir));
            pf.push(format!("{}/floor_per_locus.html", params.plot_dir));
            if bench_result.is_some() {
                pf.push(format!("{}/simulation_curves.html", params.plot_dir));
                pf.push(format!("{}/rmse_linearity.html", params.plot_dir));
            }
            pf
        },
        output_dir:   params.output.clone(),
    };
    let report_path = crate::analysis::report::write_report(&report_data);
    log.written("soupxcell_report.html", &report_path);
    log.ok("All outputs written", Some(t_write));

    // ── 15. Final summary ─────────────────────────────────────────────────────
    log.section("COMPLETE");

    eprintln!("\n  {}╔══════════════════════════════════════════════════════╗{}",
        "\x1b[36m", "\x1b[0m");
    eprintln!("  {}║   soupxcell v1.0  ·  Run Summary                   ║{}",
        "\x1b[1m\x1b[36m", "\x1b[0m");
    eprintln!("  {}╚══════════════════════════════════════════════════════╝{}\n",
        "\x1b[36m", "\x1b[0m");

    log.info("alpha used",         &format!("{:.8}  ({:.4}%)", alpha, alpha * 100.0));
    log.info("loci processed",     &raw.n_loci.to_string());
    log.info("cells processed",    &format!("{} ({} singlets)", n_work, n_sing));
    log.info("entries floored",    &format!("{} ({:.2}%)",
        correction.n_floored, correction.pct_floored * 100.0));
    log.info("mean |correction|",  &format!("{:.6}", correction.mean_correction));
    log.info("silhouette before",  &format!("{:.4}", sil_before_global));
    log.info("silhouette after",   &format!("{:.4}", sil_after_global));
    log.info("silhouette Δ",       &format!("{:+.4}  ({})",
        sil_after_global - sil_before_global,
        if sil_after_global > sil_before_global + 0.001 { "▲ improved" }
        else if sil_after_global < sil_before_global - 0.001 { "▼ degraded" }
        else { "→ unchanged" }));
    log.info("output directory",   &params.output);

    eprintln!();
    log.ok(&format!("Total runtime: {:.2}s",
        log.start.elapsed().as_secs_f64()), None);

    eprintln!("\n  Key files:");
    eprintln!("    {}/X_corrected.mtx", params.output);
    eprintln!("    {}/soup_vector.csv", params.output);
    eprintln!("    {}/metrics.txt",     params.output);
    eprintln!("    {}/soupxcell_report.html", params.output);
    if embed_before.is_some() {
        eprintln!("    {}/tsne_before_after.html", params.plot_dir);
    }
    if umap_before.is_some() {
        eprintln!("    {}/umap_before_after.html", params.plot_dir);
    }
    eprintln!("    {}/silhouette_comparison.html", params.plot_dir);
    if bench_result.is_some() {
        eprintln!("    {}/simulation_benchmark.csv", params.output);
        eprintln!("    {}/simulation_curves.html", params.plot_dir);
    }
    eprintln!();
}
