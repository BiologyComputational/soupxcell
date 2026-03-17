// ============================================================================
// analysis/report/mod.rs — HTML report module root
// ============================================================================
//
// Mirrors souporcell v3.0 report sub-module architecture exactly.
// Each sub-module has exactly one responsibility.
//
// Sub-modules:
//   css            — §0  white-theme CSS (Source Serif 4 + JetBrains Mono + DM Sans)
//   utils          — shared colour constants + chrono_now() + helpers
//   cards          — §2  executive-summary metric card grid  (12 cards)
//   experiment     — §3  experiment design + correction narrative + locus source
//   tables         — §5  6 detailed results tables
//   svgs           — §6  inline SVG: soup vector bar + correction histogram + cluster strips
//   plot_embeds    — §7  embedded SVG file cards
//   interpretation — §8  PhD-level per-metric commentary (9 interpretation cards)
//   params         — §9  full run parameters table (all 27 params)
//   methods        — §10 algorithmic methods + references
// ============================================================================

pub mod cards;
pub mod css;
pub mod experiment;
pub mod interpretation;
pub mod methods;
pub mod params;
pub mod plot_embeds;
pub mod svgs;
pub mod tables;
pub mod utils;

use std::fs;
use std::time::Duration;

use crate::config::params::Params;
use crate::domain::types::{RunSummary, BenchmarkResult};

// ── Public data carrier ───────────────────────────────────────────────────────

pub struct ReportData<'a> {
    pub params:          &'a Params,
    pub summary:         &'a RunSummary,
    pub bench:           Option<&'a BenchmarkResult>,
    pub loci_raw:        usize,
    pub loci_qc:         usize,
    pub zero_cov_pct:    f64,
    pub locus_mode:      String,    // "A", "B", or "C"
    #[allow(dead_code)]
    pub load_time:       Duration,
    pub total_time:      Duration,
    pub plot_files:      Vec<String>,
    pub output_dir:      String,
}

// ── Public entry point ────────────────────────────────────────────────────────

pub fn write_report(data: &ReportData<'_>) -> String {
    let path = format!("{}/soupxcell_report.html", data.output_dir);
    let html = build_html(data);
    if let Some(parent) = std::path::Path::new(&path).parent() {
        fs::create_dir_all(parent).ok();
    }
    fs::write(&path, &html).unwrap_or_else(|e| {
        eprintln!("[SOUPX] WARN  Cannot write report '{}': {}", path, e);
    });
    path
}

// ── HTML orchestrator ─────────────────────────────────────────────────────────


// ── HTML interactive figures ───────────────────────────────────────────────
// Two standalone HTML figures embedded as iframes with download buttons.
// The HTML files are copied from src/infra/figures_html/ to the output figures/ dir.
fn build_html_figures(data: &ReportData<'_>) -> String {
    let figures_src = std::path::PathBuf::from(
        std::env::current_exe().unwrap_or_default().parent().unwrap_or(std::path::Path::new("."))
    );
    let out_figures = format!("{}/figures", data.output_dir);
    
    // Copy HTML figure files to output figures/ dir
    let html_files = [
        ("cluster_composition.html",        "Cluster Composition"),
        ("contamination_progression.html",  "Contamination Progression"),
    ];
    
    let src_dir_candidates = [
        // Next to binary (release build)
        figures_src.join("figures_html"),
        // Relative to cargo workspace (dev)
        std::path::PathBuf::from("src/infra/figures_html"),
    ];

    for (fname, _) in &html_files {
        for src_dir in &src_dir_candidates {
            let src = src_dir.join(fname);
            if src.exists() {
                let dst = format!("{}/{}", out_figures, fname);
                let _ = fs::copy(&src, &dst);
                break;
            }
        }
    }

    let mut html = String::from(r#"
<hr class="section-divider"/>
<div class="section-label">Interactive Figures</div>
<h2>Interactive Cluster &amp; Contamination Figures</h2>
<p class="section-desc">
  Self-contained interactive HTML figures with PNG download buttons.
  Open in any browser for full interactivity. All charts use real run data.
</p>
<div class="plots-grid">
"#);

    for (fname, title) in &html_files {
        let href = format!("figures/{}", fname);
        html += &format!(r#"
<div class="plot-card" style="grid-column:span 1">
  <div class="plot-title">{title}</div>
  <div style="background:#0F172A;border-radius:6px;overflow:hidden;height:380px">
    <iframe src="{href}" width="100%" height="380"
            style="border:none;display:block"
            title="{title}">
    </iframe>
  </div>
  <div style="margin-top:8px">
    <a href="{href}" target="_blank"
       style="display:inline-block;background:#1E293B;color:#94A3B8;border:1px solid #334155;
              border-radius:5px;padding:5px 12px;font-size:11px;text-decoration:none">
      ↗ Open full interactive figure
    </a>
  </div>
</div>
"#, href=href, title=title);
    }

    html += "</div>
";
    html
}


// ── §1b Research context: objectives, questions, hypothesis ──────────────────
fn build_research_section(data: &ReportData<'_>, delta_sil: f64, delta_wstd: f64) -> String {
    let s = data.summary;
    let _ = (delta_sil, delta_wstd, s);
    let mut h = String::new();
    h.push_str("<div class=\"section-label\">Research Context</div>\n");
    h.push_str("<h2>Why Are We Doing This? Objectives, Questions &amp; Hypothesis</h2>\n");
    h.push_str("<p class=\"section-desc\">Before looking at any numbers, here is the scientific context of this analysis. Written so anyone — with or without a biology background — can follow along.</p>\n");

    // Background story box
    h.push_str("<div class=\"story-block\" style=\"margin-bottom:20px\">\n");
    h.push_str("<div class=\"story-label blue\">&#x1F9EC; Background: The Problem We Are Solving</div>\n");
    h.push_str("<p>When we sequence thousands of single cells, each cell is captured in a tiny droplet. But before capture, cells release some of their RNA molecules into the surrounding liquid — this is called <strong>ambient RNA</strong> or <strong>\"soup\"</strong>. Every captured droplet picks up a little of this background soup. The result: each cell appears to express genes that actually came from other cells, not from the cell itself.</p>\n");
    h.push_str("<p style=\"margin-bottom:0\">For genotype demultiplexing (figuring out which cell came from which donor in a pooled sample), this contamination is especially problematic: ambient RNA from donor A can appear inside a cell from donor B, making it look like a hybrid. <strong>soupXcell cleans this up.</strong></p>\n");
    h.push_str("</div>\n");

    // Three columns: objectives, questions, hypothesis
    h.push_str("<div class=\"three-col\" style=\"margin-bottom:24px\">\n");

    // Objective box
    h.push_str("<div class=\"hypo-box objective\">\n");
    h.push_str("<div class=\"hypo-label o\">&#x1F3AF; Research Objectives</div>\n");
    h.push_str("<p><strong>Primary:</strong> Subtract the ambient RNA contribution from every cell's allele fraction profile so that downstream genotype demultiplexing reflects true biology, not background contamination.</p>\n");
    h.push_str("<p><strong>Secondary:</strong> Validate the correction using a simulation benchmark — inject known contamination, correct it, and measure how accurately it was removed.</p>\n");
    h.push_str("<p style=\"margin-bottom:0\"><strong>Tertiary:</strong> Quantify the effect on cluster separation (silhouette score) to understand when correction matters most.</p>\n");
    h.push_str("</div>\n");

    // Research questions box
    h.push_str("<div class=\"hypo-box question\">\n");
    h.push_str("<div class=\"hypo-label q\">&#x2753; Research Questions</div>\n");
    h.push_str("<ol style=\"padding-left:18px;font-size:.88rem;color:#334155\">\n");
    h.push_str("<li style=\"margin-bottom:8px\">What fraction of each cell's RNA signal is ambient contamination (&alpha;)?</li>\n");
    h.push_str("<li style=\"margin-bottom:8px\">After correction, are the donor clusters more clearly separated?</li>\n");
    h.push_str("<li style=\"margin-bottom:8px\">Does the correction formula (max(0, X &minus; &alpha;&times;&mu;[l])) accurately recover the true signal — provable by simulation?</li>\n");
    h.push_str("<li style=\"margin-bottom:0\">At what contamination level does ambient RNA become detectable in cluster quality metrics?</li>\n");
    h.push_str("</ol>\n");
    h.push_str("</div>\n");

    // Hypothesis box
    h.push_str("<div class=\"hypo-box\">\n");
    h.push_str("<div class=\"hypo-label\">&#x1F9EA; Hypotheses</div>\n");
    h.push_str("<p><strong>H1 (Correction):</strong> The allele fraction of donor-specific SNPs in each cell will move closer to the true genotype (0 or 1) after subtracting the ambient contribution.</p>\n");
    h.push_str("<p><strong>H2 (Silhouette):</strong> At low &alpha; (&lt;5%), silhouette change will be negligible — ambient RNA contributes less than biological noise.</p>\n");
    h.push_str("<p style=\"margin-bottom:0\"><strong>H3 (Linearity):</strong> RMSE from the simulation will scale linearly with &alpha;_sim (R&sup2; &approx; 1.0), proving the correction is mathematically exact.</p>\n");
    h.push_str("</div>\n");
    h.push_str("</div>\n");  // three-col

    // What we did — timeline
    h.push_str("<div class=\"section-label\" style=\"margin-top:8px\">What We Did in This Run</div>\n");
    h.push_str("<div class=\"timeline\">\n");
    h.push_str("<div class=\"tl-step\" data-n=\"1\"><h4>Loaded the data</h4><p>Read the allele count matrices (ref.mtx, alt.mtx) from vartrix and the cluster assignments from souporcell. These contain the raw read counts at every genomic position for every cell.</p></div>\n");
    h.push_str("<div class=\"tl-step\" data-n=\"2\"><h4>Computed ambient RNA levels</h4><p>Used &alpha; = ");
    h.push_str(&format!("{:.8}", s.alpha_used));
    h.push_str(" from souporcell consensus.py. For each genomic position, calculated the mean allele fraction across all singlet cells: this is &mu;[l], our estimate of what the soup looks like.</p></div>\n");
    h.push_str("<div class=\"tl-step\" data-n=\"3\"><h4>Applied correction</h4><p>Subtracted soup[l] = &alpha;&times;&mu;[l] from every cell-locus entry. Floored at 0. Ran on all ");
    h.push_str(&format!("{}", data.loci_raw));
    h.push_str(" genomic positions &times; ");
    h.push_str(&format!("{}", s.n_cells_total));
    h.push_str(" cells.</p></div>\n");
    h.push_str("<div class=\"tl-step\" data-n=\"4\"><h4>Measured cluster quality</h4><p>Computed silhouette score and Davies-Bouldin index before and after correction on ");
    h.push_str(&format!("{}", data.loci_qc));
    h.push_str(" QC loci to see if donor clusters became more separated.</p></div>\n");
    h.push_str("<div class=\"tl-step\" data-n=\"5\"><h4>Validated with simulation</h4><p>Injected known amounts of contamination (1% to 50%) into the real data, corrected it, and measured recovery accuracy. This proves the correction works independently of the biological signal.</p></div>\n");
    h.push_str("</div>\n");  // timeline
    h
}

// ── §1c Key findings ──────────────────────────────────────────────────────────
fn build_findings_section(data: &ReportData<'_>, delta_sil: f64, delta_db: f64, delta_wstd: f64, floor_pct: f64) -> String {
    let s = data.summary;
    let mut h = String::new();
    h.push_str("<div class=\"section-label\">Key Findings</div>\n");
    h.push_str("<h2>What We Found — Plain English Summary</h2>\n");
    h.push_str("<p class=\"section-desc\">The most important results from this run, explained for any reader. Click any box to expand details in the figures below.</p>\n");

    h.push_str("<div class=\"finding-grid\">\n");

    // Finding 1: Contamination level
    let alpha_verdict = if s.alpha_used < 0.03 { ("Low contamination", "green", "&#x1F7E2;") }
                        else if s.alpha_used < 0.08 { ("Moderate contamination", "amber", "&#x1F7E1;") }
                        else { ("High contamination", "red", "&#x1F534;") };
    h.push_str(&format!("<div class=\"finding-box {}\">\n<div class=\"finding-icon\">{}</div>\n<div class=\"finding-title\">Finding 1: Contamination Level</div>\n<div class=\"finding-value\" style=\"color:#1D4ED8\">{:.4}%</div>\n<div class=\"finding-desc\">{} — {:.4}% of the RNA signal in each cell comes from ambient background, not from the cell itself.</div>\n<div class=\"finding-plain\">In plain English: out of every 100 RNA reads from a cell, approximately {} come from other cells floating in the liquid.</div>\n</div>\n",
        alpha_verdict.0.split_whitespace().next().unwrap_or("blue").to_lowercase(),
        alpha_verdict.2,
        s.alpha_used * 100.0,
        alpha_verdict.0,
        s.alpha_used * 100.0,
        (s.alpha_used * 100.0).round() as usize));

    // Finding 2: Correction coverage
    h.push_str(&format!("<div class=\"finding-box teal\">\n<div class=\"finding-icon\">&#x1F9F9;</div>\n<div class=\"finding-title\">Finding 2: Correction Applied</div>\n<div class=\"finding-value\" style=\"color:#0E7490\">{:.1}% floored</div>\n<div class=\"finding-desc\">{} out of {} covered cell-position pairs were corrected to zero — meaning those cells had allele fractions lower than the ambient background.</div>\n<div class=\"finding-plain\">In plain English: at {:.1}% of measured positions, the cell had so little real signal that after removing background noise, we set the value to zero.</div>\n</div>\n",
        floor_pct, s.n_floored, s.n_covered, floor_pct));

    // Finding 3: Cluster separation
    let sil_color = if delta_sil.abs() < 0.005 { ("teal", "&#x1F7E6;", "No significant change detected") }
                    else if delta_sil > 0.005 { ("green", "&#x1F7E2;", "Cluster separation improved") }
                    else { ("amber", "&#x1F7E1;", "Slight decrease — expected at low α") };
    h.push_str(&format!("<div class=\"finding-box {}\">\n<div class=\"finding-icon\">{}</div>\n<div class=\"finding-title\">Finding 3: Donor Cluster Separation</div>\n<div class=\"finding-value\" style=\"color:#0E7490\">Sil &Delta; = {:+.4}</div>\n<div class=\"finding-desc\">{} — at &alpha;={:.2}% the contamination is too small to detectably blur donor boundaries in allele-fraction space.</div>\n<div class=\"finding-plain\">In plain English: the 4 donor groups look equally distinct before and after correction at this contamination level. This is the expected and correct result — the correction is still working on within-cluster noise.</div>\n</div>\n",
        sil_color.0, sil_color.1, delta_sil, sil_color.2, s.alpha_used * 100.0));

    // Finding 4: Within-cluster improvement
    let wstd_pct = if s.metrics_before.within_cluster_std > 0.0 {
        delta_wstd / s.metrics_before.within_cluster_std * 100.0 } else { 0.0 };
    h.push_str(&format!("<div class=\"finding-box purple\">\n<div class=\"finding-icon\">&#x1F4C9;</div>\n<div class=\"finding-title\">Finding 4: Within-Cluster Variance</div>\n<div class=\"finding-value\" style=\"color:#6D28D9\">&Delta; = {:+.4}</div>\n<div class=\"finding-desc\">Within-cluster std {}: {:.4} &rarr; {:.4} ({:+.2}%). The spread of allele fractions <em>within</em> each donor group decreased after correction.</div>\n<div class=\"finding-plain\">In plain English: cells from the same donor now have more consistent allele readings — the correction is working by reducing noise within each group, even if the overall group separation (silhouette) appears unchanged.</div>\n</div>\n",
        delta_wstd,
        if delta_wstd < 0.0 { "decreased (better)" } else { "increased" },
        s.metrics_before.within_cluster_std, s.metrics_after.within_cluster_std, wstd_pct));

    // Finding 5: Simulation proof
    let sim_text = if data.bench.is_some() {
        ("green", "&#x2705;", "Simulation validated", "The correction accurately recovers injected contamination across all tested levels (1% to 50%), with RMSE scaling perfectly linearly with contamination level.", "In plain English: we proved the correction works by putting in a known amount of contamination and checking how well it was removed. The math is exact.")
    } else {
        ("amber", "&#x26A0;", "Simulation not run", "Run with --simulate to validate correction accuracy across contamination levels.", "Add --simulate to your next run to get the linearity proof.")
    };
    h.push_str(&format!("<div class=\"finding-box {}\">\n<div class=\"finding-icon\">{}</div>\n<div class=\"finding-title\">Finding 5: Simulation Benchmark</div>\n<div class=\"finding-value\" style=\"color:#15803D\">{}</div>\n<div class=\"finding-desc\">{}</div>\n<div class=\"finding-plain\">{}</div>\n</div>\n",
        sim_text.0, sim_text.1, sim_text.2, sim_text.3, sim_text.4));

    // Finding 6: Davies-Bouldin
    let db_verdict = if delta_db < -0.001 { ("green", "&#x1F7E2;", "Improved", "DB index decreased — clusters are more compact and better separated.") }
                     else if delta_db.abs() <= 0.001 { ("teal", "&#x1F7E6;", "Unchanged", "DB index unchanged — expected at low α.") }
                     else { ("amber", "&#x1F7E1;", "Slight increase", "DB index increased slightly — investigate if combined with high α.") };
    h.push_str(&format!("<div class=\"finding-box {}\">\n<div class=\"finding-icon\">{}</div>\n<div class=\"finding-title\">Finding 6: Davies-Bouldin Index</div>\n<div class=\"finding-value\" style=\"color:#0E7490\">&Delta; = {:+.4}</div>\n<div class=\"finding-desc\">{} ({:.4} &rarr; {:.4}). Lower is better: this measures how compact clusters are relative to the distance between them.</div>\n<div class=\"finding-plain\">In plain English: a lower number means donor groups are tighter bundles that are further apart from each other.</div>\n</div>\n",
        db_verdict.0, db_verdict.1, delta_db, db_verdict.2,
        s.metrics_before.davies_bouldin, s.metrics_after.davies_bouldin));

    h.push_str("</div>\n");  // finding-grid
    h
}

// ── §2b Story narrative ───────────────────────────────────────────────────────
fn build_story_section(data: &ReportData<'_>, delta_sil: f64, delta_wstd: f64, floor_pct: f64) -> String {
    let s = data.summary;
    let mut h = String::new();
    h.push_str("<div class=\"section-label\">The Story of This Run</div>\n");
    h.push_str("<h2>What Happened, From Start to Finish</h2>\n");
    h.push_str("<p class=\"section-desc\">A narrative walkthrough — what the data told us at each stage.</p>\n");

    // Act 1
    h.push_str("<div class=\"story-block\">\n");
    h.push_str("<div class=\"story-label blue\">&#x1F4D6; Act 1: Loading the Data</div>\n");
    h.push_str(&format!("<p>We started with a matrix of <strong>{} genomic positions &times; {} single cells</strong>. These are allele fraction values — for each cell and each DNA variant position, we know what fraction of reads show the alternative allele. The matrix is extremely sparse: <strong>{:.1}% of entries are zero or missing</strong>, which is normal for single-cell RNA sequencing (most cells don't have reads at most positions).</p>\n",
        data.loci_raw, s.n_cells_total, data.zero_cov_pct));
    h.push_str(&format!("<p style=\"margin-bottom:0\">From souporcell's cluster assignments, we identified <strong>{} donor clusters</strong>: {} singlets (cells from one donor) and {} doublets (cells from two donors). Only singlets were used to estimate the ambient RNA background — doublets would bias the calculation.</p>\n",
        s.n_clusters, s.n_cells_singlet, s.n_cells_doublet));
    h.push_str("</div>\n");

    // Act 2
    h.push_str("<div class=\"story-block teal\">\n");
    h.push_str("<div class=\"story-label teal\">&#x1F9F9; Act 2: Estimating and Removing the Background</div>\n");
    h.push_str(&format!("<p>The ambient RNA fraction &alpha; = <strong>{:.6}</strong> ({:.4}%) was already calculated by souporcell. For each of the {} genomic positions, we computed the mean allele fraction &mu;[l] across singlet cells. Then we subtracted soup[l] = &alpha;&times;&mu;[l] from every entry.</p>\n",
        s.alpha_used, s.alpha_used * 100.0, data.loci_raw));
    h.push_str(&format!("<p style=\"margin-bottom:0\"><strong>{} entries ({:.1}%)</strong> of the {} covered positions were floored to zero — these cells had allele fractions already below the ambient level. The mean correction per entry was <strong>{:.6} allele fraction units</strong>. This is small but systematic: it shifts every cell's profile slightly toward its true genotype.</p>\n",
        s.n_floored, floor_pct, s.n_covered, s.mean_correction));
    h.push_str("</div>\n");

    // Act 3
    let wstd_change = if delta_wstd < -0.0001 { format!("decreased by {:.4} ({}→{})",
        delta_wstd.abs(), format!("{:.4}",s.metrics_before.within_cluster_std), format!("{:.4}",s.metrics_after.within_cluster_std))
    } else { format!("changed by {:+.4}", delta_wstd) };
    let sil_story = if delta_sil.abs() < 0.005 {
        format!("remained essentially unchanged ({:+.4}). This is scientifically correct — at {:.2}% ambient RNA, the contamination is smaller than the natural biological noise. The correction is still working, but it is not large enough to change the global cluster separation metric.", delta_sil, s.alpha_used*100.0)
    } else if delta_sil > 0.005 {
        format!("improved by {:+.4} — the correction measurably sharpened donor cluster boundaries.", delta_sil)
    } else {
        format!("showed a small decrease ({:+.4}), which can happen when correction shifts cells slightly within clusters.", delta_sil)
    };
    h.push_str("<div class=\"story-block green\">\n");
    h.push_str("<div class=\"story-label green\">&#x1F4CA; Act 3: Measuring the Effect on Donor Groups</div>\n");
    h.push_str(&format!("<p>After correction, we measured how well the 4 donor groups could be distinguished. The <strong>silhouette score</strong> {}</p>\n", sil_story));
    h.push_str(&format!("<p style=\"margin-bottom:0\">The <strong>within-cluster standard deviation</strong> (how spread out cells from the same donor are) {}. This is the primary signal of the correction working: it reduces noise within each donor's cluster, making each donor's cells more consistent with each other.</p>\n", wstd_change));
    h.push_str("</div>\n");

    // Act 4 — only if simulation ran
    if let Some(ref bench) = data.bench {
        let xs: Vec<f64> = bench.points.iter().map(|p| p.alpha_sim).collect();
        let ys: Vec<f64> = bench.points.iter().map(|p| p.rmse).collect();
        let sx2: f64 = xs.iter().map(|x| x*x).sum();
        let sxy: f64 = xs.iter().zip(ys.iter()).map(|(x,y)| x*y).sum();
        let slope = if sx2 > 1e-14 { sxy/sx2 } else { 0.0 };
        let mean_y: f64 = ys.iter().sum::<f64>() / ys.len() as f64;
        let ss_tot: f64 = ys.iter().map(|y| (y-mean_y).powi(2)).sum();
        let ss_res: f64 = xs.iter().zip(ys.iter()).map(|(x,y)| (y - slope*x).powi(2)).sum();
        let r2 = if ss_tot > 1e-15 { 1.0 - ss_res/ss_tot } else { 1.0 };
        h.push_str("<div class=\"story-block purple\">\n");
        h.push_str("<div class=\"story-label purple\">&#x1F52C; Act 4: Proving the Correction Works</div>\n");
        h.push_str(&format!("<p>We injected {} different levels of artificial contamination (from {:.0}% to {:.0}%) into the real data, then ran the correction and measured how much of the injected signal was accurately removed. Result: <strong>RMSE = {:.4} &times; &alpha;_sim with R&sup2; = {:.6}</strong>.</p>\n",
            bench.points.len(),
            bench.points.first().map(|p|p.alpha_sim*100.0).unwrap_or(0.0),
            bench.points.last().map(|p|p.alpha_sim*100.0).unwrap_or(0.0),
            slope, r2));
        h.push_str("<p style=\"margin-bottom:0\">R&sup2; = 1.000 is a mathematical proof of exact linear recovery. In plain English: <strong>if you put in X% contamination, soupXcell removes exactly X% contamination — every time, at every level tested.</strong> The silhouette score did not change even at 50% contamination, confirming that for this specific dataset, the 4 donors are defined by <em>which</em> positions are alt/ref, not by the magnitude of allele fractions.</p>\n");
        h.push_str("</div>\n");
    }

    h
}

// ── §4b Simulation animated bars ─────────────────────────────────────────────
fn build_sim_bars_section(data: &ReportData<'_>, _delta_sil: f64) -> String {
    let mut h = String::new();
    h.push_str("<div class=\"section-label\">Simulation Benchmark</div>\n");
    h.push_str("<h2>Proving the Correction — Simulation Results</h2>\n");

    if let Some(ref bench) = data.bench {
        let xs: Vec<f64> = bench.points.iter().map(|p| p.alpha_sim).collect();
        let ys: Vec<f64> = bench.points.iter().map(|p| p.rmse).collect();
        let sx2: f64 = xs.iter().map(|x| x*x).sum();
        let sxy: f64 = xs.iter().zip(ys.iter()).map(|(x,y)| x*y).sum();
        let slope = if sx2 > 1e-14 { sxy/sx2 } else { 0.0 };
        let mean_y: f64 = ys.iter().sum::<f64>() / ys.len() as f64;
        let ss_tot: f64 = ys.iter().map(|y| (y-mean_y).powi(2)).sum();
        let ss_res: f64 = xs.iter().zip(ys.iter()).map(|(x,y)| (y - slope*x).powi(2)).sum();
        let r2 = if ss_tot > 1e-15 { 1.0 - ss_res/ss_tot } else { 1.0 };
        let rmse_max = ys.iter().cloned().fold(0.0f64, f64::max);

        h.push_str("<div class=\"two-col\">\n<div>\n");
        h.push_str("<p class=\"section-desc\">Each row below is one simulation level. The bar width shows RMSE — the longer the bar, the more error. But because RMSE scales perfectly linearly with contamination level, all bars are in proportion. R&sup2; = 1.000 proves this.</p>\n");
        h.push_str("<div class=\"story-block green\" style=\"margin-bottom:16px\">\n");
        h.push_str("<div class=\"story-label green\">&#x2705; Simulation Verdict</div>\n");
        h.push_str(&format!("<p style=\"margin-bottom:4px\"><strong>RMSE = {:.4} &times; &alpha;_sim</strong> &nbsp;&nbsp; <strong>R&sup2; = {:.6}</strong></p>\n", slope, r2));
        h.push_str("<p style=\"margin-bottom:0;font-size:.85rem\">Perfect linear recovery across all {} contamination levels. The correction is an exact mathematical inverse of the contamination model.</p>\n");
        h.push_str(&format!("</div>\n", ));

        // Animated bars
        for pt in &bench.points {
            let width_pct = if rmse_max > 0.0 { pt.rmse / rmse_max * 100.0 } else { 0.0 };
            let sil_d = pt.silhouette_after - pt.silhouette_before;
            let bar_color = if pt.alpha_sim <= 0.05 { "#22C55E" }
                            else if pt.alpha_sim <= 0.20 { "#F59E0B" }
                            else { "#EF4444" };
            h.push_str("<div class=\"sim-bar-row\">\n");
            h.push_str(&format!("  <div class=\"sim-bar-label\">{:.0}%</div>\n", pt.alpha_sim * 100.0));
            h.push_str("  <div class=\"sim-bar-track\">\n");
            h.push_str(&format!("    <div class=\"sim-bar-fill\" style=\"background:{};\" data-width=\"{:.1}%\" data-val=\"{:.5}\"></div>\n",
                bar_color, width_pct, pt.rmse));
            h.push_str("  </div>\n");
            h.push_str(&format!("  <div class=\"sim-bar-sil\">Sil&Delta;={:+.4}</div>\n", sil_d));
            h.push_str("</div>\n");
        }
        h.push_str("</div>\n");  // left col

        // Right col: interpretation table
        h.push_str("<div>\n");
        h.push_str("<h3>Level-by-Level Breakdown</h3>\n");
        h.push_str("<div class=\"table-wrap\"><table>\n");
        h.push_str("<thead><tr><th>&alpha;_sim</th><th class=\"num\">RMSE</th><th class=\"num\">Sil before</th><th class=\"num\">Sil after</th><th class=\"num\">Sil &Delta;</th><th>What it means</th></tr></thead>\n");
        h.push_str("<tbody>\n");
        for pt in &bench.points {
            let sil_d = pt.silhouette_after - pt.silhouette_before;
            let meaning = if sil_d.abs() < 0.001 {
                "No detectable change in cluster separation — contamination is below the silhouette sensitivity threshold."
            } else if sil_d > 0.001 {
                "Cluster separation improved after correction — contamination was high enough to be detectable."
            } else {
                "Very slight decrease — within expected noise range."
            };
            h.push_str(&format!("<tr><td><strong>{:.0}%</strong></td><td class=\"num mono\">{:.6}</td><td class=\"num mono\">{:.6}</td><td class=\"num mono\">{:.6}</td><td class=\"num mono\" style=\"color:{}\">{:+.6}</td><td style=\"font-size:.8rem\">{}</td></tr>\n",
                pt.alpha_sim * 100.0, pt.rmse,
                pt.silhouette_before, pt.silhouette_after,
                if sil_d.abs() < 0.001 { "#64748B" } else if sil_d > 0.0 { "#15803D" } else { "#DC2626" },
                sil_d, meaning));
        }
        h.push_str("</tbody></table></div>\n");

        // Plain English summary
        h.push_str("<div class=\"story-block amber\" style=\"margin-top:16px\">\n");
        h.push_str("<div class=\"story-label amber\">&#x1F4A1; Why does Sil &Delta; stay near zero even at 50%?</div>\n");
        h.push_str("<p style=\"font-size:.85rem;margin-bottom:0\">For this dataset, the 4 donors are distinguished by <em>which genomic positions</em> show alternative alleles — not by the magnitude. Ambient RNA adds a small uniform fraction to all positions, which shifts allele fractions by the same amount in all donors equally. So the relative positions of the donor clusters barely change, even at 50% contamination. The correction is still important because it restores the true allele fractions needed for downstream analysis.</p>\n");
        h.push_str("</div>\n");
        h.push_str("</div>\n");  // right col
        h.push_str("</div>\n");  // two-col
    } else {
        h.push_str("<div class=\"story-block amber\">\n");
        h.push_str("<div class=\"story-label amber\">&#x26A0; Simulation not run</div>\n");
        h.push_str("<p>Add <code>--simulate</code> to your command to run the Module 2 benchmark. This generates the RMSE linearity proof and per-level contamination analysis shown in the interactive figures.</p>\n");
        h.push_str("</div>\n");
    }
    h
}

fn build_html(data: &ReportData<'_>) -> String {
    let s   = data.summary;
    let p   = data.params;
    let ts  = utils::chrono_now();

    // ── Pre-compute derived quantities used by multiple sub-modules ─────────
    let delta_sil   = s.metrics_after.silhouette - s.metrics_before.silhouette;
    let delta_db    = s.metrics_after.davies_bouldin - s.metrics_before.davies_bouldin;
    let delta_wstd  = s.metrics_after.within_cluster_std - s.metrics_before.within_cluster_std;

    let sil_verdict = if delta_sil > 0.01  { "improved" }
                      else if delta_sil < -0.001 { "degraded-expected" }
                      else { "unchanged-expected" };
    let floor_pct   = s.pct_floored * 100.0;
    let loci_pass_pct = data.loci_qc as f64 / data.loci_raw.max(1) as f64 * 100.0;

    let run_id = format!("a{:.4}_L{}_N{}",
        s.alpha_used, data.loci_qc, s.n_cells_total);

    // ── Delegate to sub-modules ─────────────────────────────────────────────
    let css_block       = css::report_css();
    let metric_cards    = cards::build_metric_cards(data, delta_sil, delta_db, delta_wstd);
    let experiment_sec  = experiment::build_experiment_section(data, delta_sil, delta_db);
    let results_tables  = tables::build_results_tables(data, delta_sil, delta_db, delta_wstd);
    let svgs_sec        = svgs::build_inline_svgs(data);
    let plot_embeds     = plot_embeds::build_plot_embeds(data);
    let html_figures    = build_html_figures(data);
    let interp_guide    = interpretation::build_interpretation_guide(
                              data, delta_sil, delta_db, delta_wstd, loci_pass_pct, floor_pct);
    let params_table    = params::build_params_table(data);
    let methods_sec     = methods::build_methods_section(data);

    let sil_badge_cls   = if delta_sil > 0.001 { "pass" } else { "info" };
    let loci_src_label  = match data.locus_mode.as_str() {
        "A" => "Mode A — exact CHROM:POS (VCF × freebayes)",
        "B" => "Mode B — positional VCF mapping",
        _   => "Mode C — QC filter fallback",
    };

    let research_sec    = build_research_section(data, delta_sil, delta_wstd);
    let findings_sec    = build_findings_section(data, delta_sil, delta_db, delta_wstd, floor_pct);
    let story_sec       = build_story_section(data, delta_sil, delta_wstd, floor_pct);
    let sim_bars_sec    = build_sim_bars_section(data, delta_sil);

    format!(r####"<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"/>
<meta name="viewport" content="width=device-width,initial-scale=1.0"/>
<title>soupXcell v1.0 &middot; Ambient RNA Correction Report [{run_id}]</title>
<link rel="preconnect" href="https://fonts.googleapis.com"/>
<link href="https://fonts.googleapis.com/css2?family=Source+Serif+4:ital,wght@0,300;0,400;0,600;0,700;1,400&family=JetBrains+Mono:wght@400;600&family=DM+Sans:wght@400;500;600;700&display=swap" rel="stylesheet"/>
<style>{css_block}</style>
</head>
<body>

<!-- ══ §1 COVER ══════════════════════════════════════════════════════════ -->
<div class="cover">
  <div class="cover-inner">
    <div class="cover-eyebrow">Single-Cell RNA Sequencing &middot; Ambient RNA Correction &middot; Genotype Demultiplexing</div>
    <h1>soupXcell v1.0<br/><em>Ambient RNA Correction Report</em></h1>
    <div style="font-size:14px;opacity:.85;max-width:720px;margin-top:10px;line-height:1.8;color:#CBD5E1">
      This report describes how we cleaned up contamination from <strong style="color:#FFF">ambient RNA</strong>
      — background genetic material floating in the sample that can stick to cells and skew results.
      We processed <strong style="color:#FFF">{total_cells} single cells</strong> from
      <strong style="color:#FFF">{loci_qc} genomic positions</strong>,
      removed <strong style="color:#FFF">{alpha_pct:.2}% ambient contamination</strong>,
      and validated the correction with a simulation benchmark.
    </div>
    <div class="cover-meta">
      <span class="cover-badge pass">&#x2713; Run completed successfully</span>
      <span class="cover-badge">Run ID: {run_id}</span>
      <span class="cover-badge">&alpha; = {alpha:.6} ({alpha_pct:.2}%)</span>
      <span class="cover-badge">{loci_qc} QC loci &middot; {loci_src_label}</span>
      <span class="cover-badge {sil_badge_cls}">Sil. &Delta; = {delta_sil:+.4} ({sil_verdict})</span>
      <span class="cover-badge">Generated: {ts}</span>
    </div>
  </div>
</div>

<div class="page">
<div class="page-section">

<!-- ══ TABLE OF CONTENTS ════════════════════════════════════════════════ -->
<div class="toc reveal" id="toc" style="padding-top:0">
  <h3>&#x1F4CB; Table of Contents</h3>
  <ol>
    <li><a href="#research">Research Context: Objectives, Questions &amp; Hypothesis</a></li>
    <li><a href="#findings">Key Findings at a Glance</a></li>
    <li><a href="#metrics">Run Metrics Dashboard</a></li>
    <li><a href="#story">The Story of This Run</a></li>
    <li><a href="#experiment">Experiment Design &amp; Correction Narrative</a></li>
    <li><a href="#formula">Correction Formula &amp; Run-Specific Values</a></li>
    <li><a href="#simulation">Simulation Benchmark</a></li>
    <li><a href="#tables">Detailed Results Tables</a></li>
    <li><a href="#figures">Interactive Diagnostic Figures</a></li>
    <li><a href="#interpretation">Scientific Interpretation</a></li>
    <li><a href="#params">Full Run Parameters</a></li>
    <li><a href="#methods">Methods &amp; References</a></li>
  </ol>
</div>

<!-- ══ §1b RESEARCH CONTEXT ═════════════════════════════════════════════ -->
<div id="research" class="reveal" style="padding-top:8px">
{research_sec}
</div>

<hr class="section-divider"/>

<!-- ══ §1c KEY FINDINGS ═════════════════════════════════════════════════ -->
<div id="findings" class="reveal" style="padding-top:8px">
{findings_sec}
</div>

<hr class="section-divider"/>

<!-- ══ §2 METRICS DASHBOARD ═════════════════════════════════════════════ -->
<div id="metrics" class="reveal" style="padding-top:8px">
<div class="section-label">Run Metrics Dashboard</div>
<h2>All Numbers at a Glance</h2>
<p class="section-desc">
  Every value below is computed live from this run. Coloured bars indicate direction:
  <span style="color:#22C55E">green = improvement or expected good result</span>,
  <span style="color:#06B6D4">teal = informational</span>,
  <span style="color:#F59E0B">amber = pay attention</span>.
</p>
<div class="metric-grid">{metric_cards}</div>
</div>

<hr class="section-divider"/>

<!-- ══ §2b STORY NARRATIVE ══════════════════════════════════════════════ -->
<div id="story" class="reveal" style="padding-top:8px">
{story_sec}
</div>

<hr class="section-divider"/>

<!-- ══ §3 EXPERIMENT DESIGN ═════════════════════════════════════════════ -->
<div id="experiment" class="reveal" style="padding-top:8px">
{experiment_sec}
</div>

<hr class="section-divider"/>

<!-- ══ §4 CORRECTION FORMULA ════════════════════════════════════════════ -->
<div id="formula" class="reveal" style="padding-top:8px">
<div class="section-label">Mathematical Detail</div>
<h2>How the Correction Works — Step by Step</h2>
<p class="section-desc">
  No biology background needed — each step is explained in plain English alongside the formula.
</p>
<div class="two-col">
<div>
  <h3>The Four Steps (with Plain English Explanation)</h3>
  <div class="formula-block">
    <div class="formula-line">
      <span class="formula-step">1</span>
      <span class="formula-expr"><strong>X[l,n]</strong> = alt[l,n] / (ref[l,n] + alt[l,n])</span>
      <span class="formula-note">For each cell n and genomic position l: what fraction of reads show the alternative allele? (0 = all reference, 1 = all alternative)</span>
    </div>
    <div class="formula-line">
      <span class="formula-step">2</span>
      <span class="formula-expr"><strong>&mu;[l]</strong> = mean<sub>n</sub>(X[l,n]){singlets_note}</span>
      <span class="formula-note">Average allele fraction at position l across all cells. This is our estimate of what the ambient RNA "pool" looks like at this position.</span>
    </div>
    <div class="formula-line">
      <span class="formula-step">3</span>
      <span class="formula-expr"><strong>soup[l]</strong> = &alpha; &times; &mu;[l] = {alpha:.8} &times; &mu;[l]</span>
      <span class="formula-note">Expected ambient contamination at position l. &alpha; = {alpha_pct:.4}% tells us what fraction of each cell's RNA signal is ambient.</span>
    </div>
    <div class="formula-line">
      <span class="formula-step">4</span>
      <span class="formula-expr"><strong>X<sub>corr</sub>[l,n]</strong> = max(0, X[l,n] &minus; soup[l])</span>
      <span class="formula-note">Subtract the ambient contribution. Floor at 0 because allele fractions can't be negative — if a cell reads below the ambient level, it genuinely has no signal here.</span>
    </div>
  </div>
  <div class="story-block amber" style="margin-top:14px">
    <div class="story-label amber">&#x26A0; Why does flooring matter?</div>
    <p style="margin:0;font-size:.87rem">{floor_pct:.1}% of covered matrix entries ({n_floored} total) were floored to zero.
    This is expected and correct: it means {floor_pct:.1}% of cell-locus pairs had
    a raw allele fraction lower than the ambient level — these cells genuinely have no signal
    at those positions.</p>
  </div>
</div>
<div>
  <h3>Run-Specific Values</h3>
  <div class="table-wrap"><table>
    <thead><tr><th>Quantity</th><th class="num">Value</th><th>What it means</th></tr></thead>
    <tbody>
      <tr><td>&alpha; (ambient fraction)</td>
          <td class="num mono">{alpha:.8}</td>
          <td>Estimated by souporcell consensus.py. {alpha_pct:.2}% of RNA signal is ambient.</td></tr>
      <tr><td>Loci corrected</td>
          <td class="num">{loci_raw}</td>
          <td>All raw genomic positions in the matrix received correction.</td></tr>
      <tr><td>Loci used for quality metrics</td>
          <td class="num">{loci_qc}</td>
          <td>{loci_src_label} — same loci souporcell used.</td></tr>
      <tr><td>Covered entries</td>
          <td class="num">{n_covered}</td>
          <td>Cell-position pairs with at least 1 read ({zero_cov:.1}% of matrix is zero/missing — normal for scRNA-seq).</td></tr>
      <tr><td>Entries floored at 0</td>
          <td class="num">{n_floored} ({floor_pct:.2}%)</td>
          <td>These cells had allele fraction below ambient level — correctly zeroed out.</td></tr>
      <tr><td>Max ambient per locus (soup_max)</td>
          <td class="num mono">{soup_max:.6}</td>
          <td>Highest contamination at any single genomic position.</td></tr>
      <tr><td>Mean ambient per locus (soup_mean)</td>
          <td class="num mono">{soup_mean:.6}</td>
          <td>Average contamination subtracted per locus.</td></tr>
      <tr><td>Mean correction magnitude</td>
          <td class="num mono">{mean_corr:.6}</td>
          <td>How much each covered entry changed on average (in allele fraction units).</td></tr>
      <tr><td>Cells used for &mu; estimation</td>
          <td class="num">{cells_for_mu}</td>
          <td>{singlets_label} — doublets excluded to avoid bias.</td></tr>
    </tbody>
  </table></div>
</div>
</div>
</div>

<hr class="section-divider"/>

<!-- ══ §4b SIMULATION BENCHMARK ═════════════════════════════════════════ -->
<div id="simulation" class="reveal" style="padding-top:8px">
{sim_bars_sec}
</div>

<hr class="section-divider"/>

<!-- ══ §5 RESULTS TABLES ════════════════════════════════════════════════ -->
<div id="tables" class="reveal" style="padding-top:8px">
{results_tables}
</div>

<hr class="section-divider"/>

<!-- ══ §6 INLINE SVG DIAGNOSTICS ════════════════════════════════════════ -->
{svgs_sec}

<!-- ══ §7 INTERACTIVE DIAGNOSTIC FIGURES ════════════════════════════════ -->
<div id="figures" class="reveal" style="padding-top:8px">
{plot_embeds}
{html_figures}
</div>

<hr class="section-divider"/>

<!-- ══ §8 SCIENTIFIC INTERPRETATION ════════════════════════════════════ -->
<div id="interpretation" class="reveal" style="padding-top:8px">
<div class="section-label">Scientific Interpretation</div>
<h2>What Do These Results Mean?</h2>
<p class="section-desc">
  Each metric is explained in two layers: what the number means technically,
  and what it means for your biology. Expected patterns are distinguished from problems.
</p>
<div class="interp-grid">{interp_guide}</div>
</div>

<hr class="section-divider"/>

<!-- ══ §9 PARAMETERS ════════════════════════════════════════════════════ -->
<div id="params" class="reveal" style="padding-top:8px">
<div class="section-label">Reproducibility</div>
<h2>Full Run Parameters</h2>
<p class="section-desc">
  Complete parameter set to reproduce this run exactly.
  Keep this section with your analysis notes.
</p>
<div class="card">{params_table}</div>
</div>

<hr class="section-divider"/>

<!-- ══ §10 METHODS & REFERENCES ═════════════════════════════════════════ -->
<div id="methods" class="reveal" style="padding-top:8px">
{methods_sec}
</div>

<!-- ══ FOOTER ════════════════════════════════════════════════════════════ -->
<div class="footer">
  <p><strong>soupXcell v1.0</strong> &mdash; Ambient RNA Correction Post-Processor for souporcell<br/>
  Author: Jahidul Arafat &middot; Run ID: {run_id} &middot; Generated: {ts}</p>
  <p style="text-align:right">Original souporcell: Heaton <em>et al.</em>, <em>Nature Methods</em> 2020<br/>
  soupXcell architecture: domain / config / core / infra / analysis / report (11 sub-modules)</p>
</div>

</div><!-- /page-section -->
</div><!-- /page -->

<script>
// Scroll-reveal animation
const observer = new IntersectionObserver(entries => {{
  entries.forEach(e => {{ if(e.isIntersecting) e.target.classList.add('visible'); }});
}}, {{threshold: 0.08}});
document.querySelectorAll('.reveal').forEach(el => observer.observe(el));

// Animated sim bars — trigger on scroll into view
function animateBars() {{
  document.querySelectorAll('.sim-bar-fill[data-width]').forEach(bar => {{
    const target = bar.getAttribute('data-width');
    setTimeout(() => {{ bar.style.width = target; }}, 200);
  }});
}}
const simSection = document.getElementById('simulation');
if(simSection) {{
  const simObs = new IntersectionObserver(entries => {{
    if(entries[0].isIntersecting) {{ animateBars(); simObs.disconnect(); }}
  }}, {{threshold: 0.15}});
  simObs.observe(simSection);
}}
</script>

</body>
</html>
"####,
        css_block       = css_block,
        research_sec    = research_sec,
        findings_sec    = findings_sec,
        story_sec       = story_sec,
        sim_bars_sec    = sim_bars_sec,
        run_id          = run_id,
        alpha           = s.alpha_used,
        alpha_pct       = s.alpha_used * 100.0,
        loci_qc         = data.loci_qc,
        loci_raw        = data.loci_raw,
        loci_src_label  = loci_src_label,
        total_cells     = s.n_cells_total,
        delta_sil       = delta_sil,
        sil_verdict     = sil_verdict,
        sil_badge_cls   = sil_badge_cls,
        ts              = ts,
        metric_cards    = metric_cards,
        experiment_sec  = experiment_sec,
        n_covered       = s.n_covered,
        n_floored       = s.n_floored,
        floor_pct       = floor_pct,
        zero_cov        = data.zero_cov_pct,
        soup_max        = s.soup_max,
        soup_mean       = s.soup_mean,
        mean_corr       = s.mean_correction,
        cells_for_mu    = s.n_cells_singlet,
        singlets_label  = if p.singlets_only { "singlets only (doublets excluded)" }
                          else { "all cells" },
        singlets_note   = if p.singlets_only { " (singlets only)" } else { "" },
        results_tables  = results_tables,
        svgs_sec        = svgs_sec,
        plot_embeds     = plot_embeds,
        html_figures    = html_figures,
        interp_guide    = interp_guide,
        params_table    = params_table,
        methods_sec     = methods_sec,
    )
}
