// ============================================================================
// core/corrector.rs — Module 1: Ambient RNA subtraction
// ============================================================================
// Implements the supervisor's three-step correction formula:
//   Step 1: X[l,n] = alt[l,n] / (ref[l,n] + alt[l,n])   (allele fraction)
//   Step 2: soup[l] = α × mean(X[l,:])                    (per-locus soup)
//   Step 3: X_corr[l,n] = max(0, X[l,n] − soup[l])       (subtract + floor)
// ============================================================================

use crate::domain::types::{RawMatrices, CorrectionResult};
use crate::domain::math::{
    compute_allele_fractions, compute_soup_vector,
    subtract_ambient, mean_correction_magnitude,
};
use crate::infra::logger::Logger;
use std::time::Instant;

/// Run the full ambient correction pipeline (Module 1).
/// If `cell_filter` is Some, only those column indices are included in the
/// mean computation (e.g. singlets only).
pub fn run_correction(
    raw:         &RawMatrices,
    alpha:       f64,
    cell_filter: Option<&[usize]>,  // column indices to include in mean
    log:         &Logger,
) -> CorrectionResult {
    let t = Instant::now();

    // ── Step 1: allele fractions ──────────────────────────────────────────────
    log.computing(&format!("allele fraction matrix  {} loci × {} cells", raw.n_loci, raw.n_cells));
    let xmat = compute_allele_fractions(raw);
    // Count non-NaN entries for coverage stats
    let n_covered_step1 = xmat.values.iter().filter(|v| !v.is_nan()).count();
    let total_entries   = raw.n_loci * raw.n_cells;
    log.info("step 1/3 complete", "allele fraction matrix built");
    log.info("matrix dimensions",   &format!("{} loci × {} cells", raw.n_loci, raw.n_cells));
    log.info("covered entries",     &format!("{} / {} ({:.1}% non-NaN)",
        n_covered_step1, total_entries,
        n_covered_step1 as f64 / total_entries as f64 * 100.0));
    log.info("zero-coverage entries", &format!("{} ({:.1}% sparse)",
        total_entries - n_covered_step1,
        (total_entries - n_covered_step1) as f64 / total_entries as f64 * 100.0));

    // If singlets_only: zero out doublet columns for mean calculation
    // We compute soup only from cell_filter columns then apply to all cells.
    let soup = if let Some(cols) = cell_filter {
        log.info("cells used for soup estimation", &cols.len().to_string());
        // Build a view with only the filtered columns for mean computation
        let n = raw.n_cells;
        let l = raw.n_loci;
        // Compute locus means using only filtered cells
        let soup_alpha = alpha;
        let mut locus_mean = Vec::with_capacity(l);
        for li in 0..l {
            let (sum, count) = cols.iter().fold((0.0_f64, 0usize), |(s, c), &ni| {
                if ni < n {
                    let v = xmat.values[[li, ni]];
                    if !v.is_nan() { (s + v, c + 1) } else { (s, c) }
                } else { (s, c) }
            });
            locus_mean.push(if count > 0 { sum / count as f64 } else { 0.0 });
        }
        let soup_contribution: Vec<f64> = locus_mean.iter().map(|&m| alpha * m).collect();
        crate::domain::types::SoupVector {
            locus_mean, soup_contribution, alpha: soup_alpha, n_loci: l
        }
    } else {
        log.info("step", "2/3 — computing soup vector (all cells)");
        compute_soup_vector(&xmat, alpha)
    };

    log.info("alpha",         &format!("{:.8}", alpha));
    log.info("soup max",      &format!("{:.6}", soup.soup_contribution.iter().cloned().fold(f64::NEG_INFINITY, f64::max)));
    log.info("soup mean",     &format!("{:.6}", soup.soup_contribution.iter().sum::<f64>() / soup.n_loci as f64));

    // ── Step 3: subtract and floor ────────────────────────────────────────────
    log.computing("subtracting soup vector and flooring at 0");
    let (x_corrected, n_floored, n_covered) = subtract_ambient(&xmat, &soup);

    let pct_floored = if n_covered > 0 { n_floored as f64 / n_covered as f64 } else { 0.0 };
    let mean_corr   = mean_correction_magnitude(&xmat.values, &x_corrected);

    log.info("entries floored at 0", &format!("{} ({:.2}% of covered)", n_floored, pct_floored * 100.0));
    log.info("mean |correction|",    &format!("{:.6}", mean_corr));
    log.ok("Module 1 — Ambient subtraction complete", Some(t));

    CorrectionResult {
        x_raw:           xmat.values,
        x_corrected,
        soup_vector:     soup,
        alpha_used:      alpha,
        n_floored,
        pct_floored,
        mean_correction: mean_corr,
        n_covered,
    }
}
