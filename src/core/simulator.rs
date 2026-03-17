// ============================================================================
// core/simulator.rs — Module 2: Synthetic contamination injection & benchmark
// ============================================================================
// Addresses supervisor: "1.26% is small — we might need to simulate more
// to see a difference."
//
// For each α_sim level:
//   1. Inject:   X_contaminated = clamp(X_raw + α_sim × μ, 0, 1)
//   2. Recover:  X_recovered    = max(0, X_contaminated − α_sim × μ)
//   3. Measure:  RMSE(X_raw, X_recovered), silhouette before/after, floor %
// ============================================================================

use crate::domain::types::{BenchmarkResult, BenchmarkPoint};
use crate::domain::math::{inject_contamination, subtract_ambient,
                           silhouette_score, floor_fraction, rmse};
use crate::domain::types::SoupVector;
use crate::infra::logger::Logger;
use ndarray::Array2;
use std::time::Instant;

pub fn run_simulation(
    x_raw:      &Array2<f64>,
    locus_mean: &[f64],
    labels:     &[usize],        // cluster label per cell (column)
    n_clusters:  usize,
    sim_levels:  &[f64],
    sim_trials:  usize,
    true_alpha:  f64,
    log:         &Logger,
) -> BenchmarkResult {
    let t = Instant::now();
    log.info("simulation levels",  &format!("{:?}", sim_levels));
    log.info("trials per level",   &sim_trials.to_string());

    let mut points: Vec<BenchmarkPoint> = Vec::new();

    for (i, &alpha_sim) in sim_levels.iter().enumerate() {
        log.progress(&format!("α={:.4}", alpha_sim), i, sim_levels.len());
        log.computing(&format!("simulation level {}/{}: α_sim={:.4}  ({} trials)",
            i+1, sim_levels.len(), alpha_sim, sim_trials));

        // Average over sim_trials (results should be deterministic since
        // inject_contamination is pure, but we keep trials for extension)
        let mut rmse_sum = 0.0_f64;
        let mut sil_before_sum = 0.0_f64;
        let mut sil_after_sum  = 0.0_f64;
        let mut floor_before_sum = 0.0_f64;
        let mut floor_after_sum  = 0.0_f64;
        let mut mean_corr_sum = 0.0_f64;

        for _trial in 0..sim_trials {
            // Inject synthetic contamination
            let x_contaminated = inject_contamination(x_raw, locus_mean, alpha_sim);

            // Build a SoupVector for recovery using the TRUE alpha_sim
            let soup_for_recovery = SoupVector {
                locus_mean:        locus_mean.to_vec(),
                soup_contribution: locus_mean.iter().map(|&m| alpha_sim * m).collect(),
                alpha:             alpha_sim,
                n_loci:            locus_mean.len(),
            };

            // Build a minimal AlleleFractionMatrix wrapper
            let xmat_contaminated = crate::domain::types::AlleleFractionMatrix {
                values:  x_contaminated.clone(),
                n_loci:  x_contaminated.nrows(),
                n_cells: x_contaminated.ncols(),
            };

            // Recover using Module 1 subtraction
            let (x_recovered, n_floored_after, n_covered_after) =
                subtract_ambient(&xmat_contaminated, &soup_for_recovery);

            // RMSE: how well did we recover X_raw?
            rmse_sum += rmse(x_raw, &x_recovered);

            // Silhouette before (contaminated) and after (recovered)
            let (sil_before, _) = silhouette_score(&x_contaminated, labels, n_clusters, 150);
            let (sil_after, _)  = silhouette_score(&x_recovered,    labels, n_clusters, 150);
            sil_before_sum += sil_before;
            sil_after_sum  += sil_after;

            // Floor fractions
            floor_before_sum += floor_fraction(x_raw, &x_contaminated);
            floor_after_sum  += if n_covered_after > 0 {
                n_floored_after as f64 / n_covered_after as f64
            } else { 0.0 };

            // Mean correction magnitude
            let mc = crate::domain::math::mean_correction_magnitude(&x_contaminated, &x_recovered);
            mean_corr_sum += mc;
        }

        let n = sim_trials as f64;
        let rmse_avg   = rmse_sum / n;
        let sil_b_avg  = sil_before_sum / n;
        let sil_a_avg  = sil_after_sum / n;
        // Machine-parseable tab-delimited line (mirrors souporcell CONV/RESTART format)
        eprintln!("SOUPX_SIM\talpha={:.4}\trmse={:.6}\tsil_before={:.4}\tsil_after={:.4}\tdelta={:+.4}\tfloor={:.4}",
            alpha_sim, rmse_avg, sil_b_avg, sil_a_avg, sil_a_avg - sil_b_avg, floor_after_sum / n);
        log.info(
            &format!("  α={:.4} result", alpha_sim),
            &format!("rmse={:.6}  sil_before={:.4}  sil_after={:.4}  Δsil={:+.4}  floor={:.2}%",
                rmse_avg, sil_b_avg, sil_a_avg,
                sil_a_avg - sil_b_avg,
                (floor_after_sum / n) * 100.0)
        );
        points.push(BenchmarkPoint {
            alpha_sim,
            rmse:               rmse_avg,
            silhouette_before:  sil_b_avg,
            silhouette_after:   sil_a_avg,
            floor_pct_before:   floor_before_sum / n,
            floor_pct_after:    floor_after_sum  / n,
            mean_correction:    mean_corr_sum / n,
        });
    }

    log.progress("done", sim_levels.len(), sim_levels.len());
    log.ok(&format!("Module 2 — Simulation benchmark complete ({} levels × {} trials)",
        sim_levels.len(), sim_trials), Some(t));

    BenchmarkResult { points, true_alpha }
}
