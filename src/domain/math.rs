// ============================================================================
// domain/math.rs — Pure mathematical operations for soupxcell
// ============================================================================
// No I/O, no logging, no side effects. All functions are deterministic
// given their inputs.
// ============================================================================

use ndarray::{Array2, s};
use crate::domain::types::{AlleleFractionMatrix, SoupVector, RawMatrices};

// ── Allele fraction computation ───────────────────────────────────────────────

/// Compute X[l,n] = alt[l,n] / (ref[l,n] + alt[l,n]).
/// Returns NaN where total coverage is zero.
pub fn compute_allele_fractions(raw: &RawMatrices) -> AlleleFractionMatrix {
    let (l, n) = (raw.n_loci, raw.n_cells);
    let mut values = Array2::<f64>::zeros((l, n));
    for li in 0..l {
        for ni in 0..n {
            let r = raw.ref_counts[[li, ni]];
            let a = raw.alt_counts[[li, ni]];
            let tot = r + a;
            values[[li, ni]] = if tot > 0.0 { a / tot } else { f64::NAN };
        }
    }
    AlleleFractionMatrix { values, n_loci: l, n_cells: n }
}

// ── Soup vector computation ───────────────────────────────────────────────────

/// Compute the per-locus mean allele fraction (ignoring NaN/uncovered entries)
/// then multiply by alpha to produce the soup contribution vector.
///
/// Supervisor spec:
///   μ[l] = mean(X[l, :])  — mean over all cells that have coverage
///   soup[l] = α × μ[l]
pub fn compute_soup_vector(xmat: &AlleleFractionMatrix, alpha: f64) -> SoupVector {
    let l = xmat.n_loci;
    let mut locus_mean = Vec::with_capacity(l);

    for li in 0..l {
        let row = xmat.values.slice(s![li, ..]);
        let (sum, count) = row.iter().fold((0.0_f64, 0usize), |(s, c), &v| {
            if v.is_nan() { (s, c) } else { (s + v, c + 1) }
        });
        let mean = if count > 0 { sum / count as f64 } else { 0.0 };
        locus_mean.push(mean);
    }

    let soup_contribution: Vec<f64> = locus_mean.iter().map(|&m| alpha * m).collect();

    SoupVector { locus_mean, soup_contribution, alpha, n_loci: l }
}

// ── Ambient subtraction (supervisor's core formula) ──────────────────────────

/// X_corrected[l,n] = max(0.0,  X[l,n] − soup[l])
/// NaN positions stay NaN (zero coverage → no correction possible).
/// Returns (corrected matrix, n_floored, n_covered).
pub fn subtract_ambient(
    xmat:  &AlleleFractionMatrix,
    soup:  &SoupVector,
) -> (Array2<f64>, usize, usize) {
    let (l, n) = (xmat.n_loci, xmat.n_cells);
    let mut corrected = xmat.values.clone();
    let mut n_floored = 0usize;
    let mut n_covered = 0usize;

    for li in 0..l {
        let s = soup.soup_contribution[li];
        for ni in 0..n {
            let v = corrected[[li, ni]];
            if !v.is_nan() {
                n_covered += 1;
                let sub = v - s;
                if sub < 0.0 {
                    corrected[[li, ni]] = 0.0;
                    n_floored += 1;
                } else {
                    corrected[[li, ni]] = sub;
                }
            }
        }
    }
    (corrected, n_floored, n_covered)
}

/// Compute mean absolute correction over covered (non-NaN) entries.
pub fn mean_correction_magnitude(x_raw: &Array2<f64>, x_corr: &Array2<f64>) -> f64 {
    let mut sum = 0.0_f64;
    let mut count = 0usize;
    for (&r, &c) in x_raw.iter().zip(x_corr.iter()) {
        if !r.is_nan() {
            sum += (r - c).abs();
            count += 1;
        }
    }
    if count > 0 { sum / count as f64 } else { 0.0 }
}

// ── RMSE ─────────────────────────────────────────────────────────────────────

/// Root-mean-squared error between two matrices (NaN positions excluded).
pub fn rmse(a: &Array2<f64>, b: &Array2<f64>) -> f64 {
    let mut sum_sq = 0.0_f64;
    let mut count  = 0usize;
    for (&av, &bv) in a.iter().zip(b.iter()) {
        if !av.is_nan() && !bv.is_nan() {
            let d = av - bv;
            sum_sq += d * d;
            count  += 1;
        }
    }
    if count > 0 { (sum_sq / count as f64).sqrt() } else { f64::NAN }
}

// ── Silhouette score ──────────────────────────────────────────────────────────

/// Compute the global silhouette coefficient and per-cluster scores.
/// Uses Euclidean distance on the cell × locus matrix (columns = cells).
/// `labels[n]` = cluster index for cell n.
/// NaN values are treated as 0 for distance computation.
pub fn silhouette_score(
    x:              &Array2<f64>,   // L × N
    labels:         &[usize],       // cluster index per cell, length N
    n_clusters:     usize,
    max_per_cluster: usize,         // 0 = all cells (exact); >0 = subsample for speed
) -> (f64, Vec<f64>) {
    let (l, n) = (x.nrows(), x.ncols());
    assert_eq!(labels.len(), n, "labels length must equal n_cells");

    if n < 2 || n_clusters < 2 {
        return (0.0, vec![0.0; n_clusters]);
    }

    // Subsample: use max_per_cluster cells per cluster (0 = all cells, exact result).
    let effective_max = if max_per_cluster == 0 { usize::MAX } else { max_per_cluster };

    // Step 1: group original cell indices by cluster
    let mut by_cluster: Vec<Vec<usize>> = vec![Vec::new(); n_clusters];
    for (ni, &k) in labels.iter().enumerate() {
        if k < n_clusters { by_cluster[k].push(ni); }
    }

    // Step 2: take an evenly-spaced subsample from each cluster
    let mut sampled_orig:   Vec<usize> = Vec::new();
    let mut sampled_labels: Vec<usize> = Vec::new();
    for (k, members) in by_cluster.iter().enumerate() {
        if members.is_empty() { continue; }
        let take = effective_max.min(members.len());
        // stride through the member list at even intervals
        let denom = take.max(1) as f64;
        for i in 0..take {
            let idx = ((i as f64 / denom) * members.len() as f64) as usize;
            let idx = idx.min(members.len() - 1);
            sampled_orig.push(members[idx]);
            sampled_labels.push(k);
        }
    }

    let n_sub = sampled_orig.len();
    if n_sub < 2 { return (0.0, vec![0.0; n_clusters]); }

    // Step 3: build dense cell vectors for the sampled cells (NaN → 0)
    let cell_vecs: Vec<Vec<f64>> = sampled_orig.iter().map(|&orig_ni| {
        (0..l).map(|li| { let v = x[[li, orig_ni]]; if v.is_nan() { 0.0 } else { v } }).collect()
    }).collect();

    // Step 4: cluster membership lists in subsampled index space (0..n_sub)
    let mut clusters: Vec<Vec<usize>> = vec![Vec::new(); n_clusters];
    for (sub_i, &k) in sampled_labels.iter().enumerate() {
        if k < n_clusters { clusters[k].push(sub_i); }
    }

    let euclidean = |a: &[f64], b: &[f64]| -> f64 {
        a.iter().zip(b.iter()).map(|(x, y)| (x - y).powi(2)).sum::<f64>().sqrt()
    };

    // Step 5: compute per-cell silhouette in subsampled space
    let mut cell_sil = vec![0.0_f64; n_sub];
    for sub_i in 0..n_sub {
        let k = sampled_labels[sub_i];
        if k >= n_clusters { continue; }
        let my_cluster = &clusters[k];

        let a = if my_cluster.len() > 1 {
            let sum: f64 = my_cluster.iter()
                .filter(|&&j| j != sub_i)
                .map(|&j| euclidean(&cell_vecs[sub_i], &cell_vecs[j]))
                .sum();
            sum / (my_cluster.len() - 1) as f64
        } else { 0.0 };

        let b = (0..n_clusters)
            .filter(|&ck| ck != k && !clusters[ck].is_empty())
            .map(|ck| {
                let sum: f64 = clusters[ck].iter()
                    .map(|&j| euclidean(&cell_vecs[sub_i], &cell_vecs[j]))
                    .sum();
                sum / clusters[ck].len() as f64
            })
            .fold(f64::INFINITY, f64::min);

        cell_sil[sub_i] = if !b.is_finite() || (a == 0.0 && b == 0.0) { 0.0 }
                          else { (b - a) / b.max(a) };
    }

    // Step 6: per-cluster and global silhouette
    let mut per_cluster = vec![0.0_f64; n_clusters];
    for k in 0..n_clusters {
        if clusters[k].is_empty() { continue; }
        let sum: f64 = clusters[k].iter().map(|&i| cell_sil[i]).sum();
        per_cluster[k] = sum / clusters[k].len() as f64;
    }
    let global = cell_sil.iter().sum::<f64>() / n_sub as f64;
    (global, per_cluster)
}


// ── Within-cluster variance ───────────────────────────────────────────────────

/// Mean standard deviation of allele fractions within each cluster,
/// averaged across clusters. Lower = tighter clusters.
pub fn mean_within_cluster_std(
    x:          &Array2<f64>,   // L × N
    labels:     &[usize],
    n_clusters: usize,
) -> f64 {
    let (l, _n) = (x.nrows(), x.ncols());
    let mut clusters: Vec<Vec<usize>> = vec![Vec::new(); n_clusters];
    for (ni, &k) in labels.iter().enumerate() {
        if k < n_clusters { clusters[k].push(ni); }
    }

    let mut total_std = 0.0_f64;
    let mut n_non_empty = 0usize;

    for k in 0..n_clusters {
        if clusters[k].len() < 2 { continue; }
        n_non_empty += 1;
        let mut locus_stds = 0.0_f64;
        let mut locus_count = 0usize;

        for li in 0..l {
            let vals: Vec<f64> = clusters[k].iter()
                .map(|&ni| x[[li, ni]])
                .filter(|v| !v.is_nan())
                .collect();
            if vals.len() < 2 { continue; }
            locus_count += 1;
            let mean = vals.iter().sum::<f64>() / vals.len() as f64;
            let variance = vals.iter().map(|v| (v - mean).powi(2)).sum::<f64>()
                         / vals.len() as f64;
            locus_stds += variance.sqrt();
        }
        if locus_count > 0 {
            total_std += locus_stds / locus_count as f64;
        }
    }
    if n_non_empty > 0 { total_std / n_non_empty as f64 } else { 0.0 }
}

// ── Davies-Bouldin index ──────────────────────────────────────────────────────

/// Davies-Bouldin index — lower values indicate better cluster separation.
pub fn davies_bouldin(
    x:          &Array2<f64>,   // L × N
    labels:     &[usize],
    n_clusters: usize,
) -> f64 {
    let (l, _n) = (x.nrows(), x.ncols());
    if n_clusters < 2 { return f64::NAN; }

    // Cluster centroids
    let mut clusters: Vec<Vec<usize>> = vec![Vec::new(); n_clusters];
    for (ni, &k) in labels.iter().enumerate() {
        if k < n_clusters { clusters[k].push(ni); }
    }

    let mut centroids: Vec<Vec<f64>> = vec![vec![0.0; l]; n_clusters];
    for k in 0..n_clusters {
        if clusters[k].is_empty() { continue; }
        for li in 0..l {
            let vals: Vec<f64> = clusters[k].iter()
                .map(|&ni| x[[li, ni]])
                .filter(|v| !v.is_nan())
                .collect();
            if !vals.is_empty() {
                centroids[k][li] = vals.iter().sum::<f64>() / vals.len() as f64;
            }
        }
    }

    let euclidean = |a: &[f64], b: &[f64]| -> f64 {
        a.iter().zip(b).map(|(x,y)| (x-y).powi(2)).sum::<f64>().sqrt()
    };

    // Per-cluster scatter (mean distance of members to centroid)
    let scatter: Vec<f64> = (0..n_clusters).map(|k| {
        if clusters[k].is_empty() { return 0.0; }
        let sum: f64 = clusters[k].iter().map(|&ni| {
            let cell: Vec<f64> = (0..l).map(|li| {
                let v = x[[li, ni]]; if v.is_nan() { 0.0 } else { v }
            }).collect();
            euclidean(&cell, &centroids[k])
        }).sum();
        sum / clusters[k].len() as f64
    }).collect();

    // DB index
    let mut db_sum = 0.0_f64;
    for i in 0..n_clusters {
        if clusters[i].is_empty() { continue; }
        let max_r = (0..n_clusters)
            .filter(|&j| j != i && !clusters[j].is_empty())
            .map(|j| {
                let d = euclidean(&centroids[i], &centroids[j]);
                if d > 0.0 { (scatter[i] + scatter[j]) / d } else { f64::INFINITY }
            })
            .fold(f64::NEG_INFINITY, f64::max);
        if max_r.is_finite() { db_sum += max_r; }
    }
    db_sum / n_clusters as f64
}

// ── Simulation helper ─────────────────────────────────────────────────────────

/// Inject synthetic contamination at level alpha_sim.
/// X_contaminated[l,n] = clamp(X[l,n] + alpha_sim × μ[l], 0, 1)
pub fn inject_contamination(
    x:         &Array2<f64>,
    locus_mean: &[f64],
    alpha_sim:  f64,
) -> Array2<f64> {
    let (l, n) = (x.nrows(), x.ncols());
    let mut out = x.clone();
    for li in 0..l {
        let add = alpha_sim * locus_mean[li];
        for ni in 0..n {
            let v = out[[li, ni]];
            if !v.is_nan() {
                out[[li, ni]] = (v + add).min(1.0).max(0.0);
            }
        }
    }
    out
}

/// Fraction of covered entries that were floored at 0 during correction.
pub fn floor_fraction(x_raw: &Array2<f64>, x_corr: &Array2<f64>) -> f64 {
    let mut floored = 0usize;
    let mut covered = 0usize;
    for (&r, &c) in x_raw.iter().zip(x_corr.iter()) {
        if !r.is_nan() {
            covered += 1;
            if c == 0.0 && r > 0.0 { floored += 1; }
        }
    }
    if covered > 0 { floored as f64 / covered as f64 } else { 0.0 }
}
