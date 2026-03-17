// ============================================================================
// infra/io.rs — File I/O: loading all inputs and writing all outputs
// ============================================================================

use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write, BufWriter};
use std::path::Path;
use ndarray::Array2;
use crate::domain::types::{RawMatrices, CellMeta, CellStatus, CorrectionResult,
                            BenchmarkResult, SoupVector};
use crate::infra::logger::Logger;
use crate::config::params::expand_tilde;

// ── MatrixMarket (.mtx) loader ────────────────────────────────────────────────

/// Load a MatrixMarket sparse matrix into a dense Array2<f64>.
/// Handles both plain .mtx and .mtx.gz.
pub fn load_mtx(path: &str) -> Result<Array2<f64>, String> {
    let expanded = crate::config::params::expand_tilde(path);
    let file = File::open(&expanded)
        .map_err(|e| format!("Cannot open {}: {}", expanded, e))?;

    let reader: Box<dyn BufRead> = if expanded.ends_with(".gz") {
        Box::new(BufReader::new(flate2::read::GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut lines = reader.lines();
    // Skip comment lines
    let mut header_line = String::new();
    for line in lines.by_ref() {
        let l = line.map_err(|e| e.to_string())?;
        if !l.starts_with('%') {
            header_line = l;
            break;
        }
    }
    // Header: nrows ncols nnz
    let header: Vec<usize> = header_line.split_whitespace()
        .filter_map(|s| s.parse::<usize>().ok())
        .collect();
    if header.len() < 3 {
        return Err(format!("Bad MatrixMarket header in {}: {:?}", expanded, header));
    }
    let (nrows, ncols) = (header[0], header[1]);
    let mut mat = Array2::<f64>::zeros((nrows, ncols));

    for line in lines {
        let l = line.map_err(|e| e.to_string())?;
        let parts: Vec<f64> = l.split_whitespace()
            .filter_map(|s| s.parse::<f64>().ok())
            .collect();
        if parts.len() >= 3 {
            let r = parts[0] as usize - 1;  // 1-indexed → 0-indexed
            let c = parts[1] as usize - 1;
            if r < nrows && c < ncols {
                mat[[r, c]] = parts[2];
            }
        }
    }
    Ok(mat)
}

/// Load both ref.mtx and alt.mtx into a RawMatrices struct.
pub fn load_raw_matrices(ref_path: &str, alt_path: &str, log: &Logger)
    -> Result<RawMatrices, String>
{
    let t = std::time::Instant::now();
    log.info("ref matrix", ref_path);
    let ref_counts = load_mtx(ref_path)?;
    log.info("alt matrix", alt_path);
    let alt_counts = load_mtx(alt_path)?;

    if ref_counts.shape() != alt_counts.shape() {
        return Err(format!(
            "ref.mtx shape {:?} ≠ alt.mtx shape {:?} — matrices must match",
            ref_counts.shape(), alt_counts.shape()
        ));
    }
    let (n_loci, n_cells) = (ref_counts.nrows(), ref_counts.ncols());
    log.info("loci", &n_loci.to_string());
    log.info("cells (columns)", &n_cells.to_string());
    log.ok("Matrices loaded", Some(t));

    Ok(RawMatrices { ref_counts, alt_counts, n_loci, n_cells })
}

// ── clusters.tsv loader ───────────────────────────────────────────────────────

/// Load clusters.tsv produced by souporcell.
///
/// Format (tab-separated):
///   barcode  cluster  status  log_loss   [log_probs...]
/// or the earlier format:
///   barcode  status  cluster  log_loss
///
/// We detect which format by peeking at the first data line.
pub fn load_clusters(path: &str, _n_cells_in_matrix: usize, log: &Logger)
    -> Result<Vec<CellMeta>, String>
{
    let expanded = crate::config::params::expand_tilde(path);
    let file = File::open(&expanded)
        .map_err(|e| format!("Cannot open {}: {}", expanded, e))?;
    let reader = BufReader::new(file);

    let mut metas: Vec<CellMeta> = Vec::new();
    let mut barcode_to_col: std::collections::HashMap<String, usize> = std::collections::HashMap::new();

    // Build a sequential col_idx counter — clusters.tsv rows map to matrix columns
    // in the order they appear in the barcodes file.
    // Since we don't have the original barcodes file here, we assign col_idx
    // by order of appearance, consistent with how souporcell writes them.
    let mut col_counter = 0usize;

    for (_lineno, line) in reader.lines().enumerate() {
        let l = line.map_err(|e| e.to_string())?;
        if l.trim().is_empty() { continue; }
        let parts: Vec<&str> = l.split('\t').collect();
        if parts.len() < 3 { continue; }

        let barcode = parts[0].trim().to_string();

        // Detect format: souporcell v3 outputs barcode cluster log_probs...
        // where col 1 is the integer cluster index and there's no explicit status col.
        // troublet adds a status column. Handle both.
        let (status, cluster_id, log_loss) = if parts.len() >= 4 {
            // Check if parts[1] looks like singlet/doublet
            if parts[1].to_lowercase().contains("singlet")
               || parts[1].to_lowercase().contains("doublet")
               || parts[1].to_lowercase().contains("unassigned") {
                // Format: barcode status cluster log_loss
                let st = CellStatus::from_str(parts[1]);
                let cl = parts[2].trim().parse::<usize>().unwrap_or(0);
                let ll = parts[3].trim().parse::<f64>().unwrap_or(0.0);
                (st, cl, ll)
            } else {
                // Format: barcode cluster log_p0 log_p1 ... (from clusters_tmp.tsv)
                let cl = parts[1].trim().parse::<usize>().unwrap_or(0);
                let ll = parts[2].trim().parse::<f64>().unwrap_or(0.0);
                (CellStatus::Singlet, cl, ll)
            }
        } else {
            let cl = parts[1].trim().parse::<usize>().unwrap_or(0);
            (CellStatus::Singlet, cl, 0.0)
        };

        let col_idx = *barcode_to_col.entry(barcode.clone())
            .or_insert_with(|| { let c = col_counter; col_counter += 1; c });

        metas.push(CellMeta { barcode, status, cluster_id, log_loss, col_idx });
    }

    log.info("cells in clusters.tsv", &metas.len().to_string());
    let singlets = metas.iter().filter(|m| m.status.is_singlet()).count();
    let doublets = metas.len() - singlets;
    log.info("singlets", &singlets.to_string());
    log.info("doublets", &doublets.to_string());

    Ok(metas)
}

// ── ambient_rna.txt loader ────────────────────────────────────────────────────

/// Parse the single global ambient RNA fraction from ambient_rna.txt.
/// Handles: "ambient RNA estimated as 0.012619...\n"
pub fn parse_ambient_fraction(path: &str) -> Result<f64, String> {
    let expanded = crate::config::params::expand_tilde(path);
    let content = std::fs::read_to_string(&expanded)
        .map_err(|e| format!("Cannot read {}: {}", expanded, e))?;

    // Try to extract a float from the content
    for token in content.split_whitespace() {
        if let Ok(v) = token.trim_matches('%').parse::<f64>() {
            if v > 0.0 && v < 100.0 {
                // If it looks like a percentage (> 1.0), convert
                let alpha = if v > 1.0 { v / 100.0 } else { v };
                return Ok(alpha);
            }
        }
    }
    Err(format!("No valid ambient fraction found in {}.\nContent: {}", expanded, content.trim()))
}

// ── Write corrected matrix ────────────────────────────────────────────────────

/// Write a dense Array2<f64> as MatrixMarket format, skipping NaN entries.
pub fn write_corrected_mtx(
    mat:  &Array2<f64>,
    path: &str,
) -> Result<(), String> {
    let p = Path::new(path);
    if let Some(parent) = p.parent() { fs::create_dir_all(parent).ok(); }
    let file = File::create(path)
        .map_err(|e| format!("Cannot create {}: {}", path, e))?;
    let mut w = BufWriter::new(file);

    // Count non-NaN entries
    let nnz = mat.iter().filter(|v| !v.is_nan() && **v != 0.0).count();
    writeln!(w, "%%MatrixMarket matrix coordinate real general").ok();
    writeln!(w, "% Generated by soupxcell v1.0").ok();
    writeln!(w, "{} {} {}", mat.nrows(), mat.ncols(), nnz).ok();

    for ((r, c), &v) in mat.indexed_iter() {
        if !v.is_nan() && v != 0.0 {
            writeln!(w, "{} {} {:.8}", r+1, c+1, v).ok();
        }
    }
    Ok(())
}

/// Write the soup vector to CSV.
pub fn write_soup_vector(soup: &SoupVector, path: &str) -> Result<(), String> {
    let p = Path::new(path);
    if let Some(parent) = p.parent() { fs::create_dir_all(parent).ok(); }
    let file = File::create(path)
        .map_err(|e| format!("Cannot create {}: {}", path, e))?;
    let mut w = BufWriter::new(file);
    writeln!(w, "locus_idx,locus_mean_af,soup_contribution,alpha").ok();
    for i in 0..soup.n_loci {
        writeln!(w, "{},{:.8},{:.8},{:.8}",
            i, soup.locus_mean[i], soup.soup_contribution[i], soup.alpha).ok();
    }
    Ok(())
}

/// Write per-cell correction summary to CSV.
pub fn write_correction_summary(
    result: &CorrectionResult,
    cells:  &[CellMeta],
    path:   &str,
) -> Result<(), String> {
    let p = Path::new(path);
    if let Some(parent) = p.parent() { fs::create_dir_all(parent).ok(); }
    let file = File::create(path)
        .map_err(|e| format!("Cannot create {}: {}", path, e))?;
    let mut w = BufWriter::new(file);
    writeln!(w, "barcode,cluster_id,status,mean_subtracted,n_loci_floored,n_loci_covered").ok();

    let (l, _) = (result.x_raw.nrows(), result.x_raw.ncols());

    for cell in cells {
        let ni = cell.col_idx;
        if ni >= result.x_raw.ncols() { continue; }

        let mut mean_sub = 0.0_f64;
        let mut n_floor  = 0usize;
        let mut n_cover  = 0usize;

        for li in 0..l {
            let raw  = result.x_raw[[li, ni]];
            let corr = result.x_corrected[[li, ni]];
            if !raw.is_nan() {
                n_cover += 1;
                mean_sub += (raw - corr).abs();
                if raw > 0.0 && corr == 0.0 { n_floor += 1; }
            }
        }
        let mean_sub = if n_cover > 0 { mean_sub / n_cover as f64 } else { 0.0 };

        writeln!(w, "{},{},{},{:.8},{},{}",
            cell.barcode,
            cell.cluster_id,
            match cell.status {
                CellStatus::Singlet    => "singlet",
                CellStatus::Doublet    => "doublet",
                CellStatus::Unassigned => "unassigned",
            },
            mean_sub, n_floor, n_cover).ok();
    }
    Ok(())
}

/// Write benchmark results to CSV.
pub fn write_benchmark_csv(result: &BenchmarkResult, path: &str) -> Result<(), String> {
    let p = Path::new(path);
    if let Some(parent) = p.parent() { fs::create_dir_all(parent).ok(); }
    let file = File::create(path)
        .map_err(|e| format!("Cannot create {}: {}", path, e))?;
    let mut w = BufWriter::new(file);
    writeln!(w, "alpha_sim,rmse,silhouette_before,silhouette_after,silhouette_delta,floor_pct_before,floor_pct_after,mean_correction").ok();
    for pt in &result.points {
        writeln!(w, "{:.4},{:.6},{:.4},{:.4},{:.4},{:.4},{:.4},{:.6}",
            pt.alpha_sim, pt.rmse,
            pt.silhouette_before, pt.silhouette_after,
            pt.silhouette_after - pt.silhouette_before,
            pt.floor_pct_before, pt.floor_pct_after,
            pt.mean_correction).ok();
    }
    Ok(())
}

/// Write metrics.txt summary.
pub fn write_metrics_txt(
    alpha:      f64,
    n_loci:     usize,
    n_cells:    usize,
    n_floored:  usize,
    pct_floored: f64,
    mean_corr:  f64,
    sil_before: f64,
    sil_after:  f64,
    db_before:  f64,
    db_after:   f64,
    path:       &str,
) -> Result<(), String> {
    let p = Path::new(path);
    if let Some(parent) = p.parent() { fs::create_dir_all(parent).ok(); }
    let file = File::create(path)
        .map_err(|e| format!("Cannot create {}: {}", path, e))?;
    let mut w = BufWriter::new(file);
    writeln!(w, "soupxcell v1.0 — Correction Metrics Summary").ok();
    writeln!(w, "============================================").ok();
    writeln!(w, "alpha_used            : {:.8}", alpha).ok();
    writeln!(w, "n_loci                : {}", n_loci).ok();
    writeln!(w, "n_cells               : {}", n_cells).ok();
    writeln!(w, "n_entries_floored     : {}", n_floored).ok();
    writeln!(w, "pct_floored           : {:.2}%", pct_floored * 100.0).ok();
    writeln!(w, "mean_abs_correction   : {:.6}", mean_corr).ok();
    writeln!(w, "silhouette_before     : {:.4}", sil_before).ok();
    writeln!(w, "silhouette_after      : {:.4}", sil_after).ok();
    writeln!(w, "silhouette_delta      : {:+.4}", sil_after - sil_before).ok();
    writeln!(w, "davies_bouldin_before : {:.4}", db_before).ok();
    writeln!(w, "davies_bouldin_after  : {:.4}", db_after).ok();
    writeln!(w, "davies_bouldin_delta  : {:+.4}", db_after - db_before).ok();
    Ok(())
}

// ── souporcell QC locus helpers ───────────────────────────────────────────────

/// Parse "Loci passing QC: N" from a souporcell clusters.err file.
/// Returns N if found, None otherwise.
#[allow(dead_code)]
pub fn parse_qc_loci_from_err(err_path: &str) -> Option<usize> {
    let expanded = crate::config::params::expand_tilde(err_path);
    let content  = std::fs::read_to_string(&expanded).ok()?;
    for line in content.lines() {
        // Matches both:
        //   "Loci passing QC:   8219"     (souporcell v3.0 format)
        //   "total loci used 2913"         (souporcell v2.4 format)
        if line.contains("Loci passing QC") || line.contains("loci used") {
            if let Some(n) = line.split_whitespace().last()
                .and_then(|s| s.parse::<usize>().ok()) {
                return Some(n);
            }
        }
    }
    None
}

/// Count data rows (non-header, non-comment) in a VCF file.
/// Each row corresponds to one locus that souporcell kept after QC.
/// The row count = the number to pass as n_qc_loci.
///
/// Also returns the sorted list of 0-based row indices if the caller needs
/// them for position-matched filtering.
#[allow(dead_code)]
pub fn count_vcf_loci(vcf_path: &str) -> Result<usize, String> {
    use std::io::BufRead;
    let expanded = crate::config::params::expand_tilde(vcf_path);
    let file = std::fs::File::open(&expanded)
        .map_err(|e| format!("Cannot open VCF {}: {}", expanded, e))?;
    let reader: Box<dyn BufRead> = if expanded.ends_with(".gz") {
        Box::new(std::io::BufReader::new(flate2::read::GzDecoder::new(file)))
    } else {
        Box::new(std::io::BufReader::new(file))
    };
    let count = reader.lines()
        .filter_map(|l| l.ok())
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .count();
    Ok(count)
}

// ── QC-aware matrix loader (mirrors souporcell v3.0 load_cell_data logic) ────
//
// Loads ref.mtx + alt.mtx in lockstep (guaranteed same ordering from vartrix).
// Applies the same 4-parameter QC filter that souporcell uses:
//   min_ref      — min cells with ≥1 ref read   (default: 4)
//   min_alt      — min cells with ≥1 alt read   (default: 4)
//   min_ref_umis — min total ref UMIs across all cells (default: 4)
//   min_alt_umis — min total alt UMIs across all cells (default: 4)
//
// Returns:
//   x_raw:          Array2<f64>  [L_qc × N]   allele fraction matrix, NaN where uncovered
//   index_to_locus: Vec<usize>   compact_idx → original matrix row index
//   loci_raw:       usize        total loci in matrix before QC
//   loci_used:      usize        loci passing QC  (= souporcell's "Loci passing QC: N")
//   zero_cov_entries: usize      (locus,cell) pairs with no reads after QC
pub fn load_matrices_qc(
    ref_path:     &str,
    alt_path:     &str,
    min_ref:      u32,
    min_alt:      u32,
    min_ref_umis: u32,
    min_alt_umis: u32,
    log:          &crate::infra::logger::Logger,
) -> Result<(ndarray::Array2<f64>, Vec<usize>, usize, usize, usize), String> {
    use std::io::BufRead;
    use std::collections::HashMap as StdMap;

    log.computing(&format!("loading ref+alt matrices in lockstep with QC filter (min_ref={} min_alt={} min_ref_umis={} min_alt_umis={})",
        min_ref, min_alt, min_ref_umis, min_alt_umis));

    // Open files — expand ~ and handle .gz transparently
    let ref_file = std::fs::File::open(&expand_tilde(ref_path))
        .map_err(|e| format!("Cannot open ref matrix '{}': {}", ref_path, e))?;
    let alt_file = std::fs::File::open(&expand_tilde(alt_path))
        .map_err(|e| format!("Cannot open alt matrix '{}': {}", alt_path, e))?;
    let mut ref_lines = std::io::BufReader::with_capacity(256*1024, ref_file).lines();
    let mut alt_lines = std::io::BufReader::with_capacity(256*1024, alt_file).lines();

    let mut total_loci  = 0usize;
    let mut total_cells = 0usize;

    // Per-locus: [cells_with_ref, cells_with_alt]
    let mut locus_cell_counts: StdMap<usize, [u32;2]> = StdMap::new();
    // Per-locus: [total_ref_umis, total_alt_umis]
    let mut locus_umi_counts:  StdMap<usize, [u32;2]> = StdMap::new();
    // Per-locus per-cell: [ref_count, alt_count]
    let mut locus_counts: StdMap<usize, StdMap<usize,[u32;2]>> = StdMap::new();
    let mut all_loci = std::collections::BTreeSet::new();

    // ── Pass 1: aggregate per-locus stats ────────────────────────────────────
    for line_num in 0.. {
        let alt_line = match alt_lines.next() {
            Some(Ok(l)) => l, Some(Err(e)) => return Err(e.to_string()), None => break,
        };
        let ref_line = match ref_lines.next() {
            Some(Ok(l)) => l, Some(Err(e)) => return Err(e.to_string()),
            None => return Err("ref.mtx ended before alt.mtx".to_string()),
        };
        if line_num < 2 { continue; }  // MatrixMarket header lines
        if line_num == 2 {              // dimension line: rows cols nnz
            let t: Vec<&str> = alt_line.split_whitespace().collect();
            if t.len() < 2 { return Err("Malformed MTX dimension line".to_string()); }
            total_loci  = t[0].parse::<usize>().map_err(|e| e.to_string())?;
            total_cells = t[1].parse::<usize>().map_err(|e| e.to_string())?;
            continue;
        }
        // Data line: row col value  (1-indexed)
        let at: Vec<&str> = alt_line.split_whitespace().collect();
        let rt: Vec<&str> = ref_line.split_whitespace().collect();
        if at.len() < 3 || rt.len() < 3 { continue; }
        let locus: usize = at[0].parse::<usize>().map_err(|e| e.to_string())? - 1;
        let cell:  usize = at[1].parse::<usize>().map_err(|e| e.to_string())? - 1;
        let alt_c: u32   = at[2].parse::<u32>().unwrap_or(0);
        let ref_c: u32   = rt[2].parse::<u32>().unwrap_or(0);
        all_loci.insert(locus);
        let cc = locus_cell_counts.entry(locus).or_insert([0u32;2]);
        let uc = locus_umi_counts.entry(locus).or_insert([0u32;2]);
        if ref_c > 0 { cc[0] += 1; uc[0] += ref_c; }
        if alt_c > 0 { cc[1] += 1; uc[1] += alt_c; }
        locus_counts.entry(locus).or_default().insert(cell, [ref_c, alt_c]);
    }

    log.info("loci (raw in matrix)",   &total_loci.to_string());
    log.info("cells (raw in matrix)",  &total_cells.to_string());

    // ── Pass 2: QC filter (same logic as souporcell v3.0 load_cell_data) ─────
    let mut index_to_locus: Vec<usize> = Vec::new();
    for &locus in &all_loci {
        let cc = locus_cell_counts.get(&locus).copied().unwrap_or([0;2]);
        let uc = locus_umi_counts.get(&locus).copied().unwrap_or([0;2]);
        if cc[0] >= min_ref && cc[1] >= min_alt
            && uc[0] >= min_ref_umis && uc[1] >= min_alt_umis {
            index_to_locus.push(locus);
        }
    }
    // index_to_locus is already sorted because BTreeSet is sorted
    let loci_used = index_to_locus.len();
    log.info("loci (passing QC)", &format!("{} / {} ({:.1}% of raw loci)",
        loci_used, total_loci, loci_used as f64 / total_loci.max(1) as f64 * 100.0));
    log.info("QC filter applied", &format!("min_ref={} min_alt={} min_ref_umis={} min_alt_umis={}",
        min_ref, min_alt, min_ref_umis, min_alt_umis));

    if loci_used == 0 {
        return Err(format!(
            "No loci passed QC filters (min_ref={} min_alt={} min_ref_umis={} min_alt_umis={}).
             Try lowering these values or use --souporcell_vcf to import loci directly.",
            min_ref, min_alt, min_ref_umis, min_alt_umis));
    }

    // ── Build QC-filtered allele fraction matrix [L_qc × N] ─────────────────
    let mut x_raw = ndarray::Array2::<f64>::from_elem((loci_used, total_cells), f64::NAN);
    let mut zero_cov_entries = 0usize;
    for (compact_idx, &orig_locus) in index_to_locus.iter().enumerate() {
        if let Some(cell_map) = locus_counts.get(&orig_locus) {
            for cell in 0..total_cells {
                match cell_map.get(&cell) {
                    None | Some(&[0,0]) => { zero_cov_entries += 1; }
                    Some(&[ref_c, alt_c]) => {
                        let total = (ref_c + alt_c) as f64;
                        x_raw[[compact_idx, cell]] = alt_c as f64 / total;
                    }
                }
            }
        } else {
            zero_cov_entries += total_cells;
        }
    }

    let zero_pct = 100.0 * zero_cov_entries as f64 / (loci_used * total_cells).max(1) as f64;
    log.info("zero-cov entries", &format!("{} ({:.1}% of QC-loci × cells — expected >95% for sparse scRNA)",
        zero_cov_entries, zero_pct));
    if zero_pct < 80.0 {
        log.warn(&format!("Sparsity {:.1}% is lower than expected for scRNA-seq — check input matrices", zero_pct));
    }

    Ok((x_raw, index_to_locus, total_loci, loci_used, zero_cov_entries))
}

// ── VCF-based locus index resolution ─────────────────────────────────────────
//
// The cluster_genotypes.vcf produced by souporcell consensus.py is the
// DEFINITIVE record of which loci souporcell actually used.
// Each data row (non-#) = one locus. CHROM+POS identifies it.
//
// ref.mtx rows are numbered by vartrix in the same sorted CHROM:POS order
// as the freebayes/common_variants VCF that was passed to vartrix.
//
// Resolution:
//   Mode A (exact, recommended): both VCFs provided
//     - Parse CHROM+POS from cluster_genotypes.vcf → position set
//     - Parse CHROM+POS from freebayes VCF → positional index
//     - Return sorted Vec<usize> of 0-based matrix row indices
//
//   Mode B (positional, valid when VCFs share sort order):
//     - cluster_genotypes.vcf data row i → matrix row i (0-based)
//     - Valid because vartrix + consensus.py both use genomically sorted order
//
// Returns: sorted Vec<usize> of 0-based matrix row indices
pub fn parse_vcf_locus_indices(
    cluster_genotypes_vcf: &str,
    freebayes_vcf:         Option<&str>,
    log:                   &crate::infra::logger::Logger,
) -> Result<Vec<usize>, String> {
    use std::io::BufRead;

    // Helper: parse (CHROM, POS) data rows from a VCF file
    fn read_chrom_pos(path: &str) -> Result<Vec<(String, u64)>, String> {
        let expanded = expand_tilde(path);
        let file = std::fs::File::open(&expanded)
            .map_err(|e| format!("Cannot open VCF '{}': {}", path, e))?;
        let reader: Box<dyn BufRead> = if expanded.ends_with(".gz") {
            Box::new(std::io::BufReader::new(
                flate2::read::GzDecoder::new(file)
            ))
        } else {
            Box::new(std::io::BufReader::new(file))
        };
        let mut positions = Vec::new();
        for line in reader.lines() {
            let l = line.map_err(|e| e.to_string())?;
            if l.starts_with('#') || l.trim().is_empty() { continue; }
            let cols: Vec<&str> = l.splitn(5, '\t').collect();
            if cols.len() < 2 { continue; }
            let chrom = cols[0].to_string();
            let pos   = cols[1].parse::<u64>().map_err(|e| e.to_string())?;
            positions.push((chrom, pos));
        }
        Ok(positions)
    }

    let cg_positions = read_chrom_pos(cluster_genotypes_vcf)?;
    let n_vcf = cg_positions.len();
    log.info("VCF loci (cluster_genotypes.vcf)", &n_vcf.to_string());

    match freebayes_vcf {
        // ── Mode A: cross-reference with freebayes VCF for exact row indices ─
        Some(fb_path) => {
            log.info("locus resolution mode", "A — exact CHROM:POS cross-reference");
            let fb_positions = read_chrom_pos(fb_path)?;
            log.info("freebayes VCF total loci", &fb_positions.len().to_string());

            // Build position → 0-based matrix row index from freebayes VCF
            use std::collections::HashMap;
            let fb_index: HashMap<(String,u64), usize> = fb_positions
                .into_iter().enumerate()
                .map(|(i,(c,p))| ((c,p), i))
                .collect();

            // Map each cluster_genotypes position to its matrix row
            let mut indices: Vec<usize> = Vec::with_capacity(n_vcf);
            let mut not_found = 0usize;
            for (chrom, pos) in &cg_positions {
                match fb_index.get(&(chrom.clone(), *pos)) {
                    Some(&row_idx) => indices.push(row_idx),
                    None => {
                        not_found += 1;
                        // Try stripping "chr" prefix mismatch
                        let alt_chrom = if chrom.starts_with("chr") {
                            chrom[3..].to_string()
                        } else {
                            format!("chr{}", chrom)
                        };
                        if let Some(&row_idx) = fb_index.get(&(alt_chrom, *pos)) {
                            indices.push(row_idx);
                            not_found -= 1;
                        }
                    }
                }
            }
            if not_found > 0 {
                log.warn(&format!("{} loci in cluster_genotypes.vcf not found in freebayes VCF — \
                    possible chromosome naming mismatch (chr1 vs 1)", not_found));
            }
            indices.sort_unstable();
            log.info("matrix row indices resolved", &indices.len().to_string());
            Ok(indices)
        }

        // ── Mode B: positional mapping (VCF row i → matrix row i) ────────────
        None => {
            log.info("locus resolution mode",
                "B — positional (VCF row N = matrix row N, assumes same sort order)");
            log.info("tip",
                "for exact mapping provide --freebayes_vcf souporcell_merged_sorted_vcf.vcf.gz");
            // VCF data row i (0-based) = matrix row i (0-based)
            // This is valid because vartrix and consensus.py both use genomic sort order
            let indices: Vec<usize> = (0..n_vcf).collect();
            Ok(indices)
        }
    }
}
