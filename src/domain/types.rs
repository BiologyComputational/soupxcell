
// ── Per-cluster summary for report ───────────────────────────────────────────
#[derive(Clone, Default, Debug)]
pub struct PerClusterSummary {
    pub n_cells:         usize,
    pub n_singlets:      usize,
    pub n_doublets:      usize,
    pub mean_correction: f64,
    pub pct_floored:     f64,
}

// ============================================================================
// domain/types.rs — Core data structures for soupxcell
// ============================================================================
// Pure value types — no I/O, no CLI, no algorithms.
// Every other module imports from here; this module imports nothing internal.
// ============================================================================

use ndarray::Array2;

// ── Cell metadata loaded from clusters.tsv ───────────────────────────────────

#[derive(Debug, Clone, PartialEq)]
pub enum CellStatus {
    Singlet,
    Doublet,
    Unassigned,
}

impl CellStatus {
    pub fn from_str(s: &str) -> Self {
        match s.to_lowercase().as_str() {
            "singlet"  => CellStatus::Singlet,
            "doublet"  => CellStatus::Doublet,
            _          => CellStatus::Unassigned,
        }
    }
    pub fn is_singlet(&self) -> bool { matches!(self, CellStatus::Singlet) }
}

/// Per-cell metadata from clusters.tsv
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct CellMeta {
    pub barcode:    String,
    pub status:     CellStatus,
    pub cluster_id: usize,
    pub log_loss:   f64,
    /// Column index in ref.mtx / alt.mtx (0-based)
    pub col_idx:    usize,
}

// ── Raw loaded matrix pair ────────────────────────────────────────────────────

/// Sparse count matrices loaded from ref.mtx and alt.mtx.
/// Both are stored as dense f64 arrays for computation efficiency at this scale
/// (L=~2913, N=~3628 → ~84 MB — acceptable).
#[derive(Debug)]
pub struct RawMatrices {
    /// Ref allele counts  [L × N]
    pub ref_counts: Array2<f64>,
    /// Alt allele counts  [L × N]
    pub alt_counts: Array2<f64>,
    /// Number of loci (rows)
    pub n_loci: usize,
    /// Number of cells (columns)
    pub n_cells: usize,
}

// ── Allele fraction matrix + soup vector ─────────────────────────────────────

/// Allele fraction matrix X = alt / (ref + alt), NaN where coverage = 0.
/// Shape: [L × N]
#[derive(Debug)]
pub struct AlleleFractionMatrix {
    /// X[l, n] = alt[l,n] / (ref[l,n] + alt[l,n]). NaN where uncovered.
    pub values: Array2<f64>,
    pub n_loci:  usize,
    pub n_cells: usize,
}

/// Per-locus soup contribution: soup[l] = α × mean_X[l]
#[derive(Debug)]
pub struct SoupVector {
    /// Mean allele fraction per locus (ignoring NaN / zero-coverage entries)
    pub locus_mean:        Vec<f64>,
    /// soup[l] = alpha × locus_mean[l]
    pub soup_contribution: Vec<f64>,
    pub alpha:             f64,
    pub n_loci:            usize,
}

// ── Correction result (Module 1) ─────────────────────────────────────────────

#[derive(Debug)]
#[allow(dead_code)]
pub struct CorrectionResult {
    /// Original allele fraction matrix  [L × N]
    pub x_raw:          Array2<f64>,
    /// Corrected allele fraction matrix [L × N]  (floored at 0, NaN preserved)
    pub x_corrected:    Array2<f64>,
    /// The soup vector used
    pub soup_vector:    SoupVector,
    /// α that was applied
    pub alpha_used:     f64,
    /// Number of (locus, cell) pairs that were floored at 0
    pub n_floored:      usize,
    /// Fraction of covered (non-NaN) entries that hit the floor
    pub pct_floored:    f64,
    /// Mean absolute correction applied (over non-NaN entries)
    pub mean_correction: f64,
    /// Total covered (non-NaN) entries in x_raw
    pub n_covered:      usize,
}

// ── Cluster metrics (for silhouette, variance) ────────────────────────────────

#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ClusterMetrics {
    /// Global silhouette coefficient  ∈ [-1, 1]  (higher = better separation)
    pub silhouette:             f64,
    /// Davies-Bouldin index  (lower = better)
    pub davies_bouldin:         f64,
    /// Mean within-cluster standard deviation
    pub mean_within_cluster_std: f64,
    pub within_cluster_std: f64,  // alias — same value, used by report module
    /// Per-cluster silhouette scores
    pub per_cluster_silhouette: Vec<f64>,
    /// Per-cluster cell counts
    pub cluster_counts:         Vec<usize>,
}

// ── Embedding result (Module 3) ───────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq)]
#[allow(dead_code)]
pub enum EmbedMethod {
    PcaTsne,
    PcaUmap,
}

impl EmbedMethod {
    #[allow(dead_code)]
    pub fn name(&self) -> &'static str {
        match self { EmbedMethod::PcaTsne => "t-SNE", EmbedMethod::PcaUmap => "UMAP" }
    }
}

#[derive(Debug)]
#[allow(dead_code)]
pub struct EmbeddingResult {
    /// 2D coordinates  [N × 2]
    pub coords:  Array2<f64>,
    pub method:  EmbedMethod,
    pub n_cells: usize,
}

// ── Simulation benchmark result (Module 2) ────────────────────────────────────

#[derive(Debug)]
pub struct BenchmarkPoint {
    pub alpha_sim:          f64,
    pub rmse:               f64,
    pub silhouette_before:  f64,
    pub silhouette_after:   f64,
    pub floor_pct_before:   f64,
    pub floor_pct_after:    f64,
    pub mean_correction:    f64,
}

#[derive(Debug)]
pub struct BenchmarkResult {
    pub points:      Vec<BenchmarkPoint>,
    pub true_alpha:  f64,
}

// ── Summary statistics for the HTML report ───────────────────────────────────

#[derive(Debug)]
#[allow(dead_code)]
pub struct RunSummary {
    pub n_loci:              usize,
    pub n_cells_total:       usize,
    pub n_cells_singlet:     usize,
    pub n_cells_doublet:     usize,
    pub n_clusters:          usize,
    pub alpha_used:          f64,
    pub n_floored:           usize,
    pub pct_floored:         f64,
    pub mean_correction:     f64,
    pub n_covered:           usize,
    pub soup_max:             f64,
    pub soup_mean:             f64,
    pub metrics_before:      ClusterMetrics,
    pub metrics_after:       ClusterMetrics,
    pub silhouette_delta:    f64,   // after - before
    pub simulation_run:      bool,
    pub per_cluster:         Vec<PerClusterSummary>,
}
