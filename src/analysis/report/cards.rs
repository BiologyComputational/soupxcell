// ============================================================================
// report/cards.rs — §2 Executive-summary metric card grid (12 cards)
// ============================================================================

use crate::analysis::report::ReportData;

pub fn build_metric_cards(
    data:       &ReportData<'_>,
    delta_sil:  f64,
    delta_db:   f64,
    delta_wstd: f64,
) -> String {
    let s         = data.summary;
    let floor_pct = s.pct_floored * 100.0;
    let loci_pct  = data.loci_qc as f64 / data.loci_raw.max(1) as f64 * 100.0;
    let total_s   = super::utils::fmt_dur(data.total_time);

    fn sign_color(v: f64, lower_is_better: bool) -> &'static str {
        if lower_is_better {
            if v < -0.001 { "green" } else if v > 0.001 { "amber" } else { "" }
        } else {
            if v > 0.001 { "green" } else if v < -0.001 { "amber" } else { "" }
        }
    }

    fn mc(label: &str, val: &str, sub: &str, cls: &str) -> String {
        let big = val.len() <= 8;
        let val_cls = if big { "mc-value" } else { "mc-value-sm" };
        format!(r#"<div class="metric-card {cls}">
  <div class="mc-label">{label}</div>
  <div class="{val_cls}">{val}</div>
  <div class="mc-sub">{sub}</div>
</div>"#)
    }

    fn mc_delta(label: &str, val: f64, sub: &str, lower_better: bool) -> String {
        let cls = sign_color(val, lower_better);
        let sign = if val >= 0.0 { "+" } else { "" };
        mc(label, &format!("{sign}{val:.4}"), sub, cls)
    }

    let mut out = String::new();

    // Row 1: inputs
    out += &mc("Alpha (α)",
               &format!("{:.6}", s.alpha_used),
               &format!("{:.3}% ambient contamination", s.alpha_used * 100.0), "teal");
    out += &mc("QC Loci",
               &data.loci_qc.to_string(),
               &format!("{:.1}% of {} raw · {}", loci_pct, data.loci_raw, data.locus_mode.as_str()
                    .replace("A","Mode A").replace("B","Mode B").replace("C","Mode C")), "teal");
    out += &mc("Total Cells",
               &s.n_cells_total.to_string(),
               "all barcodes in clusters.tsv", "");
    out += &mc("Singlets / Doublets",
               &format!("{} / {}", s.n_cells_singlet, s.n_cells_doublet),
               "singlets used for μ estimation", "");

    // Row 2: correction
    out += &mc("Entries Floored",
               &format!("{:.1}%", floor_pct),
               &format!("{} of {} covered entries", s.n_floored, s.n_covered),
               if floor_pct > 60.0 { "amber" } else { "" });
    out += &mc("Mean |Correction|",
               &format!("{:.6}", s.mean_correction),
               "avg allele fraction shift per covered entry", "violet");
    out += &mc("Soup Max",
               &format!("{:.6}", s.soup_max),
               "max α×μ[l] across all loci", "");
    out += &mc("Zero-Cov %",
               &format!("{:.1}%", data.zero_cov_pct),
               "locus×cell pairs with 0 reads (expect >90%)",
               if data.zero_cov_pct < 80.0 { "amber" } else { "" });

    // Row 3: quality metrics Δ
    out += &mc_delta("Silhouette Δ",
                     delta_sil,
                     "after − before · positive = improved", false);
    out += &mc_delta("Within-Std Δ",
                     delta_wstd,
                     "negative = tighter clusters · primary metric", true);
    out += &mc_delta("Davies-Bouldin Δ",
                     delta_db,
                     "negative = better cluster separation", true);
    out += &mc("Runtime",
               &total_s,
               &format!("{} threads · seed {}", data.params.threads, data.params.seed), "");

    out
}
