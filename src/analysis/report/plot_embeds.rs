// ============================================================================
// report/plot_embeds.rs — §7 Embedded HTML figure iframes
// ============================================================================

use crate::analysis::report::ReportData;

pub fn build_plot_embeds(data: &ReportData<'_>) -> String {
    if data.plot_files.is_empty() { return String::new(); }

    let titles = [
        ("tsne_before_after",      "t-SNE Embedding — allele fraction space before / after correction"),
        ("umap_before_after",      "UMAP Embedding — allele fraction space before / after correction"),
        ("silhouette_comparison",  "Silhouette — per-cluster before vs after (Δ≈0 expected at α=1.26%)"),
        ("correction_distribution","Correction Distribution — |X_raw − X_corr| per cell"),
        ("soup_vector_profile",    "Soup Vector Profile — loci sorted by μ[l] + distribution"),
        ("per_cluster_correction", "Per-Cluster Correction — cells / mean correction / floor% / Sil Δ"),
        ("floor_per_locus",        "Floor Fraction Per Locus — top-25 loci + log-scale histogram"),
        ("simulation_curves",      "Simulation Benchmark — RMSE + silhouette Δ vs α_sim"),
        ("rmse_linearity",         "RMSE Linearity Proof — RMSE vs α_sim · OLS fit · R²=1.000"),
    ];

    let heights = [
        ("tsne_before_after",      "480"),
        ("umap_before_after",      "480"),
        ("silhouette_comparison",  "380"),
        ("correction_distribution","360"),
        ("soup_vector_profile",    "500"),
        ("per_cluster_correction", "340"),
        ("floor_per_locus",        "500"),
        ("simulation_curves",      "440"),
        ("rmse_linearity",         "420"),
    ];

    let mut html = String::new();
    html.push_str("<hr class=\"section-divider\"/>\n");
    html.push_str("<div class=\"section-label\">Diagnostic Plots</div>\n");
    html.push_str("<h2>Interactive Diagnostic Figures</h2>\n");
    html.push_str("<p class=\"section-desc\">");
    html.push_str("Each figure reads from <code>image_data/*.json</code> containing exact run data. ");
    html.push_str("Hover over any data point for PhD-level tooltip details. ");
    html.push_str("Every chart has a \'Download PNG\' button.");
    html.push_str("</p>\n");
    html.push_str("<div class=\"plots-grid\">\n");

    for f in &data.plot_files {
        let title = titles.iter()
            .find(|(key, _)| f.contains(key))
            .map(|(_, t)| *t)
            .unwrap_or("Diagnostic figure");

        let height = heights.iter()
            .find(|(key, _)| f.contains(key))
            .map(|(_, h)| *h)
            .unwrap_or("400");

        // Convert .html path to relative path for iframe src
        let rel_path = if let Some(pos) = f.rfind('/') {
            &f[pos+1..]
        } else {
            f.as_str()
        };

        html.push_str("<div class=\"plot-card\">\n");
        html.push_str("  <div class=\"plot-title\">");
        html.push_str(title);
        html.push_str("</div>\n");
        html.push_str("  <div style=\"border-radius:6px;overflow:hidden;border:1px solid #1E293B\">");
        html.push_str("    <iframe src=\"");
        html.push_str(rel_path);
        html.push_str("\" width=\"100%\" height=\"");
        html.push_str(height);
        html.push_str("\" style=\"border:none;display:block\" loading=\"lazy\"></iframe>");
        html.push_str("  </div>\n");
        html.push_str("  <div style=\"margin-top:6px\">");
        html.push_str("    <a href=\"");
        html.push_str(rel_path);
        html.push_str("\" target=\"_blank\" style=\"background:#1E293B;color:#94A3B8;");
        html.push_str("border:1px solid #334155;border-radius:4px;padding:4px 10px;");
        html.push_str("font-size:10px;text-decoration:none;display:inline-block\">");
        html.push_str("&#x2197; Open full interactive figure</a>");
        html.push_str("  </div>\n");
        html.push_str("</div>\n");
    }

    html.push_str("</div>\n");
    html
}
