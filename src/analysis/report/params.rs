// ============================================================================
// report/params.rs — §9 Full run parameters table (all 27 params)
// ============================================================================

use crate::analysis::report::ReportData;

pub fn build_params_table(data: &ReportData<'_>) -> String {
    let p = data.params;

    let rows: &[(&str, String, &str)] = &[
        // Inputs
        ("ref_matrix",       p.ref_matrix.clone(),  "Ref-allele UMI count matrix (.mtx) from vartrix Stage 5"),
        ("alt_matrix",       p.alt_matrix.clone(),  "Alt-allele UMI count matrix (.mtx) from vartrix Stage 5"),
        ("clusters",         p.clusters.clone(),    "Cluster assignments + singlet/doublet from souporcell/troublet"),
        ("ambient",          p.ambient.clone(),     "Ambient RNA fraction α from consensus.py Stage 8"),
        // VCF locus selection
        ("souporcell_vcf",   p.souporcell_vcf.as_deref().unwrap_or("(not set)").to_string(),
         "cluster_genotypes.vcf — primary locus source (Mode A or B)"),
        ("freebayes_vcf",    p.freebayes_vcf.as_deref().unwrap_or("(not set)").to_string(),
         "freebayes VCF — enables Mode A exact CHROM:POS cross-reference"),
        // QC fallback
        ("min_ref",          p.min_ref.to_string(),     "Min cells with ≥1 ref read (same as souporcell)"),
        ("min_alt",          p.min_alt.to_string(),     "Min cells with ≥1 alt read (same as souporcell)"),
        ("min_ref_umis",     p.min_ref_umis.to_string(),"Min total ref UMIs per locus (same as souporcell)"),
        ("min_alt_umis",     p.min_alt_umis.to_string(),"Min total alt UMIs per locus (same as souporcell)"),
        // Correction
        ("alpha_override",   p.alpha_override.map(|v| format!("{:.8}", v))
                               .unwrap_or_else(|| "from ambient_rna.txt".to_string()),
         "Override parsed α — blank = use ambient_rna.txt value"),
        ("singlets_only",    p.singlets_only.to_string(),
         "Use only singlet cells for μ estimation (doublets still corrected)"),
        // Output + runtime
        ("output",           p.output.clone(),          "Output directory for all soupXcell files"),
        ("threads",          p.threads.to_string(),     "Parallel threads (Module 2 simulation)"),
        ("seed",             p.seed.to_string(),        "RNG seed for reproducibility"),
        // Module 3 embedding
        ("embed",            format!("{:?}", p.embed),  "Dimensionality reduction method"),
        ("pca_components",   p.pca_components.to_string(), "PCA components before t-SNE/UMAP"),
        ("tsne_perplexity",  p.tsne_perplexity.to_string(), "t-SNE perplexity (neighbourhood size)"),
        ("tsne_iter",        p.tsne_iter.to_string(),   "t-SNE gradient descent iterations"),
        // Module 2 simulation
        ("simulate",         p.simulate.to_string(),    "Enable Module 2 simulation benchmark"),
        ("sim_levels",       format!("{:?}", p.sim_levels), "Contamination levels to simulate"),
        ("sim_trials",       p.sim_trials.to_string(),  "Trials per α level (results averaged)"),
        // Locus selection tuning
        ("metric_cells",     p.metric_cells.to_string(), "Cells per cluster for silhouette subsampling"),
    ];

    let mut t = "<div class='table-wrap'><table>\
        <thead><tr>\
          <th>Parameter</th>\
          <th class='mono' style='font-family:var(--font-mono)'>Value</th>\
          <th>Purpose</th>\
        </tr></thead>\
        <tbody>".to_string();

    for (key, val, desc) in rows {
        t += &format!(
            "<tr>\
             <td class='mono'>{key}</td>\
             <td class='mono'>{val}</td>\
             <td style='color:var(--muted);font-size:12px'>{desc}</td>\
             </tr>"
        );
    }
    t += "</tbody></table></div>";
    t
}
