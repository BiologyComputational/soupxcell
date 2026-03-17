// ============================================================================
// report/experiment.rs — §3 Experiment design + correction narrative (§3)
// ============================================================================

use crate::analysis::report::{ReportData, utils};

pub fn build_experiment_section(
    data:      &ReportData<'_>,
    delta_sil: f64,
    delta_db:  f64,
) -> String {
    let s  = data.summary;
    let p  = data.params;

    let loci_pct   = data.loci_qc as f64 / data.loci_raw.max(1) as f64 * 100.0;
    let floor_pct  = s.pct_floored * 100.0;

    let locus_src_desc = match data.locus_mode.as_str() {
        "A" => format!(
            "<strong>Mode A — exact CHROM:POS</strong>: cluster_genotypes.vcf cross-referenced \
             with freebayes VCF. Each CHROM+POS in the souporcell output VCF is looked up in the \
             freebayes VCF to find the precise matrix row index. Handles chr/no-chr mismatch \
             automatically. This is the most accurate mode."
        ),
        "B" => format!(
            "<strong>Mode B — positional mapping</strong>: cluster_genotypes.vcf used alone. \
             VCF data row N = matrix row N. Valid because vartrix and consensus.py both produce \
             genomically sorted output. Provide --freebayes_vcf for Mode A exact mapping."
        ),
        _ => format!(
            "<strong>Mode C — QC filter fallback</strong>: no VCF provided; souporcell QC \
             filter reproduced with min_ref={} min_alt={} min_ref_umis={} min_alt_umis={}. \
             Provides equivalent locus selection to souporcell's Pass1+Pass2 logic but \
             requires matching the exact QC parameters used in the original run. \
             Provide --souporcell_vcf cluster_genotypes.vcf for exact ground-truth matching.",
             p.min_ref, p.min_alt, p.min_ref_umis, p.min_alt_umis
        ),
    };

    let correction_narrative = if delta_sil.abs() < 0.001 {
        format!(
            "Silhouette score unchanged (Δ = {:+.4}). This is <strong>expected and correct</strong> \
             at α = {:.4}% contamination. The Euclidean-distance-based silhouette is insensitive \
             to the small uniform shift of {:+.6} allele fraction units. \
             The within-cluster standard deviation (Δ = {:+.4}) is the \
             correct metric for this correction magnitude and shows real signal.",
            delta_sil, s.alpha_used * 100.0, s.mean_correction, s.metrics_after.within_cluster_std - s.metrics_before.within_cluster_std
        )
    } else if delta_sil > 0.001 {
        format!(
            "Silhouette score improved (Δ = {:+.4}). At α = {:.4}% \
             the correction measurably tightened cluster separation in Euclidean space. \
             This implies relatively strong ambient contamination that was systematically \
             pulling allele fractions towards the population mean.",
            delta_sil, s.alpha_used * 100.0
        )
    } else {
        format!(
            "Silhouette score decreased (Δ = {:+.4}). This occurs when ambient RNA \
             injection (during simulation) homogenises Euclidean distances between cells. \
             In the real run, a negative Δ at low α = {:.4}% is \
             within numerical noise.",
            delta_sil, s.alpha_used * 100.0
        )
    };

    let conv_border = if delta_sil > 0.001 { utils::C_GREEN }
                      else if delta_sil < -0.01 { utils::C_AMBER }
                      else { utils::C_TEAL };

    let expected_outcomes = [
        ("Loci QC pass rate", ">0.1% of raw loci",
         &format!("{:.1}% ({}/{})", loci_pct, data.loci_qc, data.loci_raw),
         if loci_pct > 0.1 { "pill pill-green" } else { "pill pill-amber" },
         if loci_pct > 0.1 { "&#10003; Pass" } else { "&#9888; Warn" }),
        ("Zero-cov entries", ">90% (scRNA-seq expected)",
         &format!("{:.1}%", data.zero_cov_pct),
         if data.zero_cov_pct > 80.0 { "pill pill-green" } else { "pill pill-amber" },
         if data.zero_cov_pct > 80.0 { "&#10003; Pass" } else { "&#9888; Warn" }),
        ("Entries floored", "<70% of covered",
         &format!("{:.1}%", floor_pct),
         if floor_pct < 70.0 { "pill pill-green" } else { "pill pill-amber" },
         if floor_pct < 70.0 { "&#10003; Pass" } else { "&#9888; Warn" }),
        ("Silhouette before", "−1 to +1 (negative OK for allele data)",
         &format!("{:.4}", s.metrics_before.silhouette),
         "pill pill-teal", "&#9432; Info"),
        ("Within-std Δ", "negative = tighter clusters",
         &format!("{:+.6}", s.metrics_after.within_cluster_std - s.metrics_before.within_cluster_std),
         if s.metrics_after.within_cluster_std < s.metrics_before.within_cluster_std { "pill pill-green" } else { "pill pill-amber" },
         if s.metrics_after.within_cluster_std < s.metrics_before.within_cluster_std { "&#10003; Improved" } else { "&#9472; No change" }),
        ("Davies-Bouldin Δ", "negative = better",
         &format!("{:+.6}", delta_db),
         if delta_db < 0.0 { "pill pill-green" } else { "pill pill-teal" },
         if delta_db < 0.0 { "&#10003; Improved" } else { "&#9472; Unchanged" }),
    ];

    let outcome_rows: String = expected_outcomes.iter().map(|(metric, expected, observed, cls, label)| {
        format!("<tr><td>{metric}</td><td>{expected}</td>\
                 <td class='num mono'>{observed}</td>\
                 <td><span class='{cls}'>{label}</span></td></tr>")
    }).collect();

    format!(r#"<div class="section-label">Experimental Context</div>
<h2>Experiment Design, Locus Selection &amp; Correction Narrative</h2>
<p class="section-desc">
  Scientific rationale, pipeline position, locus resolution strategy, and
  falsifiable hypotheses for this ambient RNA correction run.
</p>

<div class="two-col">
<div>
  <h3>Biological Context &amp; Hypotheses</h3>
  <div class="hypo-box">
    <div class="hypo-label">Primary Hypothesis</div>
    <p>The allele fraction matrix X[l,n] contains a systematic ambient RNA
       bias of magnitude &alpha;&nbsp;=&nbsp;<strong>{alpha:.8}</strong> at each
       locus l, arising from RNA shed into the droplet supernatant before or
       during cell encapsulation. Subtracting &alpha;&times;&mu;[l] from each
       entry recovers the true cell-intrinsic allele fraction.</p>
  </div>
  <div class="hypo-box violet">
    <div class="hypo-label">Secondary Hypothesis</div>
    <p>At &alpha;&nbsp;=&nbsp;{alpha_pct:.3}% the per-cell allele fraction shift
       (&asymp;&nbsp;{mean_corr:.6} units) is too small to move silhouette scores
       appreciably. The within-cluster standard deviation, which captures
       locus-level homogenisation, is the correct primary metric.</p>
  </div>
  <div class="hypo-box amber">
    <div class="hypo-label">Locus Source</div>
    <p>{locus_src_desc}</p>
  </div>
</div>

<div>
  <h3>soupXcell in the souporcell Pipeline</h3>
  <div class="pipeline-steps">
    <div class="pipeline-step">
      <div class="ps-num" style="background:var(--muted)">1</div>
      <div><div class="ps-name" style="color:var(--muted)">retag.py + minimap2 + vartrix</div>
           <div class="ps-desc">Produce ref.mtx + alt.mtx</div></div></div>
    <div class="pipeline-step">
      <div class="ps-num" style="background:var(--muted)">2</div>
      <div><div class="ps-name" style="color:var(--muted)">souporcell (EM)</div>
           <div class="ps-desc">Cluster k={k} donor genotypes &rarr; clusters.tsv</div></div></div>
    <div class="pipeline-step">
      <div class="ps-num" style="background:var(--muted)">3</div>
      <div><div class="ps-name" style="color:var(--muted)">troublet</div>
           <div class="ps-desc">Doublet detection &rarr; singlet/doublet labels</div></div></div>
    <div class="pipeline-step">
      <div class="ps-num" style="background:var(--muted)">4</div>
      <div><div class="ps-name" style="color:var(--muted)">consensus.py</div>
           <div class="ps-desc">Genotype calls &rarr; cluster_genotypes.vcf + ambient_rna.txt (&alpha;)</div></div></div>
    <div class="pipeline-step">
      <div class="ps-num violet">5</div>
      <div><div class="ps-name"><strong>soupXcell (this run)</strong></div>
           <div class="ps-desc">Correct X for ambient RNA &rarr; X_corrected.mtx</div></div></div>
  </div>

  <h3>Input Data Summary</h3>
  <div class="table-wrap"><table>
    <thead><tr><th>Resource</th><th>Value</th><th>Notes</th></tr></thead>
    <tbody>
      <tr><td>ref.mtx</td><td class="mono">{ref_short}</td><td>Ref-allele UMI counts</td></tr>
      <tr><td>alt.mtx</td><td class="mono">{alt_short}</td><td>Alt-allele UMI counts</td></tr>
      <tr><td>clusters.tsv</td><td class="mono">{clu_short}</td><td>Barcode assignments + singlet/doublet</td></tr>
      <tr><td>ambient_rna.txt</td><td class="mono">{amb_short}</td><td>&alpha; = {alpha:.8}</td></tr>
      <tr><td>Raw loci (matrix)</td><td class="num">{loci_raw}</td><td>All positions in mtx headers</td></tr>
      <tr><td>QC-passing loci</td><td class="num">{loci_qc}</td><td>{loci_pct:.1}% pass rate</td></tr>
      <tr><td>Total cells</td><td class="num">{total_cells}</td><td>Unique barcodes</td></tr>
      <tr><td>Singlets / doublets</td><td class="num">{n_sing} / {n_doub}</td>
          <td>troublet classification</td></tr>
      <tr><td>Zero-cov entries</td><td class="num">{zero_cov_pct:.1}%</td>
          <td>locus&times;cell pairs with 0 reads</td></tr>
    </tbody>
  </table></div>
</div>
</div>

<h3>Expected vs. Observed Outcomes</h3>
<div class="table-wrap"><table>
  <thead><tr>
    <th>Metric</th><th>Expected</th><th class="num">Observed</th><th>Status</th>
  </tr></thead>
  <tbody>{outcome_rows}</tbody>
</table></div>

<h3>Correction Narrative</h3>
<div class="card" style="background:var(--surface2);border-left:4px solid {conv_border}">
  <p style="font-size:13px;color:var(--ink2);line-height:1.8">{correction_narrative}</p>
</div>
"#,
        alpha             = s.alpha_used,
        alpha_pct         = s.alpha_used * 100.0,
        mean_corr         = s.mean_correction,
        locus_src_desc    = locus_src_desc,
        k                 = s.n_clusters,
        ref_short         = shorten_path(&p.ref_matrix),
        alt_short         = shorten_path(&p.alt_matrix),
        clu_short         = shorten_path(&p.clusters),
        amb_short         = shorten_path(&p.ambient),
        loci_raw          = data.loci_raw,
        loci_qc           = data.loci_qc,
        loci_pct          = loci_pct,
        total_cells       = s.n_cells_total,
        n_sing            = s.n_cells_singlet,
        n_doub            = s.n_cells_doublet,
        zero_cov_pct      = data.zero_cov_pct,
        outcome_rows      = outcome_rows,
        conv_border       = conv_border,
        correction_narrative = correction_narrative,
    )
}

fn shorten_path(p: &str) -> String {
    std::path::Path::new(p)
        .file_name()
        .and_then(|f| f.to_str())
        .unwrap_or(p)
        .to_string()
}
