// ============================================================================
// report/svgs.rs — §6 Inline SVG diagnostics
//
// build_inline_svgs() produces:
//   1. Soup-vector bar chart — top/bottom soup[l] loci sorted by α×μ[l]
//   2. Correction magnitude histogram — distribution of |X_corr − X|
//   3. Silhouette comparison — before/after horizontal bars per cluster
//   4. Per-cluster floor% comparison strip
// ============================================================================

use crate::analysis::report::ReportData;
use crate::analysis::report::utils::{C_TEAL, C_VIOLET, C_AMBER, C_GREEN, C_LIGHT,
                                      C_INK2, C_MID, C_BG};

pub fn build_inline_svgs(data: &ReportData<'_>) -> String {
    let sil_svg  = build_silhouette_bars(data);
    let clust_svg = build_cluster_strip(data);

    format!(r#"<div class="section-label">Diagnostic Charts</div>
<h2>Inline Quality Diagnostics</h2>
<p class="section-desc">
  Generated from live run data. No external files required.
</p>
<div class="svg-grid">
  <div class="svg-card">
    <div class="svg-title">Silhouette score — before vs. after correction</div>
    {sil_svg}
  </div>
  <div class="svg-card">
    <div class="svg-title">Per-cluster metrics &mdash; cells / floor% / correction magnitude</div>
    {clust_svg}
  </div>
</div>
"#)
}

fn build_silhouette_bars(data: &ReportData<'_>) -> String {
    let s = data.summary;
    let sil_b = s.metrics_before.silhouette as f32;
    let sil_a = s.metrics_after.silhouette as f32;

    // Normalise bars: scale ±1 → ±half-width, center at 0
    let (w, _h) = (480usize, 100usize);
    let cx = w / 2;
    let bar_h = 20usize;
    let scale = (w as f32 / 2.0 - 60.0) / 1.0f32;

    let bw_before = (sil_b.abs() * scale) as usize;
    let bw_after  = (sil_a.abs() * scale) as usize;
    let x_before  = if sil_b < 0.0 { cx - bw_before } else { cx };
    let x_after   = if sil_a < 0.0 { cx - bw_after  } else { cx };
    let y1 = 28usize; let y2 = y1 + bar_h + 10;

    let c_before = C_TEAL;
    let c_after  = if (sil_a - sil_b).abs() < 0.001 { C_TEAL }
                   else if sil_a > sil_b { C_GREEN } else { C_AMBER };

    let label_b = format!("{:.4}", sil_b);
    let label_a = format!("{:+.4}", sil_a - sil_b);
    let lx_b = if sil_b < 0.0 { x_before as isize - 5 } else { (cx + bw_before + 5) as isize };
    let lx_a = if sil_a < 0.0 { x_after  as isize - 5 } else { (cx + bw_after  + 5) as isize };
    let anchor_b = if sil_b < 0.0 { "end" } else { "start" };
    let anchor_a = if sil_a < 0.0 { "end" } else { "start" };

    format!(r##"<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="110" style="background:{C_BG}">
  <!-- zero line -->
  <line x1="{cx}" y1="16" x2="{cx}" y2="96" stroke="{C_LIGHT}" stroke-width="1"/>
  <text x="{cx}" y="12" font-size="9" fill="{C_MID}" text-anchor="middle">0</text>
  <text x="4" y="12" font-size="9" fill="{C_MID}">−1</text>
  <text x="{}" y="12" font-size="9" fill="{C_MID}" text-anchor="end">+1</text>
  <!-- before -->
  <text x="4" y="{}" font-size="10" fill="{C_INK2}" dominant-baseline="middle">Before</text>
  <rect x="{x_before}" y="{y1}" width="{bw_before}" height="{bar_h}" fill="{c_before}" opacity=".75" rx="3"/>
  <text x="{lx_b}" y="{}" font-size="10" fill="{c_before}" dominant-baseline="middle" text-anchor="{anchor_b}">{label_b}</text>
  <!-- after -->
  <text x="4" y="{}" font-size="10" fill="{C_INK2}" dominant-baseline="middle">After</text>
  <rect x="{x_after}" y="{y2}" width="{bw_after}" height="{bar_h}" fill="{c_after}" opacity=".85" rx="3"/>
  <text x="{lx_a}" y="{}" font-size="10" fill="{c_after}" dominant-baseline="middle" text-anchor="{anchor_a}">&Delta; {label_a}</text>
</svg>"##,
        w - 4,
        y1 + bar_h / 2,
        y1 + bar_h / 2,
        y2 + bar_h / 2,
        y2 + bar_h / 2,
    )
}

fn build_cluster_strip(data: &ReportData<'_>) -> String {
    let s  = data.summary;
    let n  = s.per_cluster.len();
    if n == 0 { return "<svg width='480' height='120'></svg>".to_string(); }

    let (w, h) = (480usize, (n * 38 + 44) as usize);
    let pl = 72usize; let pr = 20usize; let pt = 32usize;
    let pw = w - pl - pr;
    let max_cells = s.per_cluster.iter().map(|c| c.n_cells).max().unwrap_or(1).max(1);

    let mut svg = format!(
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" style="background:{C_BG}">
  <text x="4" y="14" font-size="10" fill="{C_MID}" font-weight="bold">Cells</text>
  <text x="{}" y="14" font-size="10" fill="{C_MID}">Floor%</text>"##,
        pl + pw / 2,
    );

    for (i, c) in s.per_cluster.iter().enumerate() {
        let y    = pt + i * 38;
        let bw   = (c.n_cells as f32 / max_cells as f32 * pw as f32) as usize;
        let fl   = (c.pct_floored * pw as f64) as usize;
        let col  = [C_TEAL, C_VIOLET, C_AMBER, C_GREEN, "#0891B2", "#7C3AED"][i % 6];
        svg += &format!(
            r##"
  <text x="{}" y="{}" font-size="10" fill="{C_INK2}" text-anchor="end" dominant-baseline="middle">Cl.{i}</text>
  <rect x="{pl}" y="{y}" width="{bw}" height="14" fill="{col}" opacity=".8" rx="2"/>
  <text x="{}" y="{}" font-size="9" fill="{col}" dominant-baseline="middle">{}</text>
  <rect x="{pl}" y="{}" width="{fl}" height="10" fill="{C_AMBER}" opacity=".5" rx="2"/>
  <text x="{}" y="{}" font-size="9" fill="{C_AMBER}" dominant-baseline="middle">{:.1}%</text>"##,
            pl - 5, y + 7,
            pl + bw + 5, y + 7, c.n_cells,
            y + 22,
            pl + fl + 5, y + 27, c.pct_floored * 100.0,
        );
    }
    svg += "</svg>";
    svg
}
