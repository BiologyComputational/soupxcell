// ============================================================================
// analysis/plots.rs — Pure SVG generation for all soupxcell visualisations
// ============================================================================
// Five plots, all produced as self-contained SVG strings.
// No external renderer — same approach as souporcell.
//
// 1. tsne_before_after.svg      — side-by-side scatter before/after correction
// 2. umap_before_after.svg      — same structure, UMAP embedding
// 3. silhouette_comparison.svg  — grouped bar chart per-cluster, before/after
// 4. correction_distribution.svg— histogram of per-cell mean correction amount
// 5. simulation_curves.svg      — RMSE + silhouette vs α_sim level
// ============================================================================

use ndarray::Array2;
use crate::domain::types::BenchmarkResult;

// ── Colour palette ────────────────────────────────────────────────────────────

const CLUSTER_COLOURS: &[&str] = &[
    "#06B6D4", "#8B5CF6", "#10B981", "#F59E0B",
    "#EF4444", "#EC4899", "#3B82F6", "#84CC16",
    "#F97316", "#14B8A6", "#A855F7", "#E11D48",
    "#0EA5E9", "#D97706", "#059669", "#7C3AED",
];
const BG:      &str = "#05080F";
const SURFACE: &str = "#090D18";
const BORDER:  &str = "#1C2A3E";
const INK:     &str = "#E8EEF8";
const INK2:    &str = "#94A3B8";
const INK3:    &str = "#4A6080";
const VIOLET:  &str = "#8B5CF6";
const CYAN:    &str = "#06B6D4";
const AMBER:   &str = "#F59E0B";
const GREEN:   &str = "#10B981";
const RED:     &str = "#EF4444";

fn cluster_colour(k: usize) -> &'static str {
    CLUSTER_COLOURS[k % CLUSTER_COLOURS.len()]
}

// ── Axis helpers ─────────────────────────────────────────────────────────────

fn data_range(coords: &Array2<f64>, col: usize) -> (f64, f64) {
    let vals: Vec<f64> = coords.column(col).iter().cloned()
        .filter(|v| v.is_finite()).collect();
    if vals.is_empty() { return (-1.0, 1.0); }
    let min = vals.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let pad = (max - min) * 0.05;
    (min - pad, max + pad)
}

fn data_to_screen(v: f64, dmin: f64, dmax: f64, smin: f64, smax: f64) -> f64 {
    if (dmax - dmin).abs() < 1e-10 { return (smin + smax) / 2.0; }
    smin + (v - dmin) / (dmax - dmin) * (smax - smin)
}

// ── SVG base ─────────────────────────────────────────────────────────────────

fn svg_open(w: usize, h: usize) -> String {
    format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" style="background:{BG};font-family:'IBM Plex Mono',monospace">"#,
        w=w, h=h
    )
}

fn svg_close() -> &'static str { "</svg>" }

fn svg_rect(x: f64, y: f64, w: f64, h: f64, fill: &str, rx: f64) -> String {
    format!(r#"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" fill="{}" rx="{:.1}"/>"#,
        x, y, w, h, fill, rx)
}

fn svg_text(x: f64, y: f64, text: &str, size: f64, fill: &str, anchor: &str) -> String {
    format!(r#"<text x="{:.1}" y="{:.1}" font-size="{:.1}" fill="{}" text-anchor="{}" dominant-baseline="middle">{}</text>"#,
        x, y, size, fill, anchor, text)
}

fn svg_line(x1: f64, y1: f64, x2: f64, y2: f64, stroke: &str, width: f64) -> String {
    format!(r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{}" stroke-width="{:.1}"/>"#,
        x1, y1, x2, y2, stroke, width)
}

// ── PLOT 1 & 2: Before/After scatter (t-SNE or UMAP) ────────────────────────

/// Generate a side-by-side before/after scatter plot.
/// `method_name` = "t-SNE" or "UMAP"
pub fn plot_before_after(
    before:      &Array2<f64>,  // N × 2
    after:       &Array2<f64>,  // N × 2
    labels:      &[usize],
    n_clusters:   usize,
    method_name:  &str,
    sil_before:   f64,
    sil_after:    f64,
    alpha:        f64,
) -> String {
    // ── Publication-quality layout ─────────────────────────────────────────
    // White background, proper tick marks, convex hull outlines per cluster,
    // legend inside plot area, neutral delta annotation.
    let total_w = 1240usize;
    let total_h = 580usize;

    // Panel geometry
    let panel_w  = 540.0_f64;
    let panel_h  = 480.0_f64;
    let panel_gap = 60.0_f64;
    let margin_l  = 50.0_f64;
    let margin_t  = 60.0_f64;
    let pad_l = 55.0_f64;  // inside panel
    let pad_b = 45.0_f64;
    let pad_t = 48.0_f64;
    let pad_r = 20.0_f64;
    let plot_w = panel_w - pad_l - pad_r;
    let plot_h = panel_h - pad_t - pad_b;

    let ox = [margin_l, margin_l + panel_w + panel_gap];

    // ── Colour palette: publication-friendly on white ──────────────────────
    const PUB_COLORS: &[&str] = &[
        "#2563EB",  // blue
        "#DC2626",  // red
        "#16A34A",  // green
        "#D97706",  // amber
        "#7C3AED",  // purple
        "#0891B2",  // teal
        "#9D174D",  // rose
        "#065F46",  // dark green
    ];
    let pub_color = |k: usize| PUB_COLORS[k % PUB_COLORS.len()];

    // ── Hardcoded colours (Rust 2021: hex literals in format strings must be named) ─
    let C_111827: &str = "#111827";
    let C_374151: &str = "#374151";
    let C_6B7280: &str = "#6B7280";
    let _C_92400E: &str = "#92400E";
    let C_9CA3AF: &str = "#9CA3AF";
    let C_D1D5DB: &str = "#D1D5DB";
    let C_E5E7EB: &str = "#E5E7EB";
    let C_FAFAFA: &str = "#FAFAFA";


    // Background: white
    let mut svg = format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {total_w} {total_h}" width="{total_w}" height="{total_h}" style="background:#FFFFFF;font-family:'Helvetica Neue',Arial,sans-serif">"#
    );

    // ── Main title ─────────────────────────────────────────────────────────
    svg += &format!(
        r#"<text x="{:.1}" y="34" font-size="15" fill="{C_111827}" text-anchor="middle" font-weight="600">{} Embedding — Ambient RNA Correction  (α = {:.4})</text>"#,
        total_w as f64 / 2.0, method_name, alpha
    );

    let coords  = [before, after];
    let titles  = [
        format!("Before correction"),
        format!("After correction"),
    ];
    let sil_vals = [sil_before, sil_after];

    for p in 0..2 {
        let ox_p = ox[p];
        let data = coords[p];
        let sil  = sil_vals[p];

        let (xmin, xmax) = data_range(data, 0);
        let (ymin, ymax) = data_range(data, 1);

        let sx = |v: f64| ox_p + pad_l + data_to_screen(v, xmin, xmax, 0.0, plot_w);
        let sy = |v: f64| margin_t + pad_t + data_to_screen(v, ymax, ymin, 0.0, plot_h);

        // Panel border (light grey)
        svg += &format!(
            r#"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" fill="{C_FAFAFA}" stroke="{C_D1D5DB}" stroke-width="1" rx="4"/>"#,
            ox_p, margin_t, panel_w, panel_h
        );

        // Panel title (bold, above panel)
        let title_col = if p == 0 { "#2563EB" } else { "#16A34A" };
        svg += &format!(
            r#"<text x="{:.1}" y="{:.1}" font-size="12.5" fill="{title_col}" text-anchor="middle" font-weight="600">{}</text>"#,
            ox_p + panel_w / 2.0, margin_t - 8.0, titles[p]
        );

        // Silhouette score badge inside panel top-right
        let sil_col = if sil > 0.1 { "#16A34A" } else if sil > -0.05 { _C_92400E } else { "#DC2626" };
        svg += &format!(
            r#"<rect x="{:.1}" y="{:.1}" width="110" height="20" fill="{sil_col}" rx="3" opacity="0.12"/>"#,
            ox_p + panel_w - pad_r - 114.0, margin_t + 8.0
        );
        svg += &format!(
            r#"<text x="{:.1}" y="{:.1}" font-size="9.5" fill="{sil_col}" text-anchor="middle" font-weight="600">Silhouette = {:.4}</text>"#,
            ox_p + panel_w - pad_r - 59.0, margin_t + 19.0, sil
        );

        // Gridlines (very faint)
        let n_ticks = 4usize;
        for ti in 0..=n_ticks {
            let xv = xmin + (xmax - xmin) * ti as f64 / n_ticks as f64;
            let yv = ymin + (ymax - ymin) * ti as f64 / n_ticks as f64;
            let gx = sx(xv);
            let gy = sy(yv);
            // Vertical grid
            if gx >= ox_p + pad_l && gx <= ox_p + pad_l + plot_w {
                svg += &format!(
                    r#"<line x1="{gx:.1}" y1="{:.1}" x2="{gx:.1}" y2="{:.1}" stroke="{C_E5E7EB}" stroke-width="0.8"/>"#,
                    margin_t + pad_t, margin_t + pad_t + plot_h
                );
                // X tick + label
                svg += &format!(
                    r#"<line x1="{gx:.1}" y1="{:.1}" x2="{gx:.1}" y2="{:.1}" stroke="{C_6B7280}" stroke-width="1"/>"#,
                    margin_t + pad_t + plot_h, margin_t + pad_t + plot_h + 4.0
                );
                svg += &format!(
                    r#"<text x="{gx:.1}" y="{:.1}" font-size="8" fill="{C_6B7280}" text-anchor="middle">{:.1}</text>"#,
                    margin_t + pad_t + plot_h + 13.0, xv
                );
            }
            // Horizontal grid
            if gy >= margin_t + pad_t && gy <= margin_t + pad_t + plot_h {
                svg += &format!(
                    r#"<line x1="{:.1}" y1="{gy:.1}" x2="{:.1}" y2="{gy:.1}" stroke="{C_E5E7EB}" stroke-width="0.8"/>"#,
                    ox_p + pad_l, ox_p + pad_l + plot_w
                );
                svg += &format!(
                    r#"<line x1="{:.1}" y1="{gy:.1}" x2="{:.1}" y2="{gy:.1}" stroke="{C_6B7280}" stroke-width="1"/>"#,
                    ox_p + pad_l - 4.0, ox_p + pad_l
                );
                svg += &format!(
                    r#"<text x="{:.1}" y="{gy:.1}" font-size="8" fill="{C_6B7280}" text-anchor="end" dominant-baseline="middle">{:.1}</text>"#,
                    ox_p + pad_l - 7.0, yv
                );
            }
        }

        // Plot frame
        svg += &format!(
            r#"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" fill="none" stroke="{C_9CA3AF}" stroke-width="1.2"/>"#,
            ox_p + pad_l, margin_t + pad_t, plot_w, plot_h
        );

        // ── Per-cluster: filled ellipse background + outline + centroid + dots ──
        // Draw in this order: filled ellipse (back) → dots → outline → centroid (front)
        // This matches the souporcell paper style where clusters are visible even
        // when they overlap in t-SNE space.
        let n = data.nrows().min(labels.len());
        
        // Step 1: compute per-cluster stats (cx, cy, rx, ry, count)
        // stored as (cx, cy, rx, ry, count)
        let cluster_stats: Vec<(f64,f64,f64,f64,usize)> = (0..n_clusters).map(|k| {
            let pts: Vec<(f64,f64)> = (0..n)
                .filter(|&i| labels.get(i).cloned().unwrap_or(99) == k
                    && data[[i,0]].is_finite() && data[[i,1]].is_finite())
                .map(|i| (sx(data[[i,0]]), sy(data[[i,1]])))
                .collect();
            if pts.len() < 3 { return (0.0, 0.0, 0.0, 0.0, pts.len()); }
            let cxv = pts.iter().map(|p| p.0).sum::<f64>() / pts.len() as f64;
            let cyv = pts.iter().map(|p| p.1).sum::<f64>() / pts.len() as f64;
            let sxv = (pts.iter().map(|p| (p.0-cxv).powi(2)).sum::<f64>() / pts.len() as f64).sqrt();
            let syv = (pts.iter().map(|p| (p.1-cyv).powi(2)).sum::<f64>() / pts.len() as f64).sqrt();
            (cxv, cyv, (sxv * 1.8).max(8.0), (syv * 1.8).max(8.0), pts.len())
        }).collect();

        // Step 2: filled background ellipses (lowest layer — draw all clusters first)
        for k in 0..n_clusters {
            let s = &cluster_stats[k];
            if s.4 < 3 { continue; }
            let col = pub_color(k);
            svg += &format!(
                r#"<ellipse cx="{:.1}" cy="{:.1}" rx="{:.1}" ry="{:.1}" fill="{col}" fill-opacity="0.10" stroke="none"/>"#,
                s.0, s.1, s.2, s.3
            );
        }

        // Step 3: scatter dots — cluster by cluster, smallest cluster on top
        // Sort clusters by descending size so small clusters render on top
        let mut draw_order: Vec<usize> = (0..n_clusters).collect();
        draw_order.sort_by(|&a, &b| {
            cluster_stats[b].4.cmp(&cluster_stats[a].4)
        });
        for &k in &draw_order {
            let col = pub_color(k);
            for ni in 0..n {
                if labels.get(ni).cloned().unwrap_or(99) != k { continue; }
                let x = data[[ni, 0]];
                let y = data[[ni, 1]];
                if !x.is_finite() || !y.is_finite() { continue; }
                svg += &format!(
                    r#"<circle cx="{:.1}" cy="{:.1}" r="1.8" fill="{col}" fill-opacity="0.72"/>"#,
                    sx(x), sy(y)
                );
            }
        }

        // Step 4: ellipse outlines (on top of dots, below centroids)
        for k in 0..n_clusters {
            let s = &cluster_stats[k];
            if s.4 < 3 { continue; }
            let col = pub_color(k);
            svg += &format!(
                r#"<ellipse cx="{:.1}" cy="{:.1}" rx="{:.1}" ry="{:.1}" fill="none" stroke="{col}" stroke-width="2.0" stroke-dasharray="6,3"/>"#,
                s.0, s.1, s.2, s.3
            );
        }

        // Step 5: centroid markers (cross shape — topmost layer)
        for k in 0..n_clusters {
            let s = &cluster_stats[k];
            if s.4 < 3 { continue; }
            let col = pub_color(k);
            let arm = 7.0_f64;
            // Horizontal bar
            svg += &format!(r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{col}" stroke-width="2.5"/>"#,
                s.0 - arm, s.1, s.0 + arm, s.1);
            // Vertical bar  
            svg += &format!(r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{col}" stroke-width="2.5"/>"#,
                s.0, s.1 - arm, s.0, s.1 + arm);
            // White halo behind cross for visibility
            svg += &format!(r#"<circle cx="{:.1}" cy="{:.1}" r="3.5" fill="white" opacity="0.7"/>"#,
                s.0, s.1);
        }

        // ── Axis labels ────────────────────────────────────────────────────
        svg += &format!(
            r#"<text x="{:.1}" y="{:.1}" font-size="10" fill="{C_374151}" text-anchor="middle">{} dim 1</text>"#,
            ox_p + pad_l + plot_w / 2.0, margin_t + pad_t + plot_h + 32.0, method_name
        );
        svg += &format!(
            r#"<text x="{:.1}" y="{:.1}" font-size="10" fill="{C_374151}" text-anchor="middle" transform="rotate(-90,{:.1},{:.1})">{} dim 2</text>"#,
            ox_p + 16.0, margin_t + pad_t + plot_h / 2.0,
            ox_p + 16.0, margin_t + pad_t + plot_h / 2.0,
            method_name
        );
    }

    // ── Legend (inside right panel, bottom-left) ───────────────────────────
    let leg_ox = ox[1] + pad_l + 8.0;
    let leg_oy = margin_t + pad_t + plot_h - (n_clusters.min(8) as f64) * 17.0 - 12.0;
    svg += &format!(
        r#"<rect x="{:.1}" y="{:.1}" width="90" height="{:.1}" fill="rgba(255,255,255,0.85)" stroke="{C_D1D5DB}" stroke-width="0.8" rx="3"/>"#,
        leg_ox - 4.0, leg_oy - 6.0, n_clusters.min(8) as f64 * 17.0 + 12.0
    );
    for k in 0..n_clusters.min(8) {
        let ly = leg_oy + k as f64 * 17.0;
        let col = pub_color(k);
        svg += &format!(r#"<circle cx="{:.1}" cy="{:.1}" r="5" fill="{col}"/>"#,
            leg_ox + 5.0, ly);
        svg += &format!(
            r#"<text x="{:.1}" y="{:.1}" font-size="9.5" fill="{C_374151}" dominant-baseline="middle">Cluster {k}</text>"#,
            leg_ox + 14.0, ly
        );
    }

    // ── Footer annotation ─────────────────────────────────────────────────
    let delta = sil_after - sil_before;
    // Use neutral grey for |delta| < 0.005, green for positive, red for negative
    let (dcol, dtext) = if delta.abs() < 0.005 {
        (C_6B7280, format!("Silhouette Δ ≈ 0  (before = {:.4}  after = {:.4})  — Expected: α = {:.4} is small; use --simulate to test at higher contamination levels", sil_before, sil_after, alpha))
    } else if delta > 0.0 {
        ("#16A34A", format!("Silhouette Δ = {:+.4}  (before = {:.4}  after = {:.4})  — Correction improved cluster separation", delta, sil_before, sil_after))
    } else {
        ("#DC2626", format!("Silhouette Δ = {:+.4}  (before = {:.4}  after = {:.4})", delta, sil_before, sil_after))
    };
    svg += &format!(
        r#"<text x="{:.1}" y="{:.1}" font-size="9" fill="{dcol}" text-anchor="middle">{dtext}</text>"#,
        total_w as f64 / 2.0, total_h as f64 - 10.0
    );

    svg += "</svg>";
    svg
}


// ── PLOT 3: Silhouette comparison bar chart ───────────────────────────────────

pub fn plot_silhouette_comparison(
    per_cluster_before: &[f64],
    per_cluster_after:  &[f64],
    global_before:       f64,
    global_after:        f64,
    n_clusters:          usize,
) -> String {
    let w = 800usize;
    let h = 520usize;
    let pad_l = 80.0_f64;
    let pad_r = 160.0_f64;   // extra right margin for global labels
    let pad_t = 60.0_f64;
    let pad_b = 80.0_f64;
    let plot_w = w as f64 - pad_l - pad_r;
    let plot_h = h as f64 - pad_t - pad_b;

    // ── Y-axis: span ACTUAL value range, not forced 0→0.5 ─────────────────
    // Silhouette scores can be negative. Build axis around actual data range.
    let all_vals: Vec<f64> = per_cluster_before.iter()
        .chain(per_cluster_after.iter())
        .chain(std::iter::once(&global_before))
        .chain(std::iter::once(&global_after))
        .cloned().collect();
    let v_min = all_vals.iter().cloned().fold(f64::INFINITY, f64::min);
    let v_max = all_vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    // Pad by 30% of range, always include 0
    let range = (v_max - v_min).max(0.01);
    let y_min = (v_min - range * 0.30).min(-0.001);
    let y_max = (v_max + range * 0.30).max( 0.001);
    let y_range = y_max - y_min;

    let to_y = |v: f64| -> f64 {
        pad_t + plot_h * (1.0 - (v - y_min) / y_range)
    };
    let zero_y = to_y(0.0);

    let n = n_clusters.min(per_cluster_before.len()).min(per_cluster_after.len());
    let group_w = plot_w / (n + 1) as f64;
    let bar_w   = group_w * 0.35;

    let mut svg = svg_open(w, h);

    svg += &svg_text(w as f64 / 2.0, 28.0,
        "Per-Cluster Silhouette Score — Before vs After Correction",
        13.0, INK, "middle");

    svg += &svg_rect(pad_l, pad_t, plot_w, plot_h, SURFACE, 6.0);

    // Y gridlines at nice intervals
    let n_ticks = 6usize;
    for ti in 0..=n_ticks {
        let v = y_min + y_range * ti as f64 / n_ticks as f64;
        let y = to_y(v);
        if y < pad_t || y > pad_t + plot_h { continue; }
        let lw = if v.abs() < 1e-6 { 1.5 } else { 0.8 };
        let lc = if v.abs() < 1e-6 { INK2 } else { BORDER };
        svg += &svg_line(pad_l, y, pad_l + plot_w, y, lc, lw);
        svg += &svg_text(pad_l - 6.0, y, &format!("{:.3}", v), 8.5, INK3, "end");
    }

    // Zero baseline (prominent)
    svg += &svg_line(pad_l, zero_y, pad_l + plot_w, zero_y, INK2, 1.5);
    svg += &svg_text(pad_l - 6.0, zero_y, "0.000", 8.5, INK2, "end");

    // Bars — support both positive (above zero) and negative (below zero)
    for k in 0..n {
        let cx = pad_l + (k + 1) as f64 * group_w;
        let bv = per_cluster_before[k];
        let av = per_cluster_after[k];

        for (offset, val, is_after) in [(-bar_w/2.0 - 1.0, bv, false), (1.0, av, true)] {
            let bar_top = to_y(val.max(0.0_f64));
            let bar_bot = to_y(val.min(0.0_f64));
            let bar_h   = (bar_bot - bar_top).max(2.0);
            let col_base = cluster_colour(k);
            let fill = if is_after {
                col_base.to_string()
            } else {
                format!("{}88", col_base)
            };
            let x = cx + offset;
            svg += &svg_rect(x, bar_top, bar_w / 2.0, bar_h, &fill, 2.0);

            // Value label above/below bar
            let label_y = if val >= 0.0 { bar_top - 9.0 } else { bar_bot + 12.0 };
            let ink = if is_after { INK } else { INK2 };
            svg += &svg_text(x + bar_w / 4.0, label_y,
                &format!("{:.3}", val), 7.5, ink, "middle");
        }

        // Delta annotation below x-axis
        let delta = av - bv;
        let dcol = if delta.abs() < 0.001 { INK3 }
                   else if delta > 0.0 { GREEN } else { RED };
        let x_lbl = cx - bar_w / 4.0;
        svg += &svg_text(x_lbl, pad_t + plot_h + 28.0,
            &format!("{:+.4}", delta), 8.0, dcol, "middle");
        svg += &svg_text(x_lbl, pad_t + plot_h + 16.0,
            &format!("Cluster {}", k), 9.0, INK3, "middle");
    }

    // Global silhouette dashed reference lines
    for (gval, gcol, glabel) in [
        (global_before, CYAN,  format!("global before = {:.4}", global_before)),
        (global_after,  GREEN, format!("global after  = {:.4}", global_after)),
    ] {
        let gy = to_y(gval);
        if gy >= pad_t && gy <= pad_t + plot_h {
            svg += &format!(
                r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{gcol}" stroke-width="1.5" stroke-dasharray="6,3" opacity="0.9"/>"#,
                pad_l, gy, pad_l + plot_w, gy);
            svg += &svg_text(pad_l + plot_w + 6.0, gy, &glabel, 8.0, gcol, "start");
        }
    }

    // Legend
    svg += &svg_rect(pad_l + 4.0, pad_t + 4.0, 150.0, 36.0, BG, 4.0);
    svg += &svg_rect(pad_l + 10.0, pad_t + 10.0, 12.0, 10.0, &format!("{}88", CYAN), 2.0);
    svg += &svg_text(pad_l + 26.0, pad_t + 15.0, "Before correction", 8.5, INK2, "start");
    svg += &svg_rect(pad_l + 10.0, pad_t + 24.0, 12.0, 10.0, CYAN, 2.0);
    svg += &svg_text(pad_l + 26.0, pad_t + 29.0, "After correction",  8.5, INK,  "start");

    // Interpretation note at bottom
    let global_delta = global_after - global_before;
    let note = if global_delta.abs() < 0.001 {
        format!("Global Δ = {:+.4}  →  No change expected at α=1.26% (correct result — use simulation to confirm at higher α)", global_delta)
    } else if global_delta > 0.0 {
        format!("Global Δ = {:+.4}  →  Correction improved cluster separation", global_delta)
    } else {
        format!("Global Δ = {:+.4}  →  Correction decreased cluster separation", global_delta)
    };
    svg += &svg_text(w as f64 / 2.0, h as f64 - 14.0, &note, 8.0, INK3, "middle");

    svg += &format!(
        r#"<text x="14" y="{:.1}" font-size="9" fill="{INK3}" text-anchor="middle" transform="rotate(-90,14,{:.1})">Silhouette Score</text>"#,
        pad_t + plot_h / 2.0, pad_t + plot_h / 2.0);

    svg += svg_close();
    svg
}

// ── PLOT 4: Correction magnitude distribution ─────────────────────────────────

/// Histogram of per-cell mean correction amounts.
/// `corrections` = per-cell mean |X_raw - X_corrected| (non-NaN entries only)
pub fn plot_correction_distribution(
    corrections: &[f64],
    alpha:        f64,
    n_floored:    usize,
    pct_floored:  f64,
) -> String {
    let w = 720usize;
    let h = 400usize;
    let pad_l = 70.0_f64;
    let pad_r = 30.0_f64;
    let pad_t = 50.0_f64;
    let pad_b = 60.0_f64;
    let plot_w = w as f64 - pad_l - pad_r;
    let plot_h = h as f64 - pad_t - pad_b;
    let n_bins  = 40usize;

    let vals: Vec<f64> = corrections.iter().cloned().filter(|v| v.is_finite()).collect();
    if vals.is_empty() {
        let mut svg = svg_open(w, h);
        svg += &svg_text(w as f64/2.0, h as f64/2.0, "No data", 12.0, INK3, "middle");
        svg += svg_close();
        return svg;
    }

    let vmin = vals.iter().cloned().fold(f64::INFINITY, f64::min);
    let vmax = vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let bin_w = (vmax - vmin) / n_bins as f64;

    let mut counts = vec![0usize; n_bins];
    for &v in &vals {
        let b = ((v - vmin) / (bin_w + 1e-12)) as usize;
        let b = b.min(n_bins - 1);
        counts[b] += 1;
    }
    let count_max = *counts.iter().max().unwrap_or(&1) as f64;

    let to_x = |b: usize| -> f64 { pad_l + b as f64 / n_bins as f64 * plot_w };
    let to_y = |c: usize| -> f64 { pad_t + plot_h - c as f64 / count_max * plot_h };

    let mut svg = svg_open(w, h);
    svg += &svg_text(w as f64/2.0, 26.0,
        &format!("Correction Magnitude Distribution  α={:.6}", alpha),
        12.0, INK, "middle");
    svg += &svg_rect(pad_l, pad_t, plot_w, plot_h, SURFACE, 6.0);

    // Bars
    let bw = plot_w / n_bins as f64;
    for (b, &c) in counts.iter().enumerate() {
        if c == 0 { continue; }
        let x = to_x(b);
        let y = to_y(c);
        let bar_h = pad_t + plot_h - y;
        svg += &svg_rect(x + 1.0, y, bw - 2.0, bar_h, VIOLET, 2.0);
    }

    // Axes
    svg += &svg_line(pad_l, pad_t + plot_h, pad_l + plot_w, pad_t + plot_h, INK3, 1.0);
    svg += &svg_line(pad_l, pad_t, pad_l, pad_t + plot_h, INK3, 1.0);

    // X ticks
    for ti in 0..=5 {
        let v = vmin + (vmax - vmin) * ti as f64 / 5.0;
        let x = pad_l + (v - vmin) / (vmax - vmin + 1e-12) * plot_w;
        svg += &svg_line(x, pad_t + plot_h, x, pad_t + plot_h + 4.0, INK3, 1.0);
        svg += &svg_text(x, pad_t + plot_h + 14.0,
            &format!("{:.4}", v), 8.0, INK3, "middle");
    }

    svg += &svg_text(pad_l + plot_w / 2.0, h as f64 - 8.0,
        "Mean correction per cell", 9.0, INK3, "middle");
    svg += &format!(
        r#"<text x="14" y="{:.1}" font-size="9" fill="{}" text-anchor="middle" transform="rotate(-90,14,{:.1})">Cell count</text>"#,
        pad_t + plot_h / 2.0, INK3, pad_t + plot_h / 2.0);

    // Stats annotation
    let mean = vals.iter().sum::<f64>() / vals.len() as f64;
    svg += &svg_rect(pad_l + 4.0, pad_t + 4.0, 220.0, 54.0, BG, 4.0);
    svg += &svg_text(pad_l + 10.0, pad_t + 14.0,
        &format!("n cells:      {}", vals.len()), 8.5, INK2, "start");
    svg += &svg_text(pad_l + 10.0, pad_t + 26.0,
        &format!("mean corr:    {:.6}", mean), 8.5, INK2, "start");
    svg += &svg_text(pad_l + 10.0, pad_t + 38.0,
        &format!("floored:      {} ({:.1}%)", n_floored, pct_floored * 100.0),
        8.5, if pct_floored > 0.05 { AMBER } else { INK2 }, "start");
    svg += &svg_text(pad_l + 10.0, pad_t + 50.0,
        &format!("alpha used:   {:.8}", alpha), 8.5, CYAN, "start");

    svg += svg_close();
    svg
}

// ── PLOT 5: Simulation curves ─────────────────────────────────────────────────

pub fn plot_simulation_curves(result: &BenchmarkResult) -> String {
    if result.points.is_empty() { return String::new(); }

    let w = 900usize;
    let h = 520usize;
    let pad_l = 72.0_f64;
    let pad_r  = 80.0_f64;   // second y-axis
    let pad_t  = 56.0_f64;
    let pad_b  = 70.0_f64;
    let plot_w = w as f64 - pad_l - pad_r;
    let plot_h = h as f64 - pad_t - pad_b;

    let pts = &result.points;
    let alpha_vals: Vec<f64> = pts.iter().map(|p| p.alpha_sim).collect();
    let rmse_vals:  Vec<f64> = pts.iter().map(|p| p.rmse).collect();
    let sil_bef:    Vec<f64> = pts.iter().map(|p| p.silhouette_before).collect();
    let sil_aft:    Vec<f64> = pts.iter().map(|p| p.silhouette_after).collect();

    let alpha_min = alpha_vals.first().cloned().unwrap_or(0.0);
    let alpha_max = alpha_vals.last().cloned().unwrap_or(1.0);
    let rmse_max  = rmse_vals.iter().cloned().fold(0.0_f64, f64::max) * 1.15;
    let rmse_max  = if rmse_max < 1e-10 { 1.0 } else { rmse_max };
    let sil_min   = sil_bef.iter().chain(sil_aft.iter()).cloned().fold(0.0_f64, f64::min) - 0.05;
    let sil_max   = sil_bef.iter().chain(sil_aft.iter()).cloned().fold(0.0_f64, f64::max) + 0.05;

    let to_x = |a: f64| pad_l + (a - alpha_min) / (alpha_max - alpha_min + 1e-10) * plot_w;
    let to_y_rmse = |v: f64| pad_t + plot_h - v / rmse_max * plot_h;
    let to_y_sil  = |v: f64| pad_t + plot_h - (v - sil_min) / (sil_max - sil_min + 1e-10) * plot_h;

    let mut svg = svg_open(w, h);

    // Title
    svg += &svg_text(w as f64 / 2.0, 28.0,
        "Simulation Benchmark — Correction Recovery vs Contamination Level",
        13.0, INK, "middle");

    svg += &svg_rect(pad_l, pad_t, plot_w, plot_h, SURFACE, 6.0);

    // Gridlines
    for ti in 0..=5 {
        let y = pad_t + plot_h * ti as f64 / 5.0;
        svg += &svg_line(pad_l, y, pad_l + plot_w, y, BORDER, 1.0);
    }

    // True alpha vertical marker
    let tx = to_x(result.true_alpha);
    svg += &format!(
        r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{}" stroke-width="1.5" stroke-dasharray="4,3" opacity="0.8"/>"#,
        tx, pad_t, tx, pad_t + plot_h, AMBER);
    svg += &svg_text(tx + 3.0, pad_t + 10.0,
        &format!("true α={:.4}", result.true_alpha), 8.0, AMBER, "start");

    // RMSE line (primary Y axis, left)
    let rmse_path: String = pts.iter().enumerate().map(|(i, p)| {
        let x = to_x(p.alpha_sim);
        let y = to_y_rmse(p.rmse);
        if i == 0 { format!("M {:.1} {:.1}", x, y) }
        else       { format!("L {:.1} {:.1}", x, y) }
    }).collect::<Vec<_>>().join(" ");
    svg += &format!(r#"<path d="{}" fill="none" stroke="{}" stroke-width="2.5"/>"#, rmse_path, CYAN);
    for p in pts {
        svg += &format!(r#"<circle cx="{:.1}" cy="{:.1}" r="4" fill="{}"/>"#,
            to_x(p.alpha_sim), to_y_rmse(p.rmse), CYAN);
    }

    // Silhouette BEFORE line
    let sil_b_path: String = pts.iter().enumerate().map(|(i, p)| {
        let x = to_x(p.alpha_sim); let y = to_y_sil(p.silhouette_before);
        if i == 0 { format!("M {:.1} {:.1}", x, y) } else { format!("L {:.1} {:.1}", x, y) }
    }).collect::<Vec<_>>().join(" ");
    svg += &format!(r#"<path d="{}" fill="none" stroke="{}" stroke-width="2" stroke-dasharray="6,3"/>"#,
        sil_b_path, AMBER);

    // Silhouette AFTER line
    let sil_a_path: String = pts.iter().enumerate().map(|(i, p)| {
        let x = to_x(p.alpha_sim); let y = to_y_sil(p.silhouette_after);
        if i == 0 { format!("M {:.1} {:.1}", x, y) } else { format!("L {:.1} {:.1}", x, y) }
    }).collect::<Vec<_>>().join(" ");
    svg += &format!(r#"<path d="{}" fill="none" stroke="{}" stroke-width="2.5"/>"#,
        sil_a_path, GREEN);
    for p in pts {
        svg += &format!(r#"<circle cx="{:.1}" cy="{:.1}" r="3" fill="{}"/>"#,
            to_x(p.alpha_sim), to_y_sil(p.silhouette_after), GREEN);
    }

    // X axis labels
    for p in pts {
        svg += &svg_text(to_x(p.alpha_sim), pad_t + plot_h + 16.0,
            &format!("{:.2}", p.alpha_sim), 8.5, INK3, "middle");
    }
    svg += &svg_text(pad_l + plot_w / 2.0, h as f64 - 8.0,
        "Simulated ambient fraction (α)", 9.0, INK3, "middle");

    // Left Y axis label (RMSE)
    svg += &format!(
        r#"<text x="12" y="{:.1}" font-size="9" fill="{}" text-anchor="middle" transform="rotate(-90,12,{:.1})">RMSE (recovery)</text>"#,
        pad_t + plot_h / 2.0, CYAN, pad_t + plot_h / 2.0);

    // Right Y axis label (silhouette)
    let rx = w as f64 - 14.0;
    svg += &format!(
        r#"<text x="{:.1}" y="{:.1}" font-size="9" fill="{}" text-anchor="middle" transform="rotate(90,{:.1},{:.1})">Silhouette score</text>"#,
        rx, pad_t + plot_h / 2.0, GREEN, rx, pad_t + plot_h / 2.0);

    // Legend
    svg += &svg_rect(pad_l + 4.0, pad_t + 4.0, 200.0, 54.0, BG, 4.0);
    svg += &svg_line(pad_l + 10.0, pad_t + 14.0, pad_l + 28.0, pad_t + 14.0, CYAN, 2.5);
    svg += &svg_text(pad_l + 34.0, pad_t + 14.0, "RMSE (recovery quality)", 8.5, INK2, "start");
    svg += &format!(r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{}" stroke-width="2" stroke-dasharray="6,3"/>"#,
        pad_l + 10.0, pad_t + 28.0, pad_l + 28.0, pad_t + 28.0, AMBER);
    svg += &svg_text(pad_l + 34.0, pad_t + 28.0, "Silhouette (before)", 8.5, INK2, "start");
    svg += &svg_line(pad_l + 10.0, pad_t + 42.0, pad_l + 28.0, pad_t + 42.0, GREEN, 2.5);
    svg += &svg_text(pad_l + 34.0, pad_t + 42.0, "Silhouette (after)", 8.5, INK, "start");

    svg += svg_close();
    svg
}

// ── Write helpers ─────────────────────────────────────────────────────────────

pub fn write_svg(svg: &str, path: &str) -> Result<(), String> {
    let p = std::path::Path::new(path);
    if let Some(parent) = p.parent() {
        std::fs::create_dir_all(parent).ok();
    }
    std::fs::write(path, svg).map_err(|e| format!("Cannot write {}: {}", path, e))
}

// ══════════════════════════════════════════════════════════════════════════════
// PLOT 6: Soup vector profile
// ══════════════════════════════════════════════════════════════════════════════
//
// Two-panel figure:
//   Top:    Horizontal bar chart — top 30 loci by soup[l] magnitude (sorted desc)
//   Bottom: Histogram of ALL L loci's soup[l] values (distribution shape)
//
// Data: soup_contribution[l] (len = L_raw), locus_mean[l], alpha
// No structural changes needed — all data is in SoupVector.
// ─────────────────────────────────────────────────────────────────────────────

pub fn plot_soup_vector_profile(
    soup_contribution: &[f64],
    locus_mean:        &[f64],
    alpha:             f64,
    top_n:             usize,   // typically 30
) -> String {
    let top_n = top_n.min(soup_contribution.len());
    let total_w = 900usize;
    let total_h = 700usize;

    // ── Sort loci by locus_mean DESC (μ[l] = mean allele fraction per locus)
    // soup[l] = α × μ[l] — all loci with μ=1.0 tie at the same soup value,
    // making a bar chart sorted by soup uninformative (all bars equal length).
    // Sorting by μ[l] gives meaningful variation: loci range from μ≈0 to μ=1.
    // These high-μ loci are the most alt-biased and thus most soup-contaminated.
    let mut indexed: Vec<(usize, f64)> = locus_mean.iter()
        .cloned().enumerate()
        .filter(|(_, mu)| mu.is_finite() && *mu > 0.0)
        .collect();
    indexed.sort_by(|a, b| b.1.total_cmp(&a.1));

    // Deduplicate ties at μ=1.0: keep only the first 5, then skip to varied values
    // This prevents 30 identical bars and shows the actual distribution shape.
    let mut top: Vec<(usize, f64)> = Vec::with_capacity(top_n);
    let mut n_max = 0usize;
    for &(idx, mu) in &indexed {
        if mu >= 0.9999 {
            if n_max < 3 { top.push((idx, mu)); n_max += 1; }
        } else {
            top.push((idx, mu));
            if top.len() >= top_n { break; }
        }
    }
    // Fill remaining if needed
    for &(idx, mu) in &indexed {
        if top.len() >= top_n { break; }
        if !top.iter().any(|&(i, _)| i == idx) {
            top.push((idx, mu));
        }
    }
    let soup_vals: Vec<f64> = soup_contribution.iter().cloned()
        .filter(|v| v.is_finite()).collect();

    // Panel dimensions
    let bar_area_h = 440usize;
    let hist_area_y = bar_area_h + 40;
    let hist_area_h = total_h - hist_area_y - 30;
    let pad_l = 80.0_f64;
    let pad_r = 30.0_f64;
    let pad_t = 50.0_f64;
    let plot_w = total_w as f64 - pad_l - pad_r;

    // Bar scale: use μ[l] (locus mean) as the axis — ranges 0→1
    // This gives meaningful variation since μ[l] ∈ [0,1]
    let max_mu = top.first().map(|(_,v)| *v).unwrap_or(1.0).max(1e-12);
    let bar_plot_h = bar_area_h as f64 - pad_t - 20.0;
    let bar_h_each = (bar_plot_h / top_n as f64).min(18.0).max(6.0);

    let mut svg = svg_open(total_w, total_h);

    // ── PANEL 1 TITLE ──────────────────────────────────────────────────────
    svg += &svg_text(total_w as f64 / 2.0, 26.0,
        &format!("Soup Vector Profile — Top {} Loci by μ[l]  (α = {:.8}, L = {})",
                 top_n, alpha, soup_contribution.len()),
        13.0, INK, "middle");

    // Background
    svg += &svg_rect(pad_l, pad_t, plot_w, bar_plot_h + 10.0, SURFACE, 6.0);

    // Axis label
    svg += &svg_text(pad_l + plot_w / 2.0, pad_t + bar_plot_h + 26.0,
        "μ[l] = mean allele fraction per locus  (soup[l] = α × μ[l])", 9.0, INK3, "middle");
    svg += &svg_text(12.0, pad_t + bar_plot_h / 2.0,
        &format!("Top {} loci\n(by μ[l])", top_n), 9.0, INK3, "middle");

    // Bars
    for (rank, &(locus_idx, soup_val)) in top.iter().enumerate() {
        let y   = pad_t + rank as f64 * bar_h_each + 4.0;
        // mu IS the bar value now (top contains (idx, mu))
        let mu_l   = soup_val;  // top is sorted by mu, so soup_val holds mu
        let soup_l = mu_l * alpha; // soup[l] = α × μ[l]
        let bw   = (mu_l / max_mu * plot_w).max(2.0);

        // Colour: cyan=highest μ, violet=mid, muted=lower
        let intensity = mu_l / max_mu;
        let col = if intensity > 0.66 { CYAN }
                  else if intensity > 0.33 { VIOLET }
                  else { INK3 };

        svg += &svg_rect(pad_l, y, bw, bar_h_each - 1.0, col, 2.0);

        // Value label (only if bar wide enough)
        if bw > 60.0 {
            svg += &format!(
                r#"<text x="{:.1}" y="{:.1}" font-size="8" fill="{BG}"                    text-anchor="end" dominant-baseline="middle">μ={:.4}  soup={:.6}</text>"#,
                pad_l + bw - 4.0, y + bar_h_each / 2.0 - 0.5,
                mu_l, soup_l
            );
        } else {
            svg += &format!(
                r#"<text x="{:.1}" y="{:.1}" font-size="7.5" fill="{INK3}" \
                   text-anchor="start" dominant-baseline="middle">μ={:.4}</text>"#,
                pad_l + bw + 4.0, y + bar_h_each / 2.0 - 0.5,
                mu_l
            );
        }

        // Locus index label
        svg += &format!(
            r#"<text x="{:.1}" y="{:.1}" font-size="7.5" fill="{INK3}" \
               text-anchor="end" dominant-baseline="middle">l={}</text>"#,
            pad_l - 4.0, y + bar_h_each / 2.0 - 0.5, locus_idx
        );
    }

    // X-axis ticks: μ[l] scale
    for ti in 0..=4 {
        let v = max_mu * ti as f64 / 4.0;
        let x = pad_l + v / max_mu * plot_w;
        svg += &svg_line(x, pad_t + bar_plot_h + 10.0, x, pad_t + bar_plot_h + 16.0, INK3, 1.0);
        svg += &svg_text(x, pad_t + bar_plot_h + 24.0,
            &format!("μ={:.2}", v), 8.0, INK3, "middle");
        // Also show corresponding soup value
        svg += &svg_text(x, pad_t + bar_plot_h + 34.0,
            &format!("s={:.5}", v * alpha), 7.5, INK3, "middle");
    }

    // ── PANEL 2: Histogram of ALL soup values ──────────────────────────────
    let n_bins = 50usize;
    let hist_y = hist_area_y as f64;
    let hist_h = hist_area_h as f64 - 20.0;

    svg += &svg_rect(pad_l, hist_y, plot_w, hist_h, SURFACE, 6.0);
    svg += &svg_text(total_w as f64 / 2.0, hist_y - 10.0,
        &format!("Distribution of all {} loci's soup values", soup_vals.len()),
        10.0, INK2, "middle");

    if !soup_vals.is_empty() {
        let smin = soup_vals.iter().cloned().fold(f64::INFINITY, f64::min);
        let smax = soup_vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let bw_data = (smax - smin) / n_bins as f64;
        let mut counts = vec![0usize; n_bins];
        for &v in &soup_vals {
            let b = ((v - smin) / (bw_data + 1e-14)) as usize;
            counts[b.min(n_bins - 1)] += 1;
        }
        // Log-scale histogram (soup values are extremely right-skewed:
        // most loci have tiny soup, a few have values near alpha)
        let log_max = counts.iter().map(|&c| ((c + 1) as f64).log10())
            .fold(0.0_f64, f64::max).max(0.1);
        let bw_px = plot_w / n_bins as f64;

        for (b, &c) in counts.iter().enumerate() {
            if c == 0 { continue; }
            let x     = pad_l + b as f64 * bw_px;
            let log_c = ((c + 1) as f64).log10();
            let bh    = (log_c / log_max * hist_h).max(1.0);
            let y     = hist_y + hist_h - bh;
            // Colour: first bin (loci with near-zero soup) in cyan; rest in violet
            let col = if b == 0 { CYAN } else { VIOLET };
            svg += &svg_rect(x + 0.5, y, bw_px - 1.0, bh, col, 1.0);
        }

        // Log Y-axis label
        let max_pow = log_max.ceil() as usize;
        for p in 0..=max_pow {
            let y = hist_y + hist_h - (p as f64 / log_max * hist_h);
            svg += &svg_line(pad_l - 4.0, y, pad_l, y, INK3, 1.0);
            svg += &svg_text(pad_l - 6.0, y, &format!("10^{}", p), 7.0, INK3, "end");
        }

        // Median and mean lines
        let mut sv2 = soup_vals.clone();
        sv2.sort_by(|a, b| a.total_cmp(b));
        let median   = sv2[sv2.len() / 2];
        let mean_val = sv2.iter().sum::<f64>() / sv2.len() as f64;
        for (val, col, label) in [
            (median,   AMBER, format!("median={:.6}", median)),
            (mean_val, CYAN,  format!("mean={:.6}",   mean_val)),
        ] {
            let mx = pad_l + (val - smin) / (smax - smin + 1e-14) * plot_w;
            svg += &format!(
                r#"<line x1="{mx:.1}" y1="{hist_y:.1}" x2="{mx:.1}" y2="{:.1}" stroke="{col}" stroke-width="1.5" stroke-dasharray="4,3"/>"#,
                hist_y + hist_h);
            svg += &svg_text(mx + 3.0, hist_y + 10.0, &label, 8.0, col, "start");
        }

        // X axis ticks
        for ti in 0..=5 {
            let v = smin + (smax - smin) * ti as f64 / 5.0;
            let x = pad_l + (v - smin) / (smax - smin + 1e-14) * plot_w;
            svg += &svg_text(x, hist_y + hist_h + 14.0, &format!("{:.6}", v), 7.5, INK3, "middle");
        }
        svg += &svg_text(pad_l + plot_w / 2.0, hist_y + hist_h + 26.0,
            "soup[l] = α × μ[l]  (log₁₀ count y-axis)", 9.0, INK3, "middle");

        // Annotation
        let near_zero = sv2.iter().filter(|&&v| v < smax * 0.05).count();
        svg += &svg_text(pad_l + 6.0, hist_y + 14.0,
            &format!("{} loci with soup < 5% of max  ({:.1}%)",
                near_zero, near_zero as f64 / sv2.len().max(1) as f64 * 100.0),
            7.5, INK2, "start");
    }

    svg += svg_close();
    svg
}

// ══════════════════════════════════════════════════════════════════════════════
// PLOT 7: Per-cluster correction strip
// ══════════════════════════════════════════════════════════════════════════════
//
// One row per cluster showing 4 metrics as proportional segments or bars:
//   • Cell count bar       (relative width = n_cells / max_cells)
//   • Mean correction      (violet bar, scale 0→max_correction)
//   • Floor%               (amber segment)
//   • Silhouette Δ         (green/red dot with numeric label)
//
// Data: labels[], cell_corrections[], per_cluster_silhouette before+after,
//       cluster_counts[], n_clusters.  All in scope in main.rs.
// ─────────────────────────────────────────────────────────────────────────────

pub fn plot_per_cluster_correction(
    labels:            &[usize],
    cell_corrections:  &[f64],
    sil_before:        &[f64],
    sil_after:         &[f64],
    cluster_counts:    &[usize],
    n_clusters:        usize,
    alpha:             f64,
    x_raw:             &ndarray::Array2<f64>,   // [L_raw × N] — for real floor count
    soup_contribution: &[f64],                   // [L_raw]     — soup[l] = α×μ[l]
    cell_col_indices:  &[usize],                 // working-cell → matrix column mapping
) -> String {
    let w = 900usize;
    let row_h = 56usize;
    let pad_t = 70usize;
    let pad_l = 110usize;
    let pad_r = 120usize;
    let pad_b = 50usize;
    let h = pad_t + n_clusters * row_h + pad_b;
    let plot_w = w - pad_l - pad_r;

    let max_cells = cluster_counts.iter().max().cloned().unwrap_or(1).max(1);

    // ── Per-cluster mean correction (from labels × cell_corrections) ──────
    let mut per_cluster_corr  = vec![0.0f64; n_clusters];
    let mut per_cluster_count = vec![0usize;  n_clusters];

    for (&lbl, &corr) in labels.iter().zip(cell_corrections.iter()) {
        if lbl < n_clusters {
            per_cluster_corr[lbl]  += corr;
            per_cluster_count[lbl] += 1;
        }
    }
    for k in 0..n_clusters {
        if per_cluster_count[k] > 0 {
            per_cluster_corr[k] /= per_cluster_count[k] as f64;
        }
    }
    let max_corr = per_cluster_corr.iter().cloned()
        .fold(0.0_f64, f64::max).max(1e-10);

    // ── Per-cluster REAL floor fraction ───────────────────────────────────
    // A cell is "floored at a locus" when x_raw[l, col] < soup[l] (and not NaN).
    // floor_pct[k] = (total floored entries for cluster k) / (total covered entries for k)
    // This uses the actual x_raw matrix and soup_contribution vector — no approximation.
    let mut per_cluster_floored = vec![0usize; n_clusters];
    let mut per_cluster_covered = vec![0usize; n_clusters];

    let n_loci = x_raw.nrows().min(soup_contribution.len());
    for (wi, (&col_idx, &lbl)) in cell_col_indices.iter().zip(labels.iter()).enumerate() {
        if lbl >= n_clusters || col_idx >= x_raw.ncols() { continue; }
        let _ = wi;
        for li in 0..n_loci {
            let v = x_raw[[li, col_idx]];
            if !v.is_nan() {
                per_cluster_covered[lbl] += 1;
                if v < soup_contribution[li] {
                    per_cluster_floored[lbl] += 1;
                }
            }
        }
    }
    let per_cluster_floor_pct: Vec<f64> = (0..n_clusters).map(|k| {
        if per_cluster_covered[k] > 0 {
            per_cluster_floored[k] as f64 / per_cluster_covered[k] as f64
        } else { 0.0 }
    }).collect();
    let max_floor_pct = per_cluster_floor_pct.iter().cloned()
        .fold(0.0_f64, f64::max).max(1e-10);

    let mut svg = svg_open(w, h);

    // Title
    svg += &svg_text(w as f64 / 2.0, 26.0,
        &format!("Per-Cluster Correction Summary  (α = {:.6}, {} clusters)",
                 alpha, n_clusters),
        13.0, INK, "middle");

    // Column headers
    let headers = [
        (pad_l as f64 + plot_w as f64 * 0.0 / 4.0, "Cells"),
        (pad_l as f64 + plot_w as f64 * 1.0 / 4.0, "Mean |Correction|"),
        (pad_l as f64 + plot_w as f64 * 2.0 / 4.0, "Floor %"),
        (pad_l as f64 + plot_w as f64 * 3.0 / 4.0, "Silhouette Δ"),
    ];
    for &(x, label) in &headers {
        svg += &svg_text(x + plot_w as f64 / 8.0, 50.0, label, 9.0, INK3, "middle");
    }
    let col_w = plot_w as f64 / 4.0;

    for k in 0..n_clusters {
        let y      = pad_t + k * row_h;
        let col    = cluster_colour(k);
        let n_cells = cluster_counts.get(k).cloned().unwrap_or(0);
        let corr   = per_cluster_corr[k];
        let delta_sil = if k < sil_after.len() && k < sil_before.len() {
            sil_after[k] - sil_before[k]
        } else { 0.0 };

        // Row background
        svg += &svg_rect(pad_l as f64, y as f64 + 2.0,
                          plot_w as f64, row_h as f64 - 4.0,
                          SURFACE, 4.0);

        // Cluster label
        svg += &svg_text(pad_l as f64 - 8.0, y as f64 + row_h as f64 / 2.0,
            &format!("Cluster {k}"), 10.0, col, "end");

        let bar_y = y as f64 + (row_h as f64 - 16.0) / 2.0;

        // Column 0: Cell count bar
        let cx = pad_l as f64;
        let cw  = n_cells as f64 / max_cells as f64 * (col_w - 12.0);
        svg += &svg_rect(cx + 4.0, bar_y, cw, 16.0, col, 3.0);
        svg += &svg_text(cx + cw + 8.0, y as f64 + row_h as f64 / 2.0,
            &n_cells.to_string(), 9.0, INK2, "start");

        // Column 1: Mean correction bar
        let mx = pad_l as f64 + col_w;
        let mw  = corr / max_corr * (col_w - 12.0);
        svg += &svg_rect(mx + 4.0, bar_y, mw.max(2.0), 16.0, VIOLET, 3.0);
        svg += &svg_text(mx + mw + 8.0, y as f64 + row_h as f64 / 2.0,
            &format!("{:.5}", corr), 8.5, INK2, "start");

        // Column 2: Real floor fraction for this cluster
        // Computed above from x_raw[l,col] < soup[l] over all covered (non-NaN) entries
        let floor_pct = per_cluster_floor_pct[k];
        let fx = pad_l as f64 + col_w * 2.0;
        let fw  = if max_floor_pct > 0.0 { floor_pct / max_floor_pct * (col_w - 12.0) } else { 0.0 };
        svg += &svg_rect(fx + 4.0, bar_y, fw.max(2.0), 16.0, AMBER, 3.0);
        svg += &svg_text(fx + fw + 8.0, y as f64 + row_h as f64 / 2.0,
            &format!("{:.1}%", floor_pct * 100.0), 8.5, INK2, "start");

        // Column 3: Silhouette Δ
        let dx  = pad_l as f64 + col_w * 3.0 + col_w / 2.0;
        let dcol = if delta_sil > 0.001 { GREEN }
                   else if delta_sil < -0.001 { RED } else { INK3 };
        let sign = if delta_sil >= 0.0 { "+" } else { "" };
        svg += &format!(
            r#"<circle cx="{:.1}" cy="{:.1}" r="5" fill="{dcol}"/>"#,
            dx - 20.0, y as f64 + row_h as f64 / 2.0
        );
        svg += &svg_text(dx - 10.0, y as f64 + row_h as f64 / 2.0,
            &format!("{sign}{:.4}", delta_sil), 9.0, dcol, "start");
    }

    // Column dividers
    for ci in 1..4 {
        let x = pad_l as f64 + col_w * ci as f64;
        svg += &svg_line(x, pad_t as f64 - 10.0, x,
                          (pad_t + n_clusters * row_h) as f64, INK3, 0.5);
    }

    svg += svg_close();
    svg
}

// ══════════════════════════════════════════════════════════════════════════════
// PLOT 8: Floor fraction per locus
// ══════════════════════════════════════════════════════════════════════════════
//
// Two sub-panels:
//   Top:    Top-25 loci by floor count (horizontal bars, sorted desc)
//   Bottom: Histogram of floor_count across all L_qc loci
//
// Per-locus floor count is computed here from x_raw + soup_contribution.
// No change to subtract_ambient() needed.
//
// Data: x_raw[L×N], soup_contribution[L], index_to_locus[L_qc→L_raw]
//       (index_to_locus maps compact metric index → original matrix row)
// ─────────────────────────────────────────────────────────────────────────────

pub fn plot_floor_per_locus(
    x_raw:             &ndarray::Array2<f64>,  // [L_raw × N]
    soup_contribution: &[f64],                  // [L_raw]
    index_to_locus:    &[usize],               // [L_qc] compact→original
    top_n:              usize,
) -> String {
    let l_qc  = index_to_locus.len();
    let n     = x_raw.ncols();
    let top_n = top_n.min(l_qc);

    // Compute per-locus floor count for QC-passing loci only
    let mut floor_counts: Vec<(usize, usize, usize)> = Vec::with_capacity(l_qc);
    // (compact_idx, original_locus, floor_count)

    for (ci, &orig_li) in index_to_locus.iter().enumerate() {
        if orig_li >= soup_contribution.len() { continue; }
        let s = soup_contribution[orig_li];
        if orig_li >= x_raw.nrows() { continue; }
        let mut fc = 0usize;
        for ni in 0..n {
            let v = x_raw[[orig_li, ni]];
            if !v.is_nan() && v < s { fc += 1; }
        }
        floor_counts.push((ci, orig_li, fc));
    }

    // Sort descending by floor count
    floor_counts.sort_by(|a, b| b.2.cmp(&a.2));

    let w = 900usize;
    let top_area_h = 420usize;
    let hist_y = top_area_h + 40;
    let hist_area_h = 160usize;
    let total_h = hist_y + hist_area_h + 40;
    let pad_l = 90.0_f64;
    let pad_r = 30.0_f64;
    let pad_t = 50.0_f64;
    let plot_w = w as f64 - pad_l - pad_r;
    let bar_plot_h = top_area_h as f64 - pad_t - 10.0;
    let bar_h_each = (bar_plot_h / top_n as f64).min(16.0).max(5.0);

    let max_fc = floor_counts.first().map(|t| t.2).unwrap_or(1).max(1);
    let pct_of_n = |fc: usize| fc as f64 / n.max(1) as f64 * 100.0;

    let mut svg = svg_open(w, total_h);

    // Title
    svg += &svg_text(w as f64 / 2.0, 26.0,
        &format!("Floor Fraction Per Locus — Top {} of {} QC loci  (N={} cells)",
                 top_n, l_qc, n),
        13.0, INK, "middle");

    svg += &svg_rect(pad_l, pad_t, plot_w, bar_plot_h + 10.0, SURFACE, 6.0);

    // Bars for top-N loci
    for (rank, &(ci, orig_li, fc)) in floor_counts.iter().take(top_n).enumerate() {
        let y    = pad_t + rank as f64 * bar_h_each + 3.0;
        let bw   = fc as f64 / max_fc as f64 * plot_w;
        let pct  = pct_of_n(fc);

        // Colour intensity by how severely floored
        let col = if pct > 60.0 { RED }
                  else if pct > 30.0 { AMBER }
                  else { VIOLET };

        svg += &svg_rect(pad_l, y, bw.max(2.0), bar_h_each - 1.0, col, 2.0);

        // Locus label
        svg += &format!(
            r#"<text x="{:.1}" y="{:.1}" font-size="7.5" fill="{INK3}" \
               text-anchor="end" dominant-baseline="middle">l={orig_li}(q{ci})</text>"#,
            pad_l - 4.0, y + bar_h_each / 2.0 - 0.5
        );

        // Value labels
        if bw > 80.0 {
            svg += &format!(
                r#"<text x="{:.1}" y="{:.1}" font-size="8" fill="{BG}" \
                   text-anchor="end" dominant-baseline="middle">{fc} cells  ({:.1}%)</text>"#,
                pad_l + bw - 4.0, y + bar_h_each / 2.0 - 0.5, pct
            );
        } else {
            svg += &format!(
                r#"<text x="{:.1}" y="{:.1}" font-size="8" fill="{INK3}" \
                   text-anchor="start" dominant-baseline="middle">{fc} ({:.1}%)</text>"#,
                pad_l + bw + 4.0, y + bar_h_each / 2.0 - 0.5, pct
            );
        }
    }

    // X axis
    for ti in 0..=4 {
        let fc_v = max_fc as f64 * ti as f64 / 4.0;
        let x    = pad_l + fc_v / max_fc as f64 * plot_w;
        svg += &svg_line(x, pad_t + bar_plot_h + 10.0, x, pad_t + bar_plot_h + 16.0, INK3, 1.0);
        svg += &svg_text(x, pad_t + bar_plot_h + 24.0,
            &format!("{:.0}", fc_v), 8.0, INK3, "middle");
        svg += &svg_text(x, pad_t + bar_plot_h + 34.0,
            &format!("({:.0}%)", pct_of_n(fc_v as usize)), 7.5, INK3, "middle");
    }
    svg += &svg_text(pad_l + plot_w / 2.0, pad_t + bar_plot_h + 44.0,
        "Cells floored at this locus  (count, % of N)", 9.0, INK3, "middle");

    // ── Histogram: log-scale distribution of floor counts ─────────────────
    // The distribution is extremely right-skewed: most loci have 0 floors,
    // a few loci have hundreds to thousands. Log scale is essential.
    let all_fcs: Vec<f64> = floor_counts.iter().map(|t| t.2 as f64).collect();
    if !all_fcs.is_empty() {
        let hist_y_f = hist_y as f64;
        let hist_h_f = hist_area_h as f64 - 24.0;
        svg += &svg_rect(pad_l, hist_y_f, plot_w, hist_h_f, SURFACE, 6.0);
        svg += &svg_text(w as f64 / 2.0, hist_y_f - 10.0,
            &format!("Floor count distribution — all {} QC loci  (log₁₀ y-axis)", l_qc),
            10.0, INK2, "middle");

        let n_bins   = 35usize;
        let fc_max_h = all_fcs.iter().cloned().fold(0.0_f64, f64::max);
        let bw_data  = fc_max_h / n_bins as f64;
        let mut counts = vec![0usize; n_bins];
        for &v in &all_fcs {
            let b = (v / (bw_data + 1e-14)) as usize;
            counts[b.min(n_bins - 1)] += 1;
        }
        // Log scale: log10(count+1) so zero-count bins vanish
        let log_max = counts.iter().map(|&c| ((c + 1) as f64).log10())
            .fold(0.0_f64, f64::max).max(0.1);
        let bw_px = plot_w / n_bins as f64;

        for (b, &c) in counts.iter().enumerate() {
            if c == 0 { continue; }
            let x      = pad_l + b as f64 * bw_px;
            let log_c  = ((c + 1) as f64).log10();
            let bh     = (log_c / log_max * hist_h_f).max(1.0);
            // Colour: first bin (mostly-zero) in cyan, rest in violet
            let col = if b == 0 { CYAN } else { VIOLET };
            svg += &svg_rect(x + 0.5, hist_y_f + hist_h_f - bh,
                              bw_px - 1.0, bh, col, 1.0);
        }

        // Log Y-axis ticks
        let max_pow = log_max.ceil() as usize;
        for p in 0..=max_pow {
            let y = hist_y_f + hist_h_f - (p as f64 / log_max * hist_h_f);
            svg += &svg_line(pad_l - 4.0, y, pad_l, y, INK3, 1.0);
            svg += &svg_text(pad_l - 6.0, y,
                &format!("10^{}", p), 7.0, INK3, "end");
        }
        svg += &format!(
            r#"<text x="{:.1}" y="{:.1}" font-size="7.5" fill="{INK3}" text-anchor="middle" transform="rotate(-90,{:.1},{:.1})">count (log)</text>"#,
            pad_l - 18.0, hist_y_f + hist_h_f / 2.0,
            pad_l - 18.0, hist_y_f + hist_h_f / 2.0);

        // X axis ticks
        for ti in 0..=5 {
            let v = fc_max_h * ti as f64 / 5.0;
            let x = pad_l + v / (fc_max_h + 1e-14) * plot_w;
            svg += &svg_text(x, hist_y_f + hist_h_f + 14.0,
                &format!("{:.0}", v), 8.0, INK3, "middle");
        }
        svg += &svg_text(pad_l + plot_w / 2.0, hist_y_f + hist_h_f + 28.0,
            "Floor count per locus", 9.0, INK3, "middle");

        // Annotation: median + mean floor count per locus
        let zero_floor = all_fcs.iter().filter(|&&v| v == 0.0).count();
        let nonzero    = l_qc - zero_floor;
        let mut sorted_fcs = all_fcs.clone();
        sorted_fcs.sort_by(|a, b| a.total_cmp(b));
        let median_fc = sorted_fcs[sorted_fcs.len() / 2] as usize;
        let mean_fc   = all_fcs.iter().sum::<f64>() / all_fcs.len().max(1) as f64;
        svg += &svg_text(pad_l + 6.0, hist_y_f + 14.0,
            &format!("{}/{} QC loci have ≥1 floored cell ({:.1}%)  ·  median = {} cells/locus  ·  mean = {:.0}  ·  only {} loci fully above soup",
                     nonzero, l_qc,
                     nonzero as f64 / l_qc.max(1) as f64 * 100.0,
                     median_fc, mean_fc, zero_floor),
            7.0, INK2, "start");
    }

    svg += svg_close();
    svg
}

// ══════════════════════════════════════════════════════════════════════════════
// PLOT 9: RMSE linearity proof
// ══════════════════════════════════════════════════════════════════════════════
//
// Scatter plot of (α_sim, RMSE) with:
//   • Dots for each simulation level
//   • OLS fitted line through origin (RMSE = slope × α)
//   • R² annotation
//   • Slope annotation (expected ≈ 0.367)
//   • Residual panel below (RMSE − slope×α)
//
// Data: BenchmarkResult.points[].alpha_sim, .rmse — already in BenchmarkPoint.
// Regression computed here.
// ─────────────────────────────────────────────────────────────────────────────

pub fn plot_rmse_linearity(result: &BenchmarkResult) -> String {
    if result.points.len() < 2 {
        let mut s = svg_open(800, 200);
        s += &svg_text(400.0, 100.0, "Need ≥2 simulation levels for linearity proof",
                       12.0, INK3, "middle");
        s += svg_close();
        return s;
    }

    let pts = &result.points;
    let xs: Vec<f64> = pts.iter().map(|p| p.alpha_sim).collect();
    let ys: Vec<f64> = pts.iter().map(|p| p.rmse).collect();
    let n_pts = xs.len();

    // OLS regression through origin: slope = Σ(xy) / Σ(x²)
    let sum_xy: f64 = xs.iter().zip(ys.iter()).map(|(x,y)| x*y).sum();
    let sum_x2: f64 = xs.iter().map(|x| x*x).sum();
    let slope   = if sum_x2 > 1e-15 { sum_xy / sum_x2 } else { 0.0 };

    // R² (vs mean RMSE)
    let mean_y = ys.iter().sum::<f64>() / n_pts as f64;
    let ss_tot: f64 = ys.iter().map(|y| (y - mean_y).powi(2)).sum();
    let ss_res: f64 = xs.iter().zip(ys.iter()).map(|(x,y)| (y - slope*x).powi(2)).sum();
    let r2 = if ss_tot > 1e-15 { 1.0 - ss_res / ss_tot } else { 1.0 };

    // Residuals
    let residuals: Vec<f64> = xs.iter().zip(ys.iter())
        .map(|(x, y)| y - slope * x).collect();
    let res_max = residuals.iter().cloned().map(f64::abs)
        .fold(0.0_f64, f64::max).max(1e-10);

    let w = 800usize;
    let main_h = 420usize;
    let res_h  = 160usize;
    let total_h = main_h + res_h + 60;
    let pad_l = 80.0_f64; let pad_r = 30.0_f64;
    let pad_t = 56.0_f64; let pad_b = 50.0_f64;
    let plot_w = w as f64 - pad_l - pad_r;
    let plot_h = main_h as f64 - pad_t - pad_b;

    let x_min = 0.0f64;
    let x_max = xs.iter().cloned().fold(0.0_f64, f64::max) * 1.1;
    let y_min = 0.0f64;
    let y_max = ys.iter().cloned().fold(0.0_f64, f64::max) * 1.15;

    let to_x = |v: f64| pad_l + (v - x_min) / (x_max - x_min + 1e-14) * plot_w;
    let to_y = |v: f64| pad_t + plot_h - (v - y_min) / (y_max - y_min + 1e-14) * plot_h;

    let mut svg = svg_open(w, total_h);

    // Title
    svg += &svg_text(w as f64 / 2.0, 28.0,
        "RMSE Linearity Proof — Recovery Quality vs Contamination Level",
        13.0, INK, "middle");

    svg += &svg_rect(pad_l, pad_t, plot_w, plot_h, SURFACE, 6.0);

    // Grid
    for ti in 0..=4 {
        let y = pad_t + plot_h * ti as f64 / 4.0;
        svg += &svg_line(pad_l, y, pad_l + plot_w, y, BORDER, 1.0);
        let v = y_max - (y_max - y_min) * ti as f64 / 4.0;
        svg += &svg_text(pad_l - 6.0, y, &format!("{:.4}", v), 8.5, INK3, "end");
    }
    for ti in 0..=5 {
        let x = pad_l + plot_w * ti as f64 / 5.0;
        svg += &svg_line(x, pad_t, x, pad_t + plot_h, BORDER, 1.0);
        let v = x_min + (x_max - x_min) * ti as f64 / 5.0;
        svg += &svg_text(x, pad_t + plot_h + 16.0, &format!("{:.3}", v), 8.5, INK3, "middle");
    }

    // Fitted line
    let lx0 = to_x(x_min); let ly0 = to_y(slope * x_min);
    let lx1 = to_x(x_max); let ly1 = to_y(slope * x_max);
    svg += &format!(
        r#"<line x1="{lx0:.1}" y1="{ly0:.1}" x2="{lx1:.1}" y2="{ly1:.1}" \
           stroke="{GREEN}" stroke-width="2" stroke-dasharray="8,4" opacity="0.8"/>"#
    );

    // True-α reference
    let tx = to_x(result.true_alpha);
    svg += &format!(
        r#"<line x1="{tx:.1}" y1="{pad_t:.1}" x2="{tx:.1}" y2="{:.1}" \
           stroke="{AMBER}" stroke-width="1.5" stroke-dasharray="4,3" opacity="0.8"/>"#,
        pad_t + plot_h
    );
    svg += &svg_text(tx + 4.0, pad_t + 14.0,
        &format!("true α={:.4}", result.true_alpha), 8.0, AMBER, "start");

    // Scatter dots
    for (i, p) in pts.iter().enumerate() {
        let cx = to_x(p.alpha_sim);
        let cy = to_y(p.rmse);
        svg += &format!(r#"<circle cx="{cx:.1}" cy="{cy:.1}" r="5.5" fill="{CYAN}"/>"#);
        svg += &svg_text(cx, cy - 12.0,
            &format!("({:.2},{:.4})", p.alpha_sim, p.rmse), 8.0, INK2, "middle");
        let _ = i;
    }

    // Annotation box
    svg += &svg_rect(pad_l + 8.0, pad_t + 8.0, 230.0, 70.0, BG, 5.0);
    svg += &svg_text(pad_l + 16.0, pad_t + 22.0,
        &format!("slope  = {:.4}  (RMSE / α)", slope), 9.0, GREEN, "start");
    svg += &svg_text(pad_l + 16.0, pad_t + 36.0,
        &format!("R²     = {:.6}", r2), 9.0, if r2 > 0.999 { GREEN } else { AMBER }, "start");
    svg += &svg_text(pad_l + 16.0, pad_t + 50.0,
        &format!("n pts  = {} levels × {} trials avg", n_pts, 1), 9.0, INK2, "start");
    svg += &svg_text(pad_l + 16.0, pad_t + 64.0,
        if r2 > 0.999 { "✓ Perfect linearity confirmed" } else { "⚠ Check simulation levels" },
        9.0, if r2 > 0.999 { GREEN } else { AMBER }, "start");

    // Axis labels
    svg += &svg_text(pad_l + plot_w / 2.0, total_h as f64 - res_h as f64 - 30.0,
        "α (simulated contamination fraction)", 9.0, INK3, "middle");
    svg += &format!(
        r#"<text x="14" y="{:.1}" font-size="9" fill="{INK3}" text-anchor="middle" \
           transform="rotate(-90,14,{:.1})">RMSE (recovery error)</text>"#,
        pad_t + plot_h / 2.0, pad_t + plot_h / 2.0
    );

    // ── Residual panel ────────────────────────────────────────────────────
    let ry = main_h as f64 + 20.0;
    let rh = res_h as f64 - 30.0;
    svg += &svg_rect(pad_l, ry, plot_w, rh, SURFACE, 6.0);
    svg += &svg_text(w as f64 / 2.0, ry - 8.0,
        "Residuals  (RMSE − slope × α_sim)", 9.0, INK2, "middle");
    svg += &svg_line(pad_l, ry + rh / 2.0, pad_l + plot_w, ry + rh / 2.0, INK3, 1.0);

    for (i, (&xi, &ri)) in xs.iter().zip(residuals.iter()).enumerate() {
        let rx = to_x(xi);
        // Clamp dot to panel bounds to prevent rendering outside the panel
        let ry_raw = ry + rh / 2.0 - (ri / res_max * rh / 2.0);
        let ry_dot = ry_raw.max(ry + 4.0).min(ry + rh - 4.0);
        let rcol = if ri.abs() < res_max * 0.1 { GREEN } else { AMBER };
        svg += &format!(r#"<circle cx="{rx:.1}" cy="{ry_dot:.1}" r="4" fill="{rcol}"/>"#);
        svg += &svg_line(rx, ry + rh / 2.0, rx, ry_dot, rcol, 1.0);
        // Show actual residual value (not clamped)
        let label = if ri.abs() < 1e-6 { "≈0".to_string() } else { format!("{:.5}", ri) };
        svg += &svg_text(rx, ry + rh + 14.0, &label, 7.5, INK3, "middle");
        let _ = i;
    }

    svg += svg_close();
    svg
}
