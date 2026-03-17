// ============================================================================
// report/utils.rs — Shared colour constants + helper functions
// ============================================================================

// ── Colour palette (no # inside format!() strings) ───────────────────────────
pub const C_TEAL:   &str = "#006A6A";
pub const C_VIOLET: &str = "#5B21B6";
pub const C_AMBER:  &str = "#C75B00";
#[allow(dead_code)]
pub const C_RED:    &str = "#B91C1C";
pub const C_GREEN:  &str = "#15803D";
pub const C_LIGHT:  &str = "#E2E8F0";
#[allow(dead_code)]
pub const C_INK:    &str = "#0F172A";
pub const C_INK2:   &str = "#334155";
pub const C_MID:    &str = "#64748B";
pub const C_BG:     &str = "#F8F7F4";

// ── Timestamp ─────────────────────────────────────────────────────────────────
pub fn chrono_now() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let secs   = SystemTime::now().duration_since(UNIX_EPOCH)
                     .map(|d| d.as_secs()).unwrap_or(0);
    let days   = secs / 86400;
    let year   = 1970 + days / 365;
    let month  = (days % 365) / 30 + 1;
    let day    = (days % 365) % 30 + 1;
    let hour   = (secs % 86400) / 3600;
    let minute = (secs % 3600) / 60;
    format!("{:04}-{:02}-{:02} {:02}:{:02} UTC", year, month, day, hour, minute)
}

// ── Format duration ───────────────────────────────────────────────────────────
pub fn fmt_dur(d: std::time::Duration) -> String {
    let s = d.as_secs_f32();
    if s < 60.0 { format!("{:.1}s", s) }
    else         { format!("{:.0}m {:.0}s", s / 60.0, s % 60.0) }
}

// ── Pill HTML helper ──────────────────────────────────────────────────────────
#[allow(dead_code)]
pub fn pill(label: &str, cls: &str) -> String {
    format!(r#"<span class="pill pill-{cls}">{label}</span>"#)
}
