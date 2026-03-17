// ============================================================================
// infra/logger.rs — Production-grade structured stderr logging
// ============================================================================
// Every line carries a [SOUPX] prefix tag for machine parsing.
// All slow operations announce themselves BEFORE they start so the user
// always knows what is running right now.
//
// Tags (grep-friendly):
//   [SOUPX] SECTION  — major pipeline stage
//   [SOUPX] STEP     — sub-operation within a module
//   [SOUPX] INFO     — key=value parameter / stat
//   [SOUPX] METRIC   — numeric result (silhouette, RMSE, ...)
//   [SOUPX] PROGRESS — progress bar checkpoint
//   [SOUPX] OK       — completion with timing
//   [SOUPX] WRITTEN  — file written to disk
//   [SOUPX] WARN     — non-fatal warning
//   [SOUPX] ERROR    — fatal error (binary exits after)
//   [SOUPX] TIME     — arbitrary elapsed checkpoint
// ============================================================================

use std::time::{Duration, Instant};

fn is_tty() -> bool {
    std::env::var("TERM").map(|t| t != "dumb").unwrap_or(false)
}

// Colour helpers — only active in TTY (plain text in pipes/files)
fn c(code: &'static str) -> &'static str { if is_tty() { code } else { "" } }
const CYAN:   &str = "\x1b[36m";
const CYAN2:  &str = "\x1b[96m";
const YELLOW: &str = "\x1b[33m";
#[allow(dead_code)]
const GREEN:  &str = "\x1b[32m";
const GREEN2: &str = "\x1b[92m";
const RED:    &str = "\x1b[31m";
const PINK:   &str = "\x1b[35m";
const VIOLET: &str = "\x1b[95m";
const BOLD:   &str = "\x1b[1m";
const DIM:    &str = "\x1b[2m";
const RESET:  &str = "\x1b[0m";

pub fn fmt_dur(d: Duration) -> String {
    let s  = d.as_secs();
    let ms = d.subsec_millis();
    if s >= 3600 { format!("{}h {:02}m {:02}s", s/3600, (s%3600)/60, s%60) }
    else if s >= 60 { format!("{}m {:02}.{}s",  s/60, s%60, ms/100) }
    else if s >= 1  { format!("{}.{:03}s",       s, ms) }
    else            { format!("{}ms",            d.as_millis()) }
}

pub struct Logger { pub start: Instant }

impl Logger {
    pub fn new() -> Self { Logger { start: Instant::now() } }

    fn elapsed(&self) -> String { fmt_dur(self.start.elapsed()) }

    // ── Major section banner ──────────────────────────────────────────────────
    pub fn section(&self, title: &str) {
        eprintln!();
        eprintln!("{}{}╔══════════════════════════════════════════════════════════════════════════╗{}",
            c(BOLD), c(CYAN), c(RESET));
        eprintln!("{}{}║  {:<72}║{}",
            c(BOLD), c(CYAN),
            format!("{}", title),
            c(RESET));
        eprintln!("{}{}╚══════════════════════════════════════════════════════════════════════════╝{}",
            c(BOLD), c(CYAN), c(RESET));
        eprintln!("{}  [SOUPX] SECTION  {}  elapsed={}{}",
            c(DIM), title, self.elapsed(), c(RESET));
    }

    // ── Sub-section divider ───────────────────────────────────────────────────
    pub fn subsection(&self, title: &str) {
        eprintln!("\n{}{}  ┄┄ {} {}[+{}]{}",
            c(BOLD), c(CYAN2), title, c(DIM), self.elapsed(), c(RESET));
    }

    // ── Step: announces a computation BEFORE it starts ───────────────────────
    pub fn step(&self, step_num: &str, description: &str) {
        eprintln!();
        eprintln!("{}{}  ▶ STEP {}  {}{}",
            c(BOLD), c(VIOLET), step_num, description, c(RESET));
        eprintln!("{}  [SOUPX] STEP     {}  {}{}",
            c(DIM), step_num, description, c(RESET));
    }

    // ── Computing: "I am about to do X right now" ─────────────────────────────
    pub fn computing(&self, what: &str) {
        eprintln!("{}{}  ⧖  {}… {}[elapsed={}]{}",
            c(DIM), c(CYAN2), what, c(DIM), self.elapsed(), c(RESET));
        eprintln!("{}  [SOUPX] COMPUTING {}  elapsed={}{}",
            c(DIM), what, self.elapsed(), c(RESET));
    }

    // ── Info key-value ────────────────────────────────────────────────────────
    pub fn info(&self, key: &str, val: &str) {
        let dots = ".".repeat(38usize.saturating_sub(key.len()));
        eprintln!("{}  │ {}{}{}  {}{}{}",
            c(DIM), key, dots, c(DIM),
            c(YELLOW), val, c(RESET));
        eprintln!("{}  [SOUPX] INFO      {}={}{}",
            c(DIM), key, val, c(RESET));
    }

    // ── Numeric metric (result, not just a param) ─────────────────────────────
    pub fn metric(&self, key: &str, val: &str, unit: &str) {
        eprintln!("{}{}  ◆ METRIC  {:<32} = {}{}{} {}{}",
            c(BOLD), c(GREEN2), key,
            c(GREEN2), val, c(DIM), unit, c(RESET));
        eprintln!("{}  [SOUPX] METRIC    {}={}  unit={}{}",
            c(DIM), key, val, unit, c(RESET));
    }

    // ── OK / completion ───────────────────────────────────────────────────────
    pub fn ok(&self, msg: &str, since: Option<Instant>) {
        let timing = since
            .map(|t| fmt_dur(t.elapsed()))
            .unwrap_or_default();
        if timing.is_empty() {
            eprintln!("{}{}  ✔  {}{}",  c(GREEN2), c(BOLD), msg, c(RESET));
            eprintln!("{}  [SOUPX] OK       {}{}",  c(DIM), msg, c(RESET));
        } else {
            eprintln!("{}{}  ✔  {}{}  [took {}]{}",
                c(GREEN2), c(BOLD), msg, c(DIM), timing, c(RESET));
            eprintln!("{}  [SOUPX] OK       {}  took={}{}",
                c(DIM), msg, timing, c(RESET));
        }
    }

    // ── File written ──────────────────────────────────────────────────────────
    pub fn written(&self, filename: &str, path: &str) {
        let size_str = std::fs::metadata(path)
            .map(|m| {
                let b = m.len();
                if b > 1024*1024 { format!("{:.1} MB", b as f64 / 1_048_576.0) }
                else if b > 1024 { format!("{:.1} KB", b as f64 / 1024.0) }
                else { format!("{} B", b) }
            })
            .unwrap_or_else(|_| "?".to_string());
        eprintln!("{}{}  ✦  {:<35} {}{}  ({}){}",
            c(PINK), c(BOLD), filename, c(DIM), path, size_str, c(RESET));
        eprintln!("{}  [SOUPX] WRITTEN   {}  path={}  size={}{}",
            c(DIM), filename, path, size_str, c(RESET));
    }

    // ── Warning ───────────────────────────────────────────────────────────────
    pub fn warn(&self, msg: &str) {
        eprintln!("{}{}  ⚠  WARNING: {}{}",
            c(BOLD), c(YELLOW), msg, c(RESET));
        eprintln!("{}  [SOUPX] WARN      {}{}",
            c(DIM), msg, c(RESET));
    }

    // ── Fatal error ───────────────────────────────────────────────────────────
    pub fn fatal(&self, msg: &str) -> ! {
        eprintln!("{}{}  ✖  FATAL: {}{}",
            c(BOLD), c(RED), msg, c(RESET));
        eprintln!("{}  [SOUPX] ERROR     {}{}",
            c(DIM), msg, c(RESET));
        std::process::exit(1);
    }

    // ── Progress bar ──────────────────────────────────────────────────────────
    pub fn progress(&self, label: &str, current: usize, total: usize) {
        let pct    = if total > 0 { current * 100 / total } else { 0 };
        let filled = pct * 44 / 100;
        let bar    = "█".repeat(filled) + &"░".repeat(44 - filled);
        eprint!("\r{}  [{}{}{}] {:>3}%  {:<30}{}",
            c(DIM), c(CYAN), bar, c(DIM), pct, label, c(RESET));
        if current >= total {
            eprintln!("  {}✔{}", c(GREEN2), c(RESET));
            eprintln!("{}  [SOUPX] PROGRESS {}  done={}{}",
                c(DIM), label, total, c(RESET));
        }
    }

    // ── Progress bar with rate + ETA ──────────────────────────────────────────
    #[allow(dead_code)]
    pub fn progress_rate(&self, label: &str, current: usize,
                          total: usize, since: Instant) {
        let elapsed_s = since.elapsed().as_secs_f64().max(0.001);
        let rate      = current as f64 / elapsed_s;
        let eta       = if rate > 0.0 && current < total {
            format!("ETA {:.0}s", (total - current) as f64 / rate)
        } else if current >= total { "done".to_string() }
        else { "".to_string() };
        let pct    = if total > 0 { current * 100 / total } else { 0 };
        let filled = pct * 30 / 100;
        let bar    = "█".repeat(filled) + &"░".repeat(30 - filled);
        eprint!("\r{}  [{}{}{}] {:>3}%  {:<22} {:>8.0}/s  {:<10}{}",
            c(DIM), c(CYAN), bar, c(DIM), pct, label, rate, eta, c(RESET));
        if current >= total { eprintln!("  {}✔{}", c(GREEN2), c(RESET)); }
    }

    // ── Timing checkpoint (no section, just a timestamp line) ─────────────────
    pub fn checkpoint(&self, label: &str) {
        eprintln!("{}  [SOUPX] TIME      {}  elapsed={}{}",
            c(DIM), label, self.elapsed(), c(RESET));
    }

    // ── Blank separator ───────────────────────────────────────────────────────
    #[allow(dead_code)]
    pub fn blank(&self) { eprintln!(); }
}

// ── Banner ────────────────────────────────────────────────────────────────────

pub fn print_banner() {
    eprintln!("{}{}", c(CYAN), r#"
  ███████╗ ██████╗ ██╗   ██╗██████╗ ██╗  ██╗ ██████╗███████╗██╗     ██╗
  ██╔════╝██╔═══██╗██║   ██║██╔══██╗╚██╗██╔╝██╔════╝██╔════╝██║     ██║
  ███████╗██║   ██║██║   ██║██████╔╝ ╚███╔╝ ██║     █████╗  ██║     ██║
  ╚════██║██║   ██║██║   ██║██╔═══╝  ██╔██╗ ██║     ██╔══╝  ██║     ██║
  ███████║╚██████╔╝╚██████╔╝██║     ██╔╝ ██╗╚██████╗███████╗███████╗███████╗
  ╚══════╝ ╚═════╝  ╚═════╝ ╚═╝     ╚═╝  ╚═╝ ╚═════╝╚══════╝╚══════╝╚══════╝"#);
    eprintln!("{}", c(RESET));
    eprintln!("  {}{}Ambient RNA Correction for Genotype Demultiplexing{}",
        c(BOLD), c(CYAN2), c(RESET));
    eprintln!("  {}soupXcell v1.0  ·  Jahidul Arafat{}",
        c(DIM), c(RESET));
    eprintln!("  {}Post-processor for souporcell  ·  Heaton et al., Nature Methods 2020{}",
        c(DIM), c(RESET));
    let exe = std::env::current_exe()
        .map(|p| p.to_string_lossy().to_string())
        .unwrap_or_else(|_| "soupxcell".to_string());
    eprintln!("  {}Binary : {}{}",  c(DIM), exe, c(RESET));
    let ts = {
        use std::time::{SystemTime, UNIX_EPOCH};
        let s = SystemTime::now().duration_since(UNIX_EPOCH)
            .map(|d| d.as_secs()).unwrap_or(0);
        format!("UTC {:02}:{:02}:{:02}", s%86400/3600, s%3600/60, s%60)
    };
    eprintln!("  {}Started: {}{}",  c(DIM), ts, c(RESET));
    eprintln!();
}
