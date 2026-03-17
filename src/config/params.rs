// ============================================================================
// config/params.rs — CLI parsing, Params struct, 3-layer merge
// ============================================================================
// Priority chain: CLI flags > --config JSON > .soupxcell.env > built-in defaults
// ============================================================================


use crate::config::config_loader::ConfigProfile;
use crate::config::env_loader::EnvMap;

// ── Embedding method enum ─────────────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq)]
pub enum EmbedChoice {
    Tsne,
    Umap,
    Both,
    None,
}

impl EmbedChoice {
    pub fn from_str(s: &str) -> Self {
        match s.to_lowercase().as_str() {
            "tsne"  => EmbedChoice::Tsne,
            "umap"  => EmbedChoice::Umap,
            "both"  => EmbedChoice::Both,
            "none"  => EmbedChoice::None,
            _       => EmbedChoice::Tsne,
        }
    }
    pub fn run_tsne(&self) -> bool { matches!(self, EmbedChoice::Tsne | EmbedChoice::Both) }
    pub fn run_umap(&self) -> bool { matches!(self, EmbedChoice::Umap | EmbedChoice::Both) }
}

// ── Main Params struct ────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct Params {
    // ── Required inputs
    pub ref_matrix:     String,
    pub alt_matrix:     String,
    pub clusters:       String,
    pub ambient:        String,

    // ── Output
    pub output:         String,
    pub plot_dir:       String,

    // ── Correction
    pub alpha_override: Option<f64>,
    pub singlets_only:  bool,

    // ── Embedding
    pub embed:          EmbedChoice,
    pub pca_components: usize,
    pub tsne_perplexity: f64,
    pub tsne_iter:      usize,
    pub seed:           u64,

    // ── Simulation
    pub simulate:       bool,
    pub sim_levels:     Vec<f64>,
    pub sim_trials:     usize,

    // ── System
    pub threads:        usize,
    pub dry_run:        bool,
    pub dry_run_yes:    bool,

    // ── Provenance (for logging)
    pub config_path:    Option<String>,
    pub env_path:       Option<String>,
    pub param_sources:  Vec<(String, String)>,  // (param, source)
    // Locus QC filter (same 4 params as souporcell — applied during matrix loading)
    pub min_ref:            u32,
    pub min_alt:            u32,
    pub min_ref_umis:       u32,
    pub min_alt_umis:       u32,
    // Metrics tuning — locus selection (priority: vcf > err > n_qc_loci > metric_loci)
    pub souporcell_vcf:     Option<String>,  // cluster_genotypes.vcf — primary locus source
    pub freebayes_vcf:      Option<String>,  // freebayes VCF for exact CHROM:POS mapping
    pub souporcell_err:     Option<String>,  // clusters.err path (parses "Loci passing QC: N")
    pub n_qc_loci:          usize,           // explicit N (0 = not set)
    pub metric_loci:        usize,           // 0 = all covered loci
    pub metric_cells:       usize,           // cells per cluster for silhouette
    pub min_locus_variance: f64,             // 0 = use top-N by variance
    // GUI
    pub gui:            bool,
    pub gui_port:       u16,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            ref_matrix:     String::new(),
            alt_matrix:     String::new(),
            clusters:       String::new(),
            ambient:        String::new(),
            output:         "./soupxcell_output".into(),
            plot_dir:       "./soupxcell_output/figures".into(),
            alpha_override: None,
            singlets_only:  true,
            embed:          EmbedChoice::Tsne,
            pca_components: 50,
            tsne_perplexity: 30.0,
            tsne_iter:      1000,
            seed:           42,
            simulate:       false,
            sim_levels:     vec![0.01, 0.05, 0.10, 0.20, 0.30, 0.50],
            sim_trials:     3,
            threads:        1,
            dry_run:        false,
            dry_run_yes:    false,
            config_path:    None,
            env_path:       None,
            param_sources:  Vec::new(),
            min_ref:            4,
            min_alt:            4,
            min_ref_umis:       4,
            min_alt_umis:       4,
            souporcell_vcf:     None,
            freebayes_vcf:      None,
            souporcell_err:     None,
            n_qc_loci:          0,
            metric_loci:        0,
            metric_cells:       150,
            min_locus_variance: 0.0,
            gui:            false,
            gui_port:       7878,
        }
    }
}

// ── Tilde expansion ───────────────────────────────────────────────────────────

pub fn expand_tilde(path: &str) -> String {
    if path.starts_with('~') {
        if let Ok(home) = std::env::var("HOME") {
            return path.replacen('~', &home, 1);
        }
    }
    path.to_string()
}

// ── Resolve macros (same pattern as souporcell) ───────────────────────────────

macro_rules! resolve_str {
    ($cli:expr, $cfg:expr, $env:expr, $default:expr) => {
        if let Some(v) = $cli { expand_tilde(&v) }
        else if let Some(v) = $cfg { expand_tilde(&v) }
        else if let Some(v) = $env { expand_tilde(&v) }
        else { $default.to_string() }
    }
}

macro_rules! resolve_float {
    ($cli:expr, $cfg:expr, $env:expr, $default:expr) => {
        if let Some(v) = $cli { v }
        else if let Some(v) = $cfg { v }
        else if let Some(v) = $env { v }
        else { $default }
    }
}

macro_rules! resolve_int {
    ($cli:expr, $cfg:expr, $env:expr, $default:expr) => {
        if let Some(v) = $cli { v }
        else if let Some(v) = $cfg { v }
        else if let Some(v) = $env { v }
        else { $default }
    }
}

macro_rules! resolve_bool_flag {
    ($cli:expr, $env:expr, $default:expr) => {
        if $cli { true }
        else if let Some(v) = $env { v }
        else { $default }
    }
}

// ── Main load_params() ───────────────────────────────────────────────────────

pub fn load_params() -> Params {
    // ── Parse raw CLI ─────────────────────────────────────────────────────────
    let yaml   = load_yaml!("params.yml");
    let app    = clap::App::from_yaml(yaml);
    let matches = app.get_matches();

    // ── Load config profile (--config) ────────────────────────────────────────
    let config_path = matches.value_of("config").map(String::from);
    let cfg = config_path.as_deref()
        .map(|p| ConfigProfile::load(p).unwrap_or_else(|e| {
            eprintln!("  ⚠  Failed to load config {}: {}", p, e);
            ConfigProfile::empty()
        }))
        .unwrap_or_else(ConfigProfile::empty);

    // ── Load env file (.soupxcell.env) ────────────────────────────────────────
    let explicit_env = matches.value_of("env").map(String::from);
    let env_path = explicit_env.clone().or_else(|| {
        // Auto-discover: CWD → ~/.soupxcell.env
        let cwd_env = std::env::current_dir()
            .ok()
            .map(|d| d.join(".soupxcell.env"))
            .filter(|p| p.exists())
            .and_then(|p| p.to_str().map(String::from));
        if cwd_env.is_some() { return cwd_env; }
        dirs_next::home_dir()
            .map(|h| h.join(".soupxcell.env"))
            .filter(|p| p.exists())
            .and_then(|p| p.to_str().map(String::from))
    });

    let env = env_path.as_deref()
        .map(|p| EnvMap::load(p).unwrap_or_else(|e| {
            eprintln!("  ⚠  Failed to load env {}: {}", p, e);
            EnvMap::empty()
        }))
        .unwrap_or_else(EnvMap::empty);

    // ── Resolve all params ────────────────────────────────────────────────────
    let ref_matrix = resolve_str!(
        matches.value_of("ref_matrix").map(String::from),
        cfg.str("ref_matrix"), env.str("SOUPX_REF"), ""
    );
    let alt_matrix = resolve_str!(
        matches.value_of("alt_matrix").map(String::from),
        cfg.str("alt_matrix"), env.str("SOUPX_ALT"), ""
    );
    let clusters = resolve_str!(
        matches.value_of("clusters").map(String::from),
        cfg.str("clusters"), env.str("SOUPX_CLUSTERS"), ""
    );
    let ambient = resolve_str!(
        matches.value_of("ambient").map(String::from),
        cfg.str("ambient"), env.str("SOUPX_AMBIENT"), ""
    );
    let output = resolve_str!(
        matches.value_of("output").map(String::from),
        cfg.str("output"), env.str("SOUPX_OUTPUT"), "./soupxcell_output"
    );
    let plot_dir_raw = resolve_str!(
        matches.value_of("plot_dir").map(String::from),
        cfg.str("plot_dir"), env.str("SOUPX_PLOT_DIR"), ""
    );
    let plot_dir = if plot_dir_raw.is_empty() {
        format!("{}/figures", output.trim_end_matches('/'))
    } else { plot_dir_raw };

    let alpha_override_cli = matches.value_of("alpha_override")
        .and_then(|v| v.parse::<f64>().ok());
    let alpha_override = if let Some(v) = alpha_override_cli { Some(v) }
        else if let Some(v) = cfg.float("alpha") { Some(v) }
        else if let Some(v) = env.float("SOUPX_ALPHA") { Some(v) }
        else { None };

    let threads = resolve_int!(
        matches.value_of("threads").and_then(|v| v.parse::<usize>().ok()),
        cfg.int("threads").map(|v| v as usize),
        env.int("SOUPX_THREADS").map(|v| v as usize),
        1usize
    );

    let embed_str = resolve_str!(
        matches.value_of("embed").map(String::from),
        cfg.str("embed"), env.str("SOUPX_EMBED"), "tsne"
    );
    let embed = EmbedChoice::from_str(&embed_str);

    let pca_components = resolve_int!(
        matches.value_of("pca_components").and_then(|v| v.parse::<usize>().ok()),
        cfg.int("pca_components").map(|v| v as usize),
        env.int("SOUPX_PCA_COMPONENTS").map(|v| v as usize),
        50usize
    );

    let tsne_perplexity = resolve_float!(
        matches.value_of("tsne_perplexity").and_then(|v| v.parse::<f64>().ok()),
        cfg.float("tsne_perplexity"),
        env.float("SOUPX_TSNE_PERPLEXITY"),
        30.0_f64
    );

    let tsne_iter = resolve_int!(
        matches.value_of("tsne_iter").and_then(|v| v.parse::<usize>().ok()),
        cfg.int("tsne_iter").map(|v| v as usize),
        env.int("SOUPX_TSNE_ITER").map(|v| v as usize),
        1000usize
    );

    let seed = resolve_int!(
        matches.value_of("seed").and_then(|v| v.parse::<u64>().ok()),
        cfg.int("seed").map(|v| v as u64),
        env.int("SOUPX_SEED").map(|v| v as u64),
        42u64
    );

    let simulate = resolve_bool_flag!(
        matches.is_present("simulate"),
        env.bool("SOUPX_SIMULATE"),
        false
    );

    // sim_levels — multiple values
    let sim_levels_cli: Vec<f64> = matches.values_of("sim_levels")
        .map(|vals| vals.filter_map(|v| v.parse::<f64>().ok()).collect())
        .unwrap_or_default();
    let sim_levels = if !sim_levels_cli.is_empty() { sim_levels_cli }
        else if let Some(v) = env.str("SOUPX_SIM_LEVELS") {
            v.split(',').filter_map(|s| s.trim().parse::<f64>().ok()).collect()
        } else { vec![0.01, 0.05, 0.10, 0.20, 0.30, 0.50] };

    let sim_trials = resolve_int!(
        matches.value_of("sim_trials").and_then(|v| v.parse::<usize>().ok()),
        cfg.int("sim_trials").map(|v| v as usize),
        env.int("SOUPX_SIM_TRIALS").map(|v| v as usize),
        3usize
    );

    let singlets_only = resolve_bool_flag!(
        matches.is_present("singlets_only"),
        env.bool("SOUPX_SINGLETS_ONLY"),
        true
    );

    let dry_run = resolve_bool_flag!(
        matches.is_present("dry_run"),
        env.bool("SOUPX_DRY_RUN"),
        false
    );
    let dry_run_yes = resolve_bool_flag!(
        matches.is_present("dry_run_yes"),
        env.bool("SOUPX_DRY_RUN_YES"),
        false
    );

    let min_ref = resolve_int!(
        matches.value_of("min_ref").and_then(|v| v.parse::<u32>().ok()),
        cfg.int("min_ref").map(|v| v as u32),
        env.int("SOUPX_MIN_REF").map(|v| v as u32),
        4u32
    );
    let min_alt = resolve_int!(
        matches.value_of("min_alt").and_then(|v| v.parse::<u32>().ok()),
        cfg.int("min_alt").map(|v| v as u32),
        env.int("SOUPX_MIN_ALT").map(|v| v as u32),
        4u32
    );
    let min_ref_umis = resolve_int!(
        matches.value_of("min_ref_umis").and_then(|v| v.parse::<u32>().ok()),
        cfg.int("min_ref_umis").map(|v| v as u32),
        env.int("SOUPX_MIN_REF_UMIS").map(|v| v as u32),
        4u32
    );
    let min_alt_umis = resolve_int!(
        matches.value_of("min_alt_umis").and_then(|v| v.parse::<u32>().ok()),
        cfg.int("min_alt_umis").map(|v| v as u32),
        env.int("SOUPX_MIN_ALT_UMIS").map(|v| v as u32),
        4u32
    );
    let freebayes_vcf = matches.value_of("freebayes_vcf").map(String::from)
        .or_else(|| env.str("SOUPX_FREEBAYES_VCF"));
    let souporcell_vcf = matches.value_of("souporcell_vcf").map(String::from)
        .or_else(|| env.str("SOUPX_SOUPORCELL_VCF"));
    let souporcell_err = matches.value_of("souporcell_err").map(String::from)
        .or_else(|| env.str("SOUPX_SOUPORCELL_ERR"));
    let n_qc_loci = resolve_int!(
        matches.value_of("n_qc_loci").and_then(|v| v.parse::<usize>().ok()),
        cfg.int("n_qc_loci").map(|v| v as usize),
        env.int("SOUPX_N_QC_LOCI").map(|v| v as usize),
        0usize
    );
    let metric_loci = resolve_int!(
        matches.value_of("metric_loci").and_then(|v| v.parse::<usize>().ok()),
        cfg.int("metric_loci").map(|v| v as usize),
        env.int("SOUPX_METRIC_LOCI").map(|v| v as usize),
        0usize
    );
    let metric_cells = resolve_int!(
        matches.value_of("metric_cells").and_then(|v| v.parse::<usize>().ok()),
        cfg.int("metric_cells").map(|v| v as usize),
        env.int("SOUPX_METRIC_CELLS").map(|v| v as usize),
        150usize
    );
    let min_locus_variance = resolve_float!(
        matches.value_of("min_locus_variance").and_then(|v| v.parse::<f64>().ok()),
        cfg.float("min_locus_variance"),
        env.float("SOUPX_MIN_LOCUS_VARIANCE"),
        0.0_f64
    );
    let gui = matches.is_present("gui");
    let gui_port = matches.value_of("gui_port")
        .and_then(|v| v.parse::<u16>().ok())
        .unwrap_or(7878);

    // ── Build provenance list ─────────────────────────────────────────────────
    let sources = format!("CLI{}{}",
        config_path.as_deref().map(|p| format!(" > config({})", p)).unwrap_or_default(),
        env_path.as_deref().map(|p| format!(" > env({})", p)).unwrap_or_default(),
    );

    Params {
        ref_matrix, alt_matrix, clusters, ambient,
        output, plot_dir,
        alpha_override,
        singlets_only,
        embed,
        pca_components,
        tsne_perplexity,
        tsne_iter,
        seed,
        simulate,
        sim_levels,
        sim_trials,
        threads,
        dry_run,
        dry_run_yes,
        config_path,
        env_path: env_path.clone(),
        param_sources: vec![("sources".into(), sources)],
        min_ref,
        min_alt,
        min_ref_umis,
        min_alt_umis,
        souporcell_vcf,
        freebayes_vcf,
        souporcell_err,
        n_qc_loci,
        metric_loci,
        metric_cells,
        min_locus_variance,
        gui,
        gui_port,
    }
}
