// ============================================================================
// infra/gui_server.rs — Web GUI server for soupxcell --gui mode
// ============================================================================
//
// Launched when the user runs:  soupxcell --gui  [--gui_port 7878]
//
// Uses ONLY std::net — zero new Cargo.toml dependencies.
//
// Endpoints:
//   GET  /                   → serves the embedded soupxcell_gui.html
//   GET  /api/version        → {"version":"1.0.0","binary":"..."}
//   GET  /api/status         → {"status":"ready"|"running"}
//   POST /api/run            → spawns soupxcell with POSTed JSON config,
//                              streams stderr back as Server-Sent Events
//   GET  /api/files?dir=...  → lists files in the given directory
//   GET  /files/...          → serves a file from the given absolute path
//
// The GUI HTML is embedded at compile time — fully self-contained binary.
// ============================================================================

use std::io::{BufRead, BufReader, Read, Write};
use std::net::{TcpListener, TcpStream};
use std::process::{Command, Stdio};
use std::sync::{Arc, Mutex};
use std::thread;

// ── Embedded GUI HTML ─────────────────────────────────────────────────────────
static GUI_HTML: &str = include_str!("soupxcell_gui.html");

pub const VERSION:  &str = env!("CARGO_PKG_VERSION");
pub const PKG_NAME: &str = env!("CARGO_PKG_NAME");

// ── Active run state ──────────────────────────────────────────────────────────
struct RunState {
    running: bool,
    log:     Vec<String>,
}
impl RunState {
    fn new() -> Self { RunState { running: false, log: Vec::new() } }
}

// ── Entry point ───────────────────────────────────────────────────────────────

pub fn launch(port: u16) {
    let addr = format!("127.0.0.1:{}", port);
    let listener = TcpListener::bind(&addr).unwrap_or_else(|e| {
        eprintln!("GUI: cannot bind {} — {}", addr, e);
        std::process::exit(1);
    });

    let binary_path = std::env::current_exe()
        .map(|p| p.to_string_lossy().to_string())
        .unwrap_or_else(|_| "soupxcell".to_string());

    eprintln!("\n╔══════════════════════════════════════════════════════════╗");
    eprintln!("║   soupxcell v{}  —  GUI Mode                          ║", VERSION);
    eprintln!("╠══════════════════════════════════════════════════════════╣");
    eprintln!("║                                                          ║");
    eprintln!("║   Open in your browser:                                  ║");
    eprintln!("║   \x1b[36mhttp://localhost:{}\x1b[0m                              ║", port);
    eprintln!("║                                                          ║");
    eprintln!("║   Press  Ctrl+C  to stop the server.                     ║");
    eprintln!("╚══════════════════════════════════════════════════════════╝\n");

    let _ = open_browser(&format!("http://localhost:{}", port));

    let run_state = Arc::new(Mutex::new(RunState::new()));
    let binary    = Arc::new(binary_path);

    for stream in listener.incoming() {
        if let Ok(s) = stream {
            let rs  = Arc::clone(&run_state);
            let bin = Arc::clone(&binary);
            thread::spawn(move || handle_connection(s, rs, bin));
        }
    }
}

// ── Connection handler ────────────────────────────────────────────────────────

fn handle_connection(
    mut stream: TcpStream,
    run_state:  Arc<Mutex<RunState>>,
    binary:     Arc<String>,
) {
    let mut buf = [0u8; 32768];
    let n = match stream.read(&mut buf) {
        Ok(n) if n > 0 => n,
        _ => return,
    };
    let raw = String::from_utf8_lossy(&buf[..n]);
    let (method, path, body) = parse_request(&raw);

    match (method.as_str(), path.as_str()) {

        // ── Main GUI ──────────────────────────────────────────────────────────
        ("GET", "/") | ("GET", "/index.html") => {
            respond_html(&mut stream, 200, &inject_version(GUI_HTML));
        }

        // ── Version ───────────────────────────────────────────────────────────
        ("GET", "/api/version") => {
            respond_json(&mut stream, 200, &format!(
                r#"{{"version":"{}","binary":"{}","pkg":"{}"}}"#,
                VERSION, json_escape(&*binary), PKG_NAME
            ));
        }

        // ── Status ────────────────────────────────────────────────────────────
        ("GET", "/api/status") => {
            let running = run_state.lock().map(|s| s.running).unwrap_or(false);
            respond_json(&mut stream, 200, &format!(
                r#"{{"status":"{}"}}"#, if running { "running" } else { "ready" }
            ));
        }

        // ── Run soupxcell ─────────────────────────────────────────────────────
        ("POST", "/api/run") => {
            {
                let mut state = run_state.lock().unwrap();
                if state.running {
                    respond_json(&mut stream, 409,
                        r#"{"error":"A run is already in progress"}"#);
                    return;
                }
                state.running = true;
                state.log.clear();
            }

            // SSE headers
            let _ = stream.write_all(
                b"HTTP/1.1 200 OK\r\n\
                  Content-Type: text/event-stream\r\n\
                  Cache-Control: no-cache\r\n\
                  Access-Control-Allow-Origin: *\r\n\
                  Connection: keep-alive\r\n\r\n"
            );

            let args      = parse_config_to_args(&body);
            let bin_path  = (*binary).clone();
            let rs        = Arc::clone(&run_state);

            // Resolve output path so we redirect stdout to a file (not a pipe)
            // — prevents OS pipe-buffer deadlock on large datasets
            let output_dir = str_val(&body, "output")
                .unwrap_or_else(|| "./soupxcell_output".to_string());
            let _ = std::fs::create_dir_all(&output_dir);
            let stdout_path = format!("{}/soupxcell_stdout.log", output_dir.trim_end_matches('/'));
            let err_path    = format!("{}/soupxcell.err",        output_dir.trim_end_matches('/'));

            let stdout_redirect = std::fs::File::create(&stdout_path)
                .map(Stdio::from).unwrap_or_else(|_| Stdio::null());
            let err_file = std::fs::File::create(&err_path);

            let child = Command::new(&bin_path)
                .args(&args)
                .stdout(stdout_redirect)
                .stderr(Stdio::piped())
                .spawn();

            match child {
                Err(e) => {
                    let _ = stream.write_all(
                        format!("data: {{\"type\":\"error\",\"msg\":\"{}\"}}\n\n",
                            json_escape(&e.to_string())).as_bytes()
                    );
                    run_state.lock().map(|mut s| s.running = false).ok();
                }
                Ok(mut child) => {
                    // Start event
                    let _ = stream.write_all(format!(
                        "data: {{\"type\":\"start\",\"cmd\":\"{}\",\"output_dir\":\"{}\"}}\n\n",
                        json_escape(&format!("{} {}", bin_path, args.join(" "))),
                        json_escape(&output_dir)
                    ).as_bytes());

                    let mut err_writer = err_file.ok().map(std::io::BufWriter::new);

                    if let Some(stderr) = child.stderr.take() {
                        let reader = BufReader::new(stderr);
                        for line in reader.lines().flatten() {
                            // Classify the line type for the GUI to style it
                            let ev_type = classify_line(&line);
                            let event = format!(
                                "data: {{\"type\":\"{}\",\"msg\":\"{}\"}}\n\n",
                                ev_type, json_escape(&line)
                            );
                            println!("[GUI-RUN] {}", &line);
                            if let Some(ref mut w) = err_writer {
                                use std::io::Write as _;
                                let _ = w.write_all(line.as_bytes());
                                let _ = w.write_all(b"\n");
                            }
                            if stream.write_all(event.as_bytes()).is_err() { break; }
                            rs.lock().map(|mut s| s.log.push(line)).ok();
                        }
                    }

                    let exit_code = child.wait()
                        .map(|s| s.code().unwrap_or(-1)).unwrap_or(-1);
                    drop(err_writer);

                    let _ = stream.write_all(format!(
                        "data: {{\"type\":\"done\",\"exit_code\":{},\"output_dir\":\"{}\"}}\n\n",
                        exit_code, json_escape(&output_dir)
                    ).as_bytes());

                    println!("[GUI-RUN] [DONE] exit={}", exit_code);
                    run_state.lock().map(|mut s| s.running = false).ok();
                }
            }
        }

        // ── List files ────────────────────────────────────────────────────────
        ("GET", p) if p.starts_with("/api/files") => {
            let dir = extract_query_param(p, "dir")
                .unwrap_or_else(|| ".".to_string());
            respond_json(&mut stream, 200, &list_files(&dir));
        }

        // ── Serve a file ──────────────────────────────────────────────────────
        ("GET", p) if p.starts_with("/files/") => {
            serve_file(&mut stream, &percent_decode(&p[7..]));
        }

        // ── CORS preflight ────────────────────────────────────────────────────
        ("OPTIONS", _) => {
            let _ = stream.write_all(
                b"HTTP/1.1 204 No Content\r\n\
                  Access-Control-Allow-Origin: *\r\n\
                  Access-Control-Allow-Methods: GET, POST, OPTIONS\r\n\
                  Access-Control-Allow-Headers: Content-Type\r\n\r\n"
            );
        }

        _ => { respond_json(&mut stream, 404, r#"{"error":"not found"}"#); }
    }
}

// ── Line classifier (for SSE event type) ─────────────────────────────────────

fn classify_line(line: &str) -> &'static str {
    if line.starts_with("══")           { "section"  }
    else if line.starts_with("  ✔")     { "ok"       }
    else if line.starts_with("  │")     { "info"     }
    else if line.starts_with("  ⚠")     { "warn"     }
    else if line.starts_with("  ✖")     { "error"    }
    else if line.contains("Module 1")   { "module1"  }
    else if line.contains("Module 2")   { "module2"  }
    else if line.contains("Module 3")   { "module3"  }
    else if line.contains("COMPLETE")   { "complete" }
    else if line.contains("written")    { "written"  }
    else if line.contains("step")       { "step"     }
    else                                { "log"      }
}

// ── JSON body → CLI args ──────────────────────────────────────────────────────

fn str_val(body: &str, key: &str) -> Option<String> {
    let needle = format!("\"{}\":", key);
    let start  = body.find(&needle)?;
    let after  = body[start + needle.len()..].trim_start();
    if after.starts_with('"') {
        let inner = &after[1..];
        let end   = inner.find('"')?;
        let val   = inner[..end].trim().to_string();
        if val.is_empty() { None } else { Some(val) }
    } else {
        let end = after.find(|c: char| c == ',' || c == '}' || c == '\n')
                       .unwrap_or(after.len());
        let val = after[..end].trim().to_string();
        if val == "null" || val.is_empty() { None } else { Some(val) }
    }
}

fn bool_val(body: &str, key: &str) -> bool {
    str_val(body, key).map(|v| v == "true").unwrap_or(false)
}

fn parse_config_to_args(body: &str) -> Vec<String> {
    let mut args: Vec<String> = Vec::new();
    macro_rules! flag {
        ($flag:expr, $key:expr) => {
            if let Some(v) = str_val(body, $key) {
                args.push(format!("--{}", $flag));
                args.push(v);
            }
        };
    }
    flag!("ref",             "ref_matrix");
    flag!("alt",             "alt_matrix");
    flag!("clusters",        "clusters");
    flag!("ambient",         "ambient");
    flag!("output",          "output");
    flag!("alpha",           "alpha_override");
    flag!("threads",         "threads");
    flag!("embed",           "embed");
    flag!("pca_components",  "pca_components");
    flag!("tsne_perplexity", "tsne_perplexity");
    flag!("tsne_iter",       "tsne_iter");
    flag!("seed",            "seed");
    flag!("sim_levels",      "sim_levels");
    flag!("sim_trials",      "sim_trials");
    if bool_val(body, "singlets_only") { args.push("--singlets_only".into()); }
    if bool_val(body, "simulate")      { args.push("--simulate".into()); }
    args
}

// ── Version injection ─────────────────────────────────────────────────────────

fn inject_version(html: &str) -> String {
    html.replace("{{VERSION}}", VERSION)
        .replace("</head>", &format!(
            "<script>window.SOUPX_VERSION=\"{}\";window.SOUPX_SERVER=true;</script></head>",
            VERSION
        ))
}

// ── File listing ──────────────────────────────────────────────────────────────

fn list_files(dir: &str) -> String {
    let path = std::path::Path::new(dir);
    if !path.exists() || !path.is_dir() {
        return format!(r#"{{"error":"directory not found: {}","files":[]}}"#,
            json_escape(dir));
    }
    let mut entries: Vec<String> = Vec::new();
    if let Ok(rd) = std::fs::read_dir(path) {
        for entry in rd.flatten() {
            if let Ok(meta) = entry.metadata() {
                if meta.is_file() {
                    let name  = entry.file_name().to_string_lossy().to_string();
                    let fpath = entry.path().to_string_lossy().to_string();
                    entries.push(format!(
                        r#"{{"name":"{}","size":{},"path":"{}"}}"#,
                        json_escape(&name), meta.len(), json_escape(&fpath)
                    ));
                }
            }
        }
    }
    // Also check figures/ subdirectory
    let figures_dir = format!("{}/figures", dir.trim_end_matches('/'));
    if let Ok(rd) = std::fs::read_dir(&figures_dir) {
        for entry in rd.flatten() {
            if let Ok(meta) = entry.metadata() {
                if meta.is_file() {
                    let name  = format!("figures/{}", entry.file_name().to_string_lossy());
                    let fpath = entry.path().to_string_lossy().to_string();
                    entries.push(format!(
                        r#"{{"name":"{}","size":{},"path":"{}"}}"#,
                        json_escape(&name), meta.len(), json_escape(&fpath)
                    ));
                }
            }
        }
    }
    entries.sort();
    format!(r#"{{"dir":"{}","files":[{}]}}"#, json_escape(dir), entries.join(","))
}

// ── File serving ──────────────────────────────────────────────────────────────

fn serve_file(stream: &mut TcpStream, path: &str) {
    match std::fs::read(path) {
        Ok(bytes) => {
            let mime = if path.ends_with(".svg")  { "image/svg+xml" }
                  else if path.ends_with(".html") { "text/html" }
                  else if path.ends_with(".csv")  { "text/csv" }
                  else if path.ends_with(".txt")  { "text/plain" }
                  else if path.ends_with(".mtx")  { "text/plain" }
                  else { "application/octet-stream" };
            let fname = std::path::Path::new(path)
                .file_name().map(|n| n.to_string_lossy().to_string())
                .unwrap_or_else(|| "file".to_string());
            let header = format!(
                "HTTP/1.1 200 OK\r\nContent-Type: {}\r\nContent-Length: {}\r\n\
                 Content-Disposition: attachment; filename=\"{}\"\r\n\
                 Access-Control-Allow-Origin: *\r\n\r\n",
                mime, bytes.len(), fname
            );
            let _ = stream.write_all(header.as_bytes());
            let _ = stream.write_all(&bytes);
        }
        Err(_) => { respond_json(stream, 404, r#"{"error":"file not found"}"#); }
    }
}

// ── HTTP helpers ──────────────────────────────────────────────────────────────

fn respond_html(stream: &mut TcpStream, code: u16, body: &str) {
    let _ = stream.write_all(format!(
        "HTTP/1.1 {} {}\r\nContent-Type: text/html; charset=utf-8\r\n\
         Content-Length: {}\r\nAccess-Control-Allow-Origin: *\r\n\r\n{}",
        code, reason(code), body.len(), body
    ).as_bytes());
}

fn respond_json(stream: &mut TcpStream, code: u16, body: &str) {
    let _ = stream.write_all(format!(
        "HTTP/1.1 {} {}\r\nContent-Type: application/json\r\n\
         Content-Length: {}\r\nAccess-Control-Allow-Origin: *\r\n\r\n{}",
        code, reason(code), body.len(), body
    ).as_bytes());
}

fn reason(code: u16) -> &'static str {
    match code { 200 => "OK", 204 => "No Content", 404 => "Not Found",
                 409 => "Conflict", _ => "Unknown" }
}

fn parse_request(raw: &str) -> (String, String, String) {
    let mut lines = raw.lines();
    let first  = lines.next().unwrap_or("");
    let parts: Vec<&str> = first.splitn(3, ' ').collect();
    let method = parts.get(0).unwrap_or(&"GET").to_string();
    let path   = parts.get(1).unwrap_or(&"/")
        .splitn(2, '?').next().unwrap_or("/").to_string();
    let body = if let Some(i) = raw.find("\r\n\r\n") { raw[i+4..].to_string() }
               else if let Some(i) = raw.find("\n\n")     { raw[i+2..].to_string() }
               else { String::new() };
    (method, path, body)
}

fn extract_query_param(path: &str, key: &str) -> Option<String> {
    let q = path.find('?').map(|i| &path[i+1..]).unwrap_or("");
    for pair in q.split('&') {
        if let Some((k, v)) = pair.split_once('=') {
            if k == key { return Some(percent_decode(v)); }
        }
    }
    None
}

fn percent_decode(s: &str) -> String {
    let mut out = String::new();
    let mut chars = s.chars().peekable();
    while let Some(c) = chars.next() {
        if c == '%' {
            let h1 = chars.next().unwrap_or('0');
            let h2 = chars.next().unwrap_or('0');
            if let Ok(b) = u8::from_str_radix(&format!("{}{}", h1, h2), 16) {
                out.push(b as char);
            }
        } else if c == '+' { out.push(' ');
        } else { out.push(c); }
    }
    out
}

fn json_escape(s: &str) -> String {
    s.replace('\\', "\\\\").replace('"', "\\\"")
     .replace('\n', "\\n").replace('\r', "\\r").replace('\t', "\\t")
}

fn open_browser(url: &str) -> Result<(), ()> {
    #[cfg(target_os = "macos")]
    { Command::new("open").arg(url).spawn().map(|_| ()).map_err(|_| ()) }
    #[cfg(target_os = "linux")]
    { Command::new("xdg-open").arg(url).spawn().map(|_| ()).map_err(|_| ()) }
    #[cfg(target_os = "windows")]
    { Command::new("cmd").args(["/C","start",url]).spawn().map(|_| ()).map_err(|_| ()) }
    #[cfg(not(any(target_os="macos",target_os="linux",target_os="windows")))]
    { Err(()) }
}
