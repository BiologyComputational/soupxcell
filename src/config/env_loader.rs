// ============================================================================
// config/env_loader.rs — .soupxcell.env reader
// ============================================================================

use std::collections::HashMap;

pub struct EnvMap(HashMap<String, String>);

impl EnvMap {
    pub fn load(path: &str) -> Result<Self, String> {
        let content = std::fs::read_to_string(path)
            .map_err(|e| format!("Cannot read {}: {}", path, e))?;
        let mut map = HashMap::new();
        for (lineno, line) in content.lines().enumerate() {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') { continue; }
            let parts: Vec<&str> = trimmed.splitn(2, '=').collect();
            if parts.len() != 2 {
                eprintln!("  ⚠  .soupxcell.env line {}: not KEY=VALUE, skipped", lineno + 1);
                continue;
            }
            let key = parts[0].trim().to_uppercase();
            if !key.starts_with("SOUPX_") {
                eprintln!("  ⚠  .soupxcell.env line {}: key {:?} does not start with SOUPX_, skipped", lineno + 1, key);
                continue;
            }
            // Strip optional surrounding quotes
            let raw_val = parts[1].trim();
            let val = if (raw_val.starts_with('"') && raw_val.ends_with('"'))
                      || (raw_val.starts_with('\'') && raw_val.ends_with('\'')) {
                raw_val[1..raw_val.len()-1].to_string()
            } else {
                raw_val.to_string()
            };
            // Tilde expansion
            let val = if val.starts_with('~') {
                std::env::var("HOME").map(|h| val.replacen('~', &h, 1)).unwrap_or(val)
            } else { val };
            map.insert(key, val);
        }
        Ok(EnvMap(map))
    }

    pub fn empty() -> Self { EnvMap(HashMap::new()) }

    pub fn str(&self, key: &str) -> Option<String> {
        self.0.get(key).cloned()
    }
    pub fn float(&self, key: &str) -> Option<f64> {
        self.0.get(key)?.parse().ok()
    }
    pub fn int(&self, key: &str) -> Option<i64> {
        self.0.get(key)?.parse().ok()
    }
    pub fn bool(&self, key: &str) -> Option<bool> {
        self.0.get(key).map(|v| matches!(v.to_lowercase().as_str(), "true" | "1" | "yes"))
    }
    #[allow(dead_code)]
    pub fn len(&self) -> usize { self.0.len() }
}
