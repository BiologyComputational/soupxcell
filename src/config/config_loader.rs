// ============================================================================
// config/config_loader.rs — JSON parameter profile reader
// ============================================================================

use serde_json::Value;

pub struct ConfigProfile(Value);

impl ConfigProfile {
    pub fn load(path: &str) -> Result<Self, String> {
        let content = std::fs::read_to_string(path)
            .map_err(|e| format!("Cannot read config {}: {}", path, e))?;
        let v: Value = serde_json::from_str(&content)
            .map_err(|e| format!("Invalid JSON in {}: {}", path, e))?;
        // Accept both {params:{...}} and flat {...}
        let params = if let Some(p) = v.get("params") { p.clone() } else { v };
        Ok(ConfigProfile(params))
    }

    pub fn empty() -> Self { ConfigProfile(Value::Object(serde_json::Map::new())) }

    pub fn str(&self, key: &str) -> Option<String> {
        self.0.get(key)?.as_str().map(|s| {
            // Tilde expansion
            if s.starts_with('~') {
                std::env::var("HOME").map(|h| s.replacen('~', &h, 1)).unwrap_or_else(|_| s.to_string())
            } else { s.to_string() }
        })
    }
    pub fn float(&self, key: &str) -> Option<f64> {
        self.0.get(key)?.as_f64()
    }
    pub fn int(&self, key: &str) -> Option<i64> {
        self.0.get(key)?.as_i64()
    }
    #[allow(dead_code)]
    pub fn bool(&self, key: &str) -> Option<bool> {
        self.0.get(key)?.as_bool()
    }
}
