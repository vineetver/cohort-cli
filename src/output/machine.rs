use serde_json::json;

use crate::error::FavorError;

use super::Output;

pub struct MachineOutput;

impl MachineOutput {
    pub fn new() -> Self {
        Self
    }
}

impl Output for MachineOutput {
    fn status(&self, msg: &str) {
        eprintln!("{}", json!({"level": "info", "message": msg}));
    }

    fn success(&self, msg: &str) {
        eprintln!("{}", json!({"level": "info", "message": msg}));
    }

    fn warn(&self, msg: &str) {
        eprintln!("{}", json!({"level": "warn", "message": msg}));
    }

    fn error(&self, err: &FavorError) {
        eprintln!(
            "{}",
            json!({"error": err.code_name(), "message": err.to_string(), "exit_code": err.exit_code()})
        );
    }

    fn result_json(&self, data: &serde_json::Value) {
        if let Ok(json) = serde_json::to_string(data) {
            println!("{json}");
        }
    }

    fn table(&self, headers: &[&str], rows: &[Vec<String>]) {
        let result: Vec<serde_json::Value> = rows
            .iter()
            .map(|row| {
                let mut obj = serde_json::Map::new();
                for (i, cell) in row.iter().enumerate() {
                    let key = headers.get(i).unwrap_or(&"");
                    obj.insert(key.to_string(), serde_json::Value::String(cell.clone()));
                }
                serde_json::Value::Object(obj)
            })
            .collect();
        if let Ok(json) = serde_json::to_string(&result) {
            println!("{json}");
        }
    }

}
