use crate::error::FavorError;

use super::Output;

pub struct HumanOutput;

impl HumanOutput {
    pub fn new() -> Self {
        Self
    }
}

impl Output for HumanOutput {
    fn status(&self, msg: &str) {
        eprintln!("  > {msg}");
    }

    fn success(&self, msg: &str) {
        eprintln!("  + {msg}");
    }

    fn warn(&self, msg: &str) {
        eprintln!("  ! {msg}");
    }

    fn error(&self, err: &FavorError) {
        eprintln!("  x {err}");
    }

    fn result_json(&self, data: &serde_json::Value) {
        if let Ok(json) = serde_json::to_string_pretty(data) {
            println!("{json}");
        }
    }

    fn table(&self, headers: &[&str], rows: &[Vec<String>]) {
        if headers.is_empty() {
            return;
        }
        let mut widths: Vec<usize> = headers.iter().map(|h| h.len()).collect();
        for row in rows {
            for (i, cell) in row.iter().enumerate() {
                if i < widths.len() {
                    widths[i] = widths[i].max(cell.len());
                }
            }
        }
        let header: String = headers
            .iter()
            .zip(&widths)
            .map(|(h, w)| format!("  {h:<w$}"))
            .collect::<Vec<_>>()
            .join("");
        eprintln!("{header}");
        let sep: String = widths.iter().map(|w| format!("  {}", "-".repeat(*w))).collect::<Vec<_>>().join("");
        eprintln!("{sep}");
        for row in rows {
            let line: String = row
                .iter()
                .zip(&widths)
                .map(|(cell, w)| format!("  {cell:<w$}"))
                .collect::<Vec<_>>()
                .join("");
            eprintln!("{line}");
        }
    }

}
