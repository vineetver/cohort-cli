pub mod human;
pub mod machine;

use crate::error::FavorError;
use crate::mode::OutputMode;

/// All commands use this trait for output — never write stdout directly.
/// Progress bars are handled directly by download/verify code with byte-specific
/// formatting — the Output trait deliberately does not abstract over them (YAGNI).
pub trait Output {
    fn status(&self, msg: &str);
    fn success(&self, msg: &str);
    fn warn(&self, msg: &str);
    fn error(&self, err: &FavorError);
    fn result_json(&self, data: &serde_json::Value);
    fn table(&self, headers: &[&str], rows: &[Vec<String>]);
}

pub fn create(mode: &OutputMode) -> Box<dyn Output> {
    match mode {
        OutputMode::Human => Box::new(human::HumanOutput::new()),
        OutputMode::Machine => Box::new(machine::MachineOutput::new()),
    }
}
