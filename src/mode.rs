use std::io::IsTerminal;

use clap::ValueEnum;

#[derive(Clone, Debug, ValueEnum)]
pub enum Format {
    Auto,
    Json,
    Human,
}

#[derive(Clone, Debug, PartialEq)]
pub enum OutputMode {
    Human,
    Machine,
}

impl OutputMode {
    pub fn detect(format: &Format) -> Self {
        if std::env::var("FAVOR_MACHINE").is_ok_and(|v| !v.is_empty()) {
            return Self::Machine;
        }
        match format {
            Format::Json => Self::Machine,
            Format::Human => Self::Human,
            Format::Auto => {
                if std::io::stdout().is_terminal() {
                    Self::Human
                } else {
                    Self::Machine
                }
            }
        }
    }

    pub fn is_machine(&self) -> bool {
        *self == Self::Machine
    }
}
