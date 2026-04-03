//! Ingest: auto-normalize variant inputs to canonical parquet.
//!
//! Two modes:
//! - Happy path: `favor ingest <file>` — auto-detect, normalize, write parquet
//! - Escape hatch: `favor ingest <file> --emit-sql` — generate editable DuckDB script

pub mod columns;
pub mod sql_gen;
pub mod build_detect;
pub mod vcf;

use std::io::{BufRead, BufReader};
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::error::FavorError;

/// Detected input format — sum type, no invalid variant possible.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum InputFormat {
    Vcf,
    Tabular,
    Parquet,
}

/// Detected delimiter for tabular files.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[derive(Serialize)]
pub enum Delimiter {
    Tab,
    Comma,
    Space,
}

impl Delimiter {
    pub fn sql_literal(self) -> &'static str {
        match self {
            Delimiter::Tab => "'\\t'",
            Delimiter::Comma => "','",
            Delimiter::Space => "' '",
        }
    }
}

/// How the ingested data can join against annotations.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum JoinKey {
    ChromPosRefAlt,
    Rsid,
    ChromPos,
}

/// Detected genome build.
#[derive(Debug, Clone, Serialize)]
#[serde(tag = "verdict")]
pub enum BuildGuess {
    #[serde(rename = "hg38")]
    Hg38,
    #[serde(rename = "hg19")]
    Hg19 { match_rate_hg38: f64, match_rate_hg19: f64 },
    #[serde(rename = "unknown")]
    Unknown,
}

/// Detected coordinate base.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum CoordBase {
    OneBased,
    ZeroBased,
    Unknown,
}

/// A resolved column mapping.
#[derive(Debug, Clone, Serialize)]
pub struct ColumnMapping {
    pub input_name: String,
    pub canonical: &'static str,
}

/// Something the CLI couldn't resolve with >95% confidence.
#[derive(Debug, Clone, Serialize)]
pub struct Ambiguity {
    pub column: String,
    pub candidates: Vec<&'static str>,
    pub reason: &'static str,
}

/// What we learned from analyzing the input — query result, not mutable state (Cmd III).
#[derive(Debug, Serialize)]
pub struct Analysis {
    pub format: InputFormat,
    pub delimiter: Option<Delimiter>,
    pub raw_columns: Vec<String>,
    pub columns: Vec<ColumnMapping>,
    pub ambiguous: Vec<Ambiguity>,
    pub unmapped: Vec<String>,
    pub join_key: JoinKey,
    pub build_guess: BuildGuess,
    pub coord_base: CoordBase,
    /// The column name in the input that maps to "chromosome"
    pub chr_col: Option<String>,
    /// The column name in the input that maps to "position"
    pub pos_col: Option<String>,
    /// The column name in the input that maps to "ref"
    pub ref_col: Option<String>,
    /// The column name in the input that maps to "alt"
    pub alt_col: Option<String>,
    /// The column name in the input that maps to "rsid"
    pub rsid_col: Option<String>,
}

impl Analysis {
    /// Can we proceed automatically, or does the user/agent need to intervene?
    pub fn needs_intervention(&self) -> bool {
        !self.ambiguous.is_empty()
            || matches!(self.build_guess, BuildGuess::Hg19 { .. })
            || matches!(self.coord_base, CoordBase::ZeroBased)
    }

    /// Status string for machine output.
    pub fn status(&self) -> &'static str {
        if self.needs_intervention() { "needs_edit" } else { "ok" }
    }
}

/// Detect input format from file extension and first line.
/// This is the boundary — after this, InputFormat is a branded type.
pub fn detect_format(path: &Path) -> Result<(InputFormat, Option<Delimiter>), FavorError> {
    let name = path.file_name()
        .ok_or_else(|| FavorError::Input("No file name".into()))?
        .to_string_lossy()
        .to_lowercase();

    // VCF: extension-only detection
    if name.ends_with(".vcf") || name.ends_with(".vcf.gz")
        || name.ends_with(".vcf.bgz") || name.ends_with(".bcf") {
        return Ok((InputFormat::Vcf, None));
    }

    // Parquet: extension-only detection
    if name.ends_with(".parquet") {
        return Ok((InputFormat::Parquet, None));
    }

    // Tabular: detect delimiter from first line
    if name.ends_with(".tsv") || name.ends_with(".tsv.gz")
        || name.ends_with(".csv") || name.ends_with(".csv.gz")
        || name.ends_with(".txt") || name.ends_with(".txt.gz") {
        let delim = sniff_delimiter(path)?;
        return Ok((InputFormat::Tabular, Some(delim)));
    }

    // Unknown extension — try to sniff
    if path.is_file() {
        let delim = sniff_delimiter(path);
        if let Ok(d) = delim {
            return Ok((InputFormat::Tabular, Some(d)));
        }
    }

    Err(FavorError::Input(format!(
        "Cannot detect format for '{}'. Supported: .vcf, .vcf.gz, .tsv, .csv, .parquet, .txt",
        path.display()
    )))
}

/// Sniff delimiter from the first line of a text file.
fn sniff_delimiter(path: &Path) -> Result<Delimiter, FavorError> {
    let file = std::fs::File::open(path)
        .map_err(|e| FavorError::Input(format!("Cannot open '{}': {e}", path.display())))?;

    if path.to_string_lossy().ends_with(".gz") {
        // For .gz files, DuckDB handles decompression — default to tab
        return Ok(Delimiter::Tab);
    }

    let mut reader = BufReader::new(file);
    let mut first_line = String::new();
    reader.read_line(&mut first_line)
        .map_err(|e| FavorError::Input(format!("Cannot read '{}': {e}", path.display())))?;

    // Skip VCF header lines
    if first_line.starts_with("##") || first_line.starts_with("#CHROM") {
        return Ok(Delimiter::Tab); // VCF is always tab
    }

    let tab_count = first_line.matches('\t').count();
    let comma_count = first_line.matches(',').count();

    if tab_count >= 2 {
        Ok(Delimiter::Tab)
    } else if comma_count >= 2 {
        Ok(Delimiter::Comma)
    } else if first_line.split_whitespace().count() >= 3 {
        Ok(Delimiter::Space)
    } else {
        Ok(Delimiter::Tab) // default
    }
}

/// Read column headers from a tabular file (first non-comment line).
pub fn read_headers(path: &Path, delimiter: Delimiter) -> Result<Vec<String>, FavorError> {
    let file = std::fs::File::open(path)
        .map_err(|e| FavorError::Input(format!("Cannot open '{}': {e}", path.display())))?;

    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line = line.map_err(|e| FavorError::Input(format!("Read error: {e}")))?;
        let trimmed = line.trim();

        // Skip empty lines and VCF meta-headers
        if trimmed.is_empty() || trimmed.starts_with("##") {
            continue;
        }

        // Strip leading # from #CHROM or #chr
        let header_line = trimmed.strip_prefix('#').unwrap_or(trimmed);

        let cols: Vec<String> = match delimiter {
            Delimiter::Tab => header_line.split('\t').map(|s| s.trim().to_string()).collect(),
            Delimiter::Comma => header_line.split(',').map(|s| s.trim().to_string()).collect(),
            Delimiter::Space => header_line.split_whitespace().map(|s| s.to_string()).collect(),
        };

        return Ok(cols);
    }

    Err(FavorError::Input(format!("No header line found in '{}'", path.display())))
}

/// Full analysis: detect format, map columns, detect join key.
/// Build detection is separate (requires DuckDB + annotation parquets).
pub fn analyze(path: &Path) -> Result<Analysis, FavorError> {
    let (format, delimiter) = detect_format(path)?;

    match format {
        InputFormat::Vcf => Ok(Analysis {
            format,
            delimiter: None,
            raw_columns: vec!["CHROM".into(), "POS".into(), "ID".into(), "REF".into(), "ALT".into()],
            columns: vec![
                ColumnMapping { input_name: "CHROM".into(), canonical: "chromosome" },
                ColumnMapping { input_name: "POS".into(), canonical: "position" },
                ColumnMapping { input_name: "REF".into(), canonical: "ref" },
                ColumnMapping { input_name: "ALT".into(), canonical: "alt" },
                ColumnMapping { input_name: "ID".into(), canonical: "rsid" },
            ],
            ambiguous: vec![],
            unmapped: vec![],
            join_key: JoinKey::ChromPosRefAlt,
            build_guess: BuildGuess::Unknown,
            coord_base: CoordBase::OneBased, // VCF is always 1-based
            chr_col: Some("CHROM".into()),
            pos_col: Some("POS".into()),
            ref_col: Some("REF".into()),
            alt_col: Some("ALT".into()),
            rsid_col: Some("ID".into()),
        }),

        InputFormat::Parquet => {
            // For parquet, we'll inspect schema via DuckDB in the command handler
            Ok(Analysis {
                format,
                delimiter: None,
                raw_columns: vec![],
                columns: vec![],
                ambiguous: vec![],
                unmapped: vec![],
                join_key: JoinKey::ChromPosRefAlt,
                build_guess: BuildGuess::Unknown,
                coord_base: CoordBase::Unknown,
                chr_col: None,
                pos_col: None,
                ref_col: None,
                alt_col: None,
                rsid_col: None,
            })
        }

        InputFormat::Tabular => {
            let delim = delimiter.unwrap_or(Delimiter::Tab);
            let raw_columns = read_headers(path, delim)?;

            let (mapped, ambiguous, unmapped) = columns::map_columns(&raw_columns);

            // Extract key column names from the mapping
            let find_col = |canonical: &str| -> Option<String> {
                mapped.iter()
                    .find(|m| m.canonical == canonical)
                    .map(|m| m.input_name.clone())
            };

            let chr_col = find_col("chromosome");
            let pos_col = find_col("position");
            let ref_col = find_col("ref");
            let alt_col = find_col("alt");
            let rsid_col = find_col("rsid");

            // Determine join key from what columns are available
            let join_key = if chr_col.is_some() && pos_col.is_some()
                && ref_col.is_some() && alt_col.is_some() {
                JoinKey::ChromPosRefAlt
            } else if rsid_col.is_some() {
                JoinKey::Rsid
            } else if chr_col.is_some() && pos_col.is_some() {
                JoinKey::ChromPos
            } else {
                // Check if there's a composite variant_id column (chr:pos:ref:alt)
                let has_variant_id = raw_columns.iter()
                    .any(|c| {
                        let lower = c.to_lowercase();
                        lower == "variant_id" || lower == "varid" || lower == "variant"
                    });
                if has_variant_id {
                    JoinKey::ChromPosRefAlt // will be parsed from composite ID
                } else {
                    JoinKey::Rsid // fallback — hopefully rsid exists
                }
            };

            Ok(Analysis {
                format,
                delimiter: Some(delim),
                raw_columns,
                columns: mapped,
                ambiguous,
                unmapped,
                join_key,
                build_guess: BuildGuess::Unknown, // filled in by build_detect
                coord_base: CoordBase::Unknown,   // filled in by build_detect
                chr_col,
                pos_col,
                ref_col,
                alt_col,
                rsid_col,
            })
        }
    }
}
