use serde_json::json;

use crate::config::Config;
use crate::error::FavorError;
use crate::output::Output;

pub fn run(output: &dyn Output) -> Result<(), FavorError> {
    let config = Config::load_configured()?;

    let has_data = config.has_annotations();
    let has_tissue = config.has_tissue();

    let manifest = json!({
        "commands": [
            {"name": "ingest",    "status": "available",   "description": "Normalize VCF/TSV to canonical parquet with vid"},
            {"name": "annotate",  "status": if has_data { "available" } else { "unavailable" }, "requires": "annotation data", "reason": if !has_data { "run `favor setup` then `favor data pull`" } else { "" }},
            {"name": "enrich",    "status": if has_tissue { "available" } else { "unavailable" }, "requires": "tissue data", "reason": if !has_tissue { "no tissue/ directory found in root" } else { "" }},
            {"name": "interpret", "status": if has_data { "available" } else { "unavailable" }, "requires": "annotated variants"},
            {"name": "staar",     "status": if has_data { "available" } else { "unavailable" }, "requires": "genotypes + annotated variants"},
        ],
        "data": {
            "root": config.data.root_dir,
            "tier": config.data.tier.as_str(),
            "annotations_present": has_data,
            "tissue_present": has_tissue,
        }
    });

    output.result_json(&manifest);
    Ok(())
}
