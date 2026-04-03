use crate::cli::DataAction;
use crate::config::{Config, DirProbe, Environment, ResourceConfig, Tier};
use crate::error::FavorError;
use crate::mode::OutputMode;
use crate::output::Output;
use crate::packs::Pack;
use crate::resource::Resources;
use crate::tui;

pub fn run(
    output: &dyn Output,
    mode: &OutputMode,
    cli_environment: Option<String>,
    cli_memory_budget: Option<String>,
) -> Result<(), FavorError> {
    if mode.is_machine() {
        return Err(FavorError::Input(
            "setup requires interactive mode — run without --format json".to_string(),
        ));
    }

    // 1. Pick tier — returns Tier directly, no index mapping
    let tier = match tui::select_tier().map_err(|e| FavorError::Internal(e.into()))? {
        Some(t) => t,
        None => { output.warn("Setup cancelled"); return Ok(()); }
    };

    // 2. Pick root directory — browser shows live data probe
    let cwd = std::env::current_dir().unwrap_or_else(|_| Config::default_root_dir());
    let root = match tui::select_directory("Select FAVOR root directory", &cwd)
        .map_err(|e| FavorError::Internal(e.into()))?
    {
        Some(p) => p,
        None => { output.warn("Setup cancelled"); return Ok(()); }
    };

    // 3. Pick add-on packs — detect what's already installed on disk
    let optional_packs = Pack::optional();
    let installed_pack_ids: Vec<String> = optional_packs.iter()
        .filter(|p| p.tables.iter().any(|t| p.local_dir(&root).join(t).is_dir()))
        .map(|p| p.id.to_string())
        .collect();
    let selected_packs = tui::select_packs(&optional_packs, &installed_pack_ids)
        .map_err(|e| FavorError::Internal(e.into()))?
        .unwrap_or_default(); // Esc = skip, not cancel

    // 4. Environment selection — HPC or workstation?
    let environment = if let Some(env_str) = &cli_environment {
        Some(env_str.parse::<Environment>()?)
    } else {
        tui::select_environment().map_err(|e| FavorError::Internal(e.into()))?
    };

    // Validate HPC choice: check srun/sbatch availability
    if environment == Some(Environment::Hpc) {
        let has_srun = std::process::Command::new("which").arg("srun")
            .stdout(std::process::Stdio::null()).stderr(std::process::Stdio::null())
            .status().map(|s| s.success()).unwrap_or(false);
        let has_sbatch = std::process::Command::new("which").arg("sbatch")
            .stdout(std::process::Stdio::null()).stderr(std::process::Stdio::null())
            .status().map(|s| s.success()).unwrap_or(false);

        if !has_srun && !has_sbatch {
            output.warn("srun/sbatch not found in PATH — are you sure this is HPC?");
            output.warn("Continuing anyway. You can re-run `favor setup` to change.");
        }
    }

    // 5. Memory budget selection
    let res_detect = Resources::detect();
    let memory_budget = if let Some(budget) = &cli_memory_budget {
        // Validate the budget string
        if ResourceConfig::parse_memory_bytes(budget).is_none() {
            return Err(FavorError::Input(format!(
                "Cannot parse memory budget '{budget}'. Use format like '16GB', '64G', '8192MB'."
            )));
        }
        Some(budget.clone())
    } else {
        tui::select_memory_budget(&res_detect)
            .map_err(|e| FavorError::Internal(e.into()))?
    };

    // 6. Probe the chosen root to show what's already there
    let probe = DirProbe::scan(&root);

    output.status("FAVOR Configuration");
    output.status(&format!("  Root:       {}", root.display()));
    output.status(&format!("  Tier:       {tier}"));

    // Annotation status from probe
    let chrom_count = match tier {
        Tier::Full => probe.full_chroms,
        Tier::Base => probe.base_chroms,
    };
    if chrom_count == 24 {
        output.success(&format!("  {tier}: 24/24 chromosomes found"));
    } else if chrom_count > 0 {
        output.status(&format!("  {tier}: {chrom_count}/24 found, will download remaining"));
    } else {
        output.status(&format!("  {tier}: will download ({tier_size})", tier_size = tier.size_human()));
    }

    if !probe.tissue_tables.is_empty() {
        output.success(&format!("  Tissue: {} tables found", probe.tissue_tables.len()));
    }
    if !selected_packs.is_empty() {
        output.status(&format!("  Packs:  {}", selected_packs.join(", ")));
    }
    if let Some(env) = environment {
        output.status(&format!("  Env:    {}", env));
    }
    if let Some(budget) = &memory_budget {
        output.status(&format!("  Budget: {}", budget));
    }
    output.status(&format!("  System: {} memory, {} threads ({})",
        res_detect.memory_human(), res_detect.threads, res_detect.environment()));

    // 7. Save config — preserve existing resource customizations, overlay new ones
    let mut existing_resources = Config::load()
        .map(|c| c.resources)
        .unwrap_or_default();

    // Update environment and budget from this setup run
    if environment.is_some() {
        existing_resources.environment = environment;
    }
    if memory_budget.is_some() {
        existing_resources.memory_budget = memory_budget;
    }

    let config = Config {
        data: crate::config::DataConfig {
            tier,
            root_dir: root.to_string_lossy().to_string(),
            packs: selected_packs.clone(),
        },
        resources: existing_resources,
    };
    config.save()?;
    output.success(&format!("Saved to {}", Config::config_path().display()));

    // 8. Download annotations (skips if already complete)
    output.status("Checking annotations...");
    if let Err(e) = super::data::run(
        DataAction::Pull { full: tier == Tier::Full, dry_run: false, yes: true, pack: None },
        output,
    ) {
        output.warn("Config saved but annotation download incomplete. Run `favor data pull` to retry.");
        return Err(e);
    }

    // 9. Download always-installed packs
    for pack in Pack::required() {
        if let Err(e) = super::data::pull_pack(pack.id, false, true, output) {
            output.warn(&format!("{} failed: {e}. Run `favor data pull --pack {}` to retry.",
                pack.name, pack.id));
        }
    }

    // 10. Download selected packs
    for pack_id in &selected_packs {
        if let Err(e) = super::data::pull_pack(pack_id, false, true, output) {
            output.warn(&format!("Pack '{pack_id}' failed: {e}. Run `favor data pull --pack {pack_id}` to retry."));
        }
    }

    output.success("Setup complete. Run: favor annotate <input.vcf>");
    Ok(())
}
