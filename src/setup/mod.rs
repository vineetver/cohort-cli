use std::io;
use std::path::{Path, PathBuf};

use crossterm::event::{self, Event, KeyCode, KeyEventKind};
use crossterm::terminal::{disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen};
use crossterm::ExecutableCommand;
use ratatui::prelude::*;
use ratatui::widgets::{Block, Borders, List, ListItem, ListState, Paragraph, Padding};
use serde_json::json;

use crate::cli::DataAction;
use crate::config::{Config, DirProbe, Environment, ProbeStatus, ResourceConfig, Tier};
use crate::data::Pack;
use crate::error::FavorError;
use crate::output::OutputMode;
use crate::output::Output;
use crate::resource::Resources;

// ---------------------------------------------------------------------------
// Init command (was commands/init.rs)
// ---------------------------------------------------------------------------

pub fn init(
    path: Option<PathBuf>,
    force: bool,
    out: &dyn Output,
    mode: &OutputMode,
) -> Result<(), FavorError> {
    let config = Config::load_configured()?;

    let project_dir = match path {
        Some(p) => {
            std::fs::create_dir_all(&p)?;
            p.canonicalize().map_err(|e| {
                FavorError::Input(format!("Cannot resolve '{}': {e}", p.display()))
            })?
        }
        None => std::env::current_dir()?,
    };

    let claude_path = project_dir.join("CLAUDE.md");
    let codex_dir = project_dir.join(".codex");
    let codex_path = codex_dir.join("instructions.md");

    let is_refresh = claude_path.exists();

    if is_refresh && !force && !mode.is_machine() {
        out.warn(
            "CLAUDE.md already exists. Use --force to overwrite, or re-run to refresh after installing new packs.",
        );
        return Ok(());
    }

    let probe = DirProbe::scan(&config.root_dir());
    let content = render_init_template(&config, &probe);

    std::fs::write(&claude_path, &content)?;
    std::fs::create_dir_all(&codex_dir)?;
    std::fs::write(&codex_path, &content)?;

    let action = if is_refresh { "Refreshed" } else { "Initialized" };
    out.success(&format!("{action} agent context in {}", project_dir.display()));
    out.status("  CLAUDE.md (Claude Code)");
    out.status("  .codex/instructions.md (Codex)");
    out.status(&format!(
        "Open your coding agent in {} and ask a research question.",
        project_dir.display()
    ));

    out.result_json(&json!({
        "project_dir": project_dir.to_string_lossy(),
        "files": [claude_path.to_string_lossy(), codex_path.to_string_lossy()],
        "tier": config.data.tier.as_str(),
        "packs": config.data.packs,
    }));

    Ok(())
}

fn render_init_template(config: &Config, probe: &DirProbe) -> String {
    let pack_list = build_pack_list(config, probe);

    let environment = config
        .resources
        .environment
        .map(|e| e.as_str().to_string())
        .unwrap_or_else(|| "auto-detect".to_string());

    let memory_budget = config
        .resources
        .memory_budget
        .clone()
        .unwrap_or_else(|| "auto-detect".to_string());

    let chrom_count = match config.data.tier {
        Tier::Base => probe.base_chroms,
        Tier::Full => probe.full_chroms,
    };

    INIT_TEMPLATE
        .replace("{root_dir}", &config.data.root_dir)
        .replace("{tier}", config.data.tier.as_str())
        .replace("{tier_size}", config.data.tier.size_human())
        .replace("{chrom_count}", &chrom_count.to_string())
        .replace("{environment}", &environment)
        .replace("{memory_budget}", &memory_budget)
        .replace("{pack_list}", &pack_list)
}

fn build_pack_list(config: &Config, probe: &DirProbe) -> String {
    let mut lines = Vec::new();
    for pack in Pack::all() {
        if pack.always_installed {
            continue;
        }
        let installed = config.data.packs.contains(&pack.id.to_string())
            || pack
                .tables
                .iter()
                .any(|t| probe.tissue_tables.contains(&t.to_string()));
        if installed {
            lines.push(format!(
                "- **{}** ({}): {}",
                pack.id, pack.size_human, pack.description
            ));
        }
    }
    if lines.is_empty() {
        lines.push(
            "- No optional packs. Run `favor data pull --pack <id>` then `favor init --force`."
                .to_string(),
        );
    }
    lines.join("\n")
}

const INIT_TEMPLATE: &str = r#"# FAVOR Genomic Analysis

Always pass `--format json` to every favor command. Use `--dry-run` before heavy computation.

## System

- Data root: {root_dir}
- Tier: {tier} ({tier_size}, {chrom_count}/24 chromosomes)
- Environment: {environment}
- Memory budget: {memory_budget}

## Installed packs

{pack_list}

## Pipeline

```
variant file -> favor ingest -> favor annotate -> favor enrich --tissue X
                                                -> favor staar (rare-variant association)
                                                -> favor interpret (variant-to-gene)
```

## Commands

| Command | What it does |
|---------|-------------|
| `favor ingest <file>` | Normalize VCF/TSV/CSV to parquet with variant ID |
| `favor annotate <file>` | Add CADD, ClinVar, gnomAD, REVEL, aPC scores, regulatory marks |
| `favor enrich <file> --tissue <name>` | Add tissue eQTL, ChromBPNet, enhancer-gene links |
| `favor staar --genotypes <vcf> --phenotype <tsv> --trait-name <col> --annotations <parquet>` | Rare-variant burden test (STAAR-O) |
| `favor schema [table]` | Show column names and types |
| `favor manifest` | List installed data and available commands |

## STAAR usage

```bash
favor staar --dry-run --format json \
  --genotypes cohort.vcf.gz \
  --phenotype pheno.tsv \
  --trait-name LDL \
  --covariates age,sex,PC1,PC2 \
  --annotations annotated.parquet \
  --masks coding
```

Use `--dry-run` first to check memory. On HPC: `srun --mem=64G -c 8 favor staar ...`

## Interpreting results

- **CADD phred > 20**: top 1% most deleterious variants genome-wide
- **REVEL > 0.5**: likely pathogenic missense
- **STAAR-O < 2.5e-6**: genome-wide significant (Bonferroni for ~20K genes)
- **pLoF mask significant**: protein-destroying variants cause the trait (strongest evidence)
- **Synonymous mask significant**: suspicious — likely confounding, not biology

STAAR runs 6 tests (Burden, SKAT, ACAT-V x two beta weights) across 11 annotation channels:
- Burden: all variants push trait same direction
- SKAT: mixed-effect variants
- ACAT-V: one or two variants drive everything

## Querying parquet output

Use DuckDB:

```sql
-- Damaging variants
SELECT vid, gencode.genes[1] AS gene, main.cadd.phred
FROM 'annotated.parquet' WHERE main.cadd.phred > 20;

-- Significant genes
SELECT * FROM 'staar_results/coding_pLoF_missense.parquet'
WHERE "STAAR-O" < 2.5e-6 ORDER BY "STAAR-O";

-- Tissue enrichment
SELECT a.vid, e.* FROM 'enriched/annotated.parquet' a
JOIN 'enriched/eqtl.parquet' e ON a.vid = e.vid
WHERE e.tissue_name LIKE '%Liver%';
```

## Output conventions

- Parquet files: main results
- `.meta.json`: parameters and counts
- stdout (--format json): structured result summary
- stderr: progress/status messages
- Exit codes: 0=ok, 1=input, 2=data missing, 3=resource, 4=analysis
"#;

// ---------------------------------------------------------------------------
// Setup command (was commands/setup.rs)
// ---------------------------------------------------------------------------

pub fn setup(
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
    let tier = match select_tier().map_err(|e| FavorError::Internal(e.into()))? {
        Some(t) => t,
        None => { output.warn("Setup cancelled"); return Ok(()); }
    };

    // 2. Pick root directory — browser shows live data probe
    let cwd = std::env::current_dir().unwrap_or_else(|_| Config::default_root_dir());
    let root = match select_directory("Select FAVOR root directory", &cwd)
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
    let selected_packs = select_packs(&optional_packs, &installed_pack_ids)
        .map_err(|e| FavorError::Internal(e.into()))?
        .unwrap_or_default(); // Esc = skip, not cancel

    // 4. Environment selection — HPC or workstation?
    let environment = if let Some(env_str) = &cli_environment {
        Some(env_str.parse::<Environment>()?)
    } else {
        select_environment().map_err(|e| FavorError::Internal(e.into()))?
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
        select_memory_budget(&res_detect)
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
    if let Err(e) = crate::data::transfer::run(
        DataAction::Pull { full: tier == Tier::Full, dry_run: false, yes: true, pack: None },
        output,
    ) {
        output.warn("Config saved but annotation download incomplete. Run `favor data pull` to retry.");
        return Err(e);
    }

    // 9. Download always-installed packs
    for pack in Pack::required() {
        if let Err(e) = crate::data::transfer::pull_pack(pack.id, false, true, output) {
            output.warn(&format!("{} failed: {e}. Run `favor data pull --pack {}` to retry.",
                pack.name, pack.id));
        }
    }

    // 10. Download selected packs
    for pack_id in &selected_packs {
        if let Err(e) = crate::data::transfer::pull_pack(pack_id, false, true, output) {
            output.warn(&format!("Pack '{pack_id}' failed: {e}. Run `favor data pull --pack {pack_id}` to retry."));
        }
    }

    output.success("Setup complete. Run: favor annotate <input.vcf>");
    Ok(())
}

// ---------------------------------------------------------------------------
// Uninstall command (was commands/uninstall.rs)
// ---------------------------------------------------------------------------

pub fn uninstall(out: &dyn Output) -> Result<(), FavorError> {
    let binary = std::env::current_exe().unwrap_or_default();
    let config_dir = Config::config_dir();

    out.status("This will remove:");
    out.status(&format!("  Binary: {}", binary.display()));
    out.status(&format!("  Config: {}", config_dir.display()));
    out.warn("Data packs at your configured root directory will NOT be deleted.");

    // Remove config directory
    if config_dir.exists() {
        std::fs::remove_dir_all(&config_dir)?;
        out.status("  Removed config directory");
    }

    // Remove PATH entry from shell rc
    let home = dirs::home_dir().unwrap_or_default();
    let install_dir = binary.parent().unwrap_or(std::path::Path::new(""));
    let install_str = install_dir.to_string_lossy();

    for rc in &[".bashrc", ".zshrc"] {
        let rc_path = home.join(rc);
        if rc_path.exists() {
            if let Ok(content) = std::fs::read_to_string(&rc_path) {
                let filtered: Vec<&str> = content.lines()
                    .filter(|line| !line.contains(&*install_str) || !line.contains("PATH"))
                    .collect();
                if filtered.len() < content.lines().count() {
                    let _ = std::fs::write(&rc_path, filtered.join("\n") + "\n");
                    out.status(&format!("  Cleaned PATH from {}", rc));
                }
            }
        }
    }

    // Remove binary last (we're running it)
    if binary.exists() {
        std::fs::remove_file(&binary)?;
    }

    out.success("Uninstalled. Data packs remain at your configured root directory.");
    Ok(())
}

// ---------------------------------------------------------------------------
// TUI widgets (was tui.rs)
// ---------------------------------------------------------------------------

/// RAII guard for raw mode + alternate screen.
/// Terminal state is always restored, even on panic.
struct TermGuard;

impl TermGuard {
    fn enter() -> io::Result<Self> {
        enable_raw_mode()?;
        io::stdout().execute(EnterAlternateScreen)?;
        Ok(TermGuard)
    }
}

impl Drop for TermGuard {
    fn drop(&mut self) {
        let _ = io::stdout().execute(LeaveAlternateScreen);
        let _ = disable_raw_mode();
    }
}

struct TierOption {
    tier: Tier,
    summary: &'static str,
    details: &'static [&'static str],
}

const TIERS: &[TierOption] = &[
    TierOption {
        tier: Tier::Base,
        summary: "Curated annotations for most analyses",
        details: &[
            "Genes        GENCODE consequence, transcripts, hgvsc/hgvsp",
            "Clinical     ClinVar significance, disease, review status",
            "Frequency    gnomAD genome+exome (AF only), TOPMed, 1000G 6 pops",
            "Coding       CADD, REVEL, SpliceAI, AlphaMissense, MaveDB",
            "Noncoding    LINSIGHT, FATHMM-XF, GPN-MSA, JARVIS, ReMM, ncER",
            "Conservation PhyloP, PhastCons (3 levels each), GERP, B-statistic",
            "Constraint   gnomAD constraint score + phred",
            "Integrative  13 aPC annotation principal component scores",
            "Regulatory   MACIE, cV2F, cCRE, GeneHancer, CAGE, super enhancers",
            "dbNSFP       15 of 30 predictors (REVEL, MetaSVM, BayesDel, ...)",
            "Other        Mutation rate, distance to TSS/TSE",
            "",
            "Not in base: UCSC, RefSeq, COSMIC, ChromHMM, ENCODE histone marks,",
            "  full gnomAD pops, 15 extra dbNSFP, mappability, variant density",
        ],
    },
    TierOption {
        tier: Tier::Full,
        summary: "Complete FAVOR database - all annotations",
        details: &[
            "Genes        GENCODE + UCSC + RefSeq",
            "Clinical     ClinVar, COSMIC cancer mutations",
            "Frequency    gnomAD 9 pops + sex-stratified + FAF, TOPMed, 1000G",
            "Coding       CADD, REVEL, SpliceAI, AlphaMissense, MaveDB",
            "Noncoding    LINSIGHT, FATHMM-XF, GPN-MSA, JARVIS, ReMM, ncER,",
            "             ncBoost, PGBoost, FunSeq, ALoFT",
            "Conservation PhyloP, PhastCons (3 levels each), GERP, B-statistic",
            "Constraint   gnomAD constraint score + phred",
            "Integrative  13 aPC annotation principal component scores",
            "Regulatory   MACIE, cV2F, cCRE, GeneHancer, CAGE, super enhancers",
            "Epigenomics  ChromHMM 25 states, ENCODE 13 histone marks (raw+phred)",
            "dbNSFP       All 30 predictors",
            "Sequence     Variant density, ReMap, GC/CpG, mappability (4 levels)",
            "Other        Mutation rate, distance, recombination, nucleotide div",
        ],
    },
];

/// Interactive tier selector with live preview panel.
/// Returns the selected Tier, or None if user cancelled.
fn select_tier() -> io::Result<Option<Tier>> {
    let _guard = TermGuard::enter()?;

    let backend = CrosstermBackend::new(io::stdout());
    let mut terminal = Terminal::new(backend)?;
    let mut selected: usize = 0;

    loop {
        terminal.draw(|frame| draw_tier(frame, selected))?;

        if let Event::Key(key) = event::read()? {
            if key.kind != KeyEventKind::Press { continue; }
            match key.code {
                KeyCode::Up | KeyCode::Char('k') => {
                    if selected > 0 { selected -= 1; }
                }
                KeyCode::Down | KeyCode::Char('j') | KeyCode::Tab => {
                    if selected < TIERS.len() - 1 { selected += 1; }
                }
                KeyCode::Enter => return Ok(Some(TIERS[selected].tier)),
                KeyCode::Esc | KeyCode::Char('q') => return Ok(None),
                _ => {}
            }
        }
    }
}

fn draw_tier(frame: &mut Frame, selected: usize) {
    let area = frame.area();

    let layout = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3),
            Constraint::Length(5),
            Constraint::Min(10),
            Constraint::Length(2),
        ])
        .split(area);

    let title = Paragraph::new("  FAVOR Setup — Annotation Tier")
        .style(Style::default().fg(Color::Cyan).bold())
        .block(Block::default().padding(Padding::top(1)));
    frame.render_widget(title, layout[0]);

    let mut selector_lines = Vec::new();
    for (i, tier_opt) in TIERS.iter().enumerate() {
        let marker = if i == selected { " > " } else { "   " };
        let style = if i == selected {
            Style::default().fg(Color::Cyan).bold()
        } else {
            Style::default().fg(Color::DarkGray)
        };
        selector_lines.push(Line::from(vec![
            Span::raw(marker),
            Span::styled(format!("{:<6}", tier_opt.tier.as_str()), style),
            Span::styled(
                format!("{:<12}", tier_opt.tier.size_human()),
                if i == selected { Style::default().fg(Color::Yellow) }
                else { Style::default().fg(Color::DarkGray) },
            ),
            Span::styled(
                tier_opt.summary,
                if i == selected { Style::default().fg(Color::White) }
                else { Style::default().fg(Color::DarkGray) },
            ),
        ]));
    }
    let selector = Paragraph::new(selector_lines)
        .block(Block::default()
            .borders(Borders::ALL)
            .title(" Select annotation tier ")
            .border_style(Style::default().fg(Color::Cyan)));
    frame.render_widget(selector, layout[1]);

    let tier_opt = &TIERS[selected];
    let detail_lines: Vec<Line> = tier_opt.details.iter().map(|l| {
        if l.is_empty() {
            Line::from("")
        } else if l.starts_with("Not in") {
            Line::from(Span::styled(
                format!("  {l}"),
                Style::default().fg(Color::DarkGray).italic(),
            ))
        } else {
            let trimmed = l.trim_start();
            if let Some(pos) = trimmed.find("  ") {
                let (label, rest) = trimmed.split_at(pos);
                Line::from(vec![
                    Span::styled(format!("  {:<13}", label), Style::default().fg(Color::Yellow)),
                    Span::styled(rest.trim_start(), Style::default().fg(Color::White)),
                ])
            } else {
                Line::from(Span::styled(format!("  {trimmed}"), Style::default().fg(Color::White)))
            }
        }
    }).collect();

    let detail_title = format!(" {} - {} ", tier_opt.tier.as_str(), tier_opt.tier.size_human());
    let detail = Paragraph::new(detail_lines)
        .block(Block::default()
            .borders(Borders::ALL)
            .title(detail_title)
            .border_style(Style::default().fg(Color::DarkGray)));
    frame.render_widget(detail, layout[2]);

    let help = Paragraph::new("  up/down navigate    enter select    esc cancel")
        .style(Style::default().fg(Color::DarkGray));
    frame.render_widget(help, layout[3]);
}

const SELECT_ENTRY: &str = "[ Use this directory ]";

/// Directory browser state — everything needed to render one frame.
struct DirBrowserState {
    prompt: String,
    current_dir: PathBuf,
    entries: Vec<String>,
    list_state: ListState,
    typing_path: bool,
    input_buf: String,
    probe: DirProbe,
}

impl DirBrowserState {
    fn new(prompt: &str, start: &Path) -> Self {
        let current_dir = if start.is_dir() {
            start.to_path_buf()
        } else {
            start.parent().map(|p| p.to_path_buf()).unwrap_or_else(|| PathBuf::from("/"))
        };
        let entries = list_dirs(&current_dir);
        let mut list_state = ListState::default();
        if !entries.is_empty() { list_state.select(Some(0)); }
        let probe = DirProbe::scan(&current_dir);

        Self {
            prompt: prompt.to_string(),
            current_dir,
            entries,
            list_state,
            typing_path: false,
            input_buf: String::new(),
            probe,
        }
    }

    fn navigate_to(&mut self, dir: PathBuf) {
        self.current_dir = dir;
        self.entries = list_dirs(&self.current_dir);
        self.list_state.select(if self.entries.is_empty() { None } else { Some(0) });
        self.input_buf = self.current_dir.to_string_lossy().to_string();
        self.probe = DirProbe::scan(&self.current_dir);
    }

    fn go_parent(&mut self) {
        if let Some(parent) = self.current_dir.parent() {
            self.navigate_to(parent.to_path_buf());
        }
    }

    fn enter_selected(&mut self) -> Option<PathBuf> {
        let i = self.list_state.selected()?;
        let name = &self.entries[i];

        if name == SELECT_ENTRY {
            return Some(self.current_dir.clone());
        }
        if name == ".." {
            self.go_parent();
            return None;
        }

        let target = self.current_dir.join(name);
        if target.is_dir() {
            self.navigate_to(target);
        }
        None
    }
}

/// Interactive directory browser with live FAVOR data probe panel.
/// Returns the selected path, or None if cancelled.
fn select_directory(prompt: &str, default_path: &Path) -> io::Result<Option<PathBuf>> {
    let _guard = TermGuard::enter()?;

    let backend = CrosstermBackend::new(io::stdout());
    let mut terminal = Terminal::new(backend)?;
    let mut state = DirBrowserState::new(prompt, default_path);

    loop {
        terminal.draw(|frame| draw_dir_browser(frame, &mut state))?;

        if let Event::Key(key) = event::read()? {
            if key.kind != KeyEventKind::Press { continue; }

            if state.typing_path {
                match key.code {
                    KeyCode::Enter => {
                        let p = PathBuf::from(&state.input_buf);
                        if p.is_dir() {
                            state.navigate_to(p);
                            state.typing_path = false;
                        }
                        // Invalid path: stay in typing mode (red highlight tells them why)
                    }
                    KeyCode::Esc => {
                        state.input_buf = state.current_dir.to_string_lossy().to_string();
                        state.typing_path = false;
                    }
                    KeyCode::Backspace => { state.input_buf.pop(); }
                    KeyCode::Char(c) => { state.input_buf.push(c); }
                    KeyCode::Tab => {
                        // Tab-complete: find first matching directory
                        if let Some(completed) = tab_complete(&state.input_buf) {
                            state.input_buf = completed;
                        }
                    }
                    _ => {}
                }
                // Eagerly update probe while typing — if path is a valid dir,
                // show what's there BEFORE the user hits enter
                let typed = Path::new(&state.input_buf);
                if typed.is_dir() {
                    state.probe = DirProbe::scan(typed);
                }
                continue;
            }

            match key.code {
                KeyCode::Up | KeyCode::Char('k') => {
                    if let Some(i) = state.list_state.selected() {
                        if i > 0 { state.list_state.select(Some(i - 1)); }
                    }
                }
                KeyCode::Down | KeyCode::Char('j') => {
                    if let Some(i) = state.list_state.selected() {
                        if i + 1 < state.entries.len() {
                            state.list_state.select(Some(i + 1));
                        }
                    }
                }
                KeyCode::Right | KeyCode::Enter => {
                    if let Some(selected) = state.enter_selected() {
                        return Ok(Some(selected));
                    }
                }
                KeyCode::Left | KeyCode::Backspace => {
                    state.go_parent();
                }
                KeyCode::Char(' ') => {
                    // Space always confirms current directory (from any position)
                    return Ok(Some(state.current_dir.clone()));
                }
                KeyCode::Char('c') => {
                    // 'c' = create directory at current path (useful for empty root)
                    // Silently no-op if it already exists
                    let _ = std::fs::create_dir_all(&state.current_dir);
                    return Ok(Some(state.current_dir.clone()));
                }
                KeyCode::Char('/') | KeyCode::Char('g') => {
                    state.typing_path = true;
                    state.input_buf = state.current_dir.to_string_lossy().to_string();
                }
                KeyCode::Esc | KeyCode::Char('q') => return Ok(None),
                _ => {}
            }
        }
    }
}

/// Tab-complete a partial path. Returns the completed path if exactly one match,
/// or the longest common prefix if multiple matches.
fn tab_complete(partial: &str) -> Option<String> {
    let path = Path::new(partial);
    let (parent, prefix) = if partial.ends_with('/') {
        // User typed a full directory — list its contents
        if path.is_dir() {
            return list_first_child(path);
        }
        return None;
    } else {
        // Split into parent dir + partial filename
        let parent = path.parent()?;
        let prefix = path.file_name()?.to_string_lossy().to_string();
        (parent, prefix)
    };

    if !parent.is_dir() { return None; }

    let matches: Vec<String> = std::fs::read_dir(parent).ok()?
        .flatten()
        .filter(|e| e.path().is_dir())
        .filter_map(|e| {
            let name = e.file_name().into_string().ok()?;
            if name.starts_with(&prefix) { Some(name) } else { None }
        })
        .collect();

    match matches.len() {
        0 => None,
        1 => {
            let completed = parent.join(&matches[0]);
            Some(format!("{}/", completed.to_string_lossy()))
        }
        _ => {
            // Find longest common prefix among matches
            let lcp = longest_common_prefix(&matches);
            let completed = parent.join(&lcp);
            Some(completed.to_string_lossy().to_string())
        }
    }
}

fn list_first_child(dir: &Path) -> Option<String> {
    let mut entries: Vec<String> = std::fs::read_dir(dir).ok()?
        .flatten()
        .filter(|e| e.path().is_dir())
        .filter_map(|e| e.file_name().into_string().ok())
        .collect();
    entries.sort();
    if entries.len() == 1 {
        Some(format!("{}/{}/", dir.to_string_lossy(), entries[0]))
    } else {
        None
    }
}

fn longest_common_prefix(strings: &[String]) -> String {
    if strings.is_empty() { return String::new(); }
    let first = &strings[0];
    let mut len = first.len();
    for s in &strings[1..] {
        len = len.min(s.len());
        for (i, (a, b)) in first.bytes().zip(s.bytes()).enumerate() {
            if a != b { len = len.min(i); break; }
        }
    }
    first[..len].to_string()
}

fn list_dirs(dir: &Path) -> Vec<String> {
    let mut dirs = vec![SELECT_ENTRY.to_string(), "..".to_string()];
    if let Ok(read) = std::fs::read_dir(dir) {
        let mut entries: Vec<String> = read
            .flatten()
            .filter(|e| e.path().is_dir())
            .filter_map(|e| e.file_name().into_string().ok())
            .collect();
        entries.sort();
        dirs.extend(entries);
    }
    dirs
}

fn draw_dir_browser(frame: &mut Frame, state: &mut DirBrowserState) {
    let area = frame.area();

    // Split: left (browser) | right (probe panel)
    let h_split = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(60), Constraint::Percentage(40)])
        .split(area);

    let left_layout = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3),  // title
            Constraint::Length(3),  // current path / input
            Constraint::Min(6),    // directory list
            Constraint::Length(2), // help
        ])
        .split(h_split[0]);

    // Title
    let title = Paragraph::new(format!("  {}", state.prompt))
        .style(Style::default().fg(Color::Cyan).bold())
        .block(Block::default().padding(Padding::top(1)));
    frame.render_widget(title, left_layout[0]);

    // Path bar — validate on every keystroke (F4)
    if state.typing_path {
        let path_valid = Path::new(&state.input_buf).is_dir();
        let path_color = if path_valid { Color::Green } else { Color::Red };
        let path_line = Line::from(vec![
            Span::styled(" > ", Style::default().fg(Color::Yellow)),
            Span::styled(&state.input_buf, Style::default().fg(path_color)),
            Span::styled("_", Style::default().fg(Color::Cyan)),
        ]);
        let hint = if path_valid { " Valid directory — enter to go " }
                   else { " Not a directory " };
        let path_block = Paragraph::new(path_line)
            .block(Block::default()
                .borders(Borders::ALL)
                .title(hint)
                .border_style(Style::default().fg(if path_valid { Color::Green } else { Color::Yellow })));
        frame.render_widget(path_block, left_layout[1]);
    } else {
        let path_line = Line::from(vec![
            Span::styled(" ", Style::default()),
            Span::styled(
                state.current_dir.to_string_lossy().to_string(),
                Style::default().fg(Color::White).bold(),
            ),
        ]);
        let path_block = Paragraph::new(path_line)
            .block(Block::default()
                .borders(Borders::ALL)
                .title(" Current directory ")
                .border_style(Style::default().fg(Color::Cyan)));
        frame.render_widget(path_block, left_layout[1]);
    }

    // Directory listing
    let items: Vec<ListItem> = state.entries.iter().enumerate().map(|(i, name)| {
        let is_selected = state.list_state.selected() == Some(i);
        let (prefix, style) = match name.as_str() {
            s if s == SELECT_ENTRY => (" ", if is_selected {
                Style::default().fg(Color::Green).bold()
            } else {
                Style::default().fg(Color::Green)
            }),
            ".." => (" ..", if is_selected {
                Style::default().fg(Color::Cyan).bold()
            } else {
                Style::default().fg(Color::DarkGray)
            }),
            s if s.starts_with('.') => (&name[..], if is_selected {
                Style::default().fg(Color::DarkGray).bold()
            } else {
                Style::default().fg(Color::DarkGray)
            }),
            _ => (&name[..], if is_selected {
                Style::default().fg(Color::Cyan).bold()
            } else {
                Style::default().fg(Color::White)
            }),
        };
        let display = if name == SELECT_ENTRY || name == ".." {
            format!(" {prefix}")
        } else {
            format!(" /{name}")
        };
        ListItem::new(display).style(style)
    }).collect();

    let list = List::new(items)
        .block(Block::default()
            .borders(Borders::ALL)
            .title(" Directories ")
            .border_style(Style::default().fg(Color::DarkGray)))
        .highlight_style(Style::default().bg(Color::DarkGray).fg(Color::White))
        .highlight_symbol(" > ");
    frame.render_stateful_widget(list, left_layout[2], &mut state.list_state);

    // Help
    let help_text = if state.typing_path {
        "  type path    enter go    esc cancel"
    } else {
        "  enter open    space select    / type path    esc cancel"
    };
    let help = Paragraph::new(help_text)
        .style(Style::default().fg(Color::DarkGray));
    frame.render_widget(help, left_layout[3]);

    let right_layout = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3),  // header
            Constraint::Min(6),    // probe details
            Constraint::Length(3), // verdict
        ])
        .split(h_split[1]);

    // Header
    let probe_header = if state.probe.has_any_data() {
        Paragraph::new("  FAVOR data detected")
            .style(Style::default().fg(Color::Green).bold())
            .block(Block::default().padding(Padding::top(1)))
    } else {
        Paragraph::new("  No FAVOR data at this path")
            .style(Style::default().fg(Color::DarkGray))
            .block(Block::default().padding(Padding::top(1)))
    };
    frame.render_widget(probe_header, right_layout[0]);

    // Probe detail lines
    let summary = state.probe.summary_lines();
    let probe_lines: Vec<Line> = summary.iter().map(|(text, status)| {
        let (icon, color) = match status {
            ProbeStatus::Good    => ("  +  ", Color::Green),
            ProbeStatus::Partial => ("  ~  ", Color::Yellow),
            ProbeStatus::Missing => ("  -  ", Color::DarkGray),
            ProbeStatus::Info    => ("     ", Color::DarkGray),
        };
        Line::from(vec![
            Span::styled(icon, Style::default().fg(color)),
            Span::styled(text, Style::default().fg(color)),
        ])
    }).collect();

    let probe_detail = Paragraph::new(probe_lines)
        .block(Block::default()
            .borders(Borders::ALL)
            .title(" Data scan ")
            .border_style(if state.probe.has_any_data() {
                Style::default().fg(Color::Green)
            } else {
                Style::default().fg(Color::DarkGray)
            }));
    frame.render_widget(probe_detail, right_layout[1]);

    // Verdict
    let verdict = if state.probe.has_any_data() {
        let tier = state.probe.detected_tier()
            .map(|t| t.as_str())
            .unwrap_or("?");
        Paragraph::new(Line::from(vec![
            Span::styled("  Ready — ", Style::default().fg(Color::Green)),
            Span::styled(format!("{tier} tier detected"), Style::default().fg(Color::Green).bold()),
        ]))
    } else {
        Paragraph::new(Line::from(vec![
            Span::styled("  Will download data after setup", Style::default().fg(Color::Yellow)),
        ]))
    };
    frame.render_widget(verdict, right_layout[2]);
}

/// Interactive multi-select for add-on packs.
/// `installed` contains pack IDs already present on disk — shown as pre-checked with "(installed)".
/// Returns list of selected pack IDs, or None if cancelled.
fn select_packs(packs: &[&Pack], installed: &[String]) -> io::Result<Option<Vec<String>>> {
    let _guard = TermGuard::enter()?;

    let backend = CrosstermBackend::new(io::stdout());
    let mut terminal = Terminal::new(backend)?;

    let mut cursor: usize = 0;
    // Pre-check packs that are already installed
    let mut checked: Vec<bool> = packs.iter()
        .map(|p| installed.iter().any(|id| id == p.id))
        .collect();

    loop {
        let draw_checked = &checked;
        terminal.draw(|frame| {
            draw_pack_selector(frame, packs, draw_checked, installed, cursor);
        })?;

        if let Event::Key(key) = event::read()? {
            if key.kind != KeyEventKind::Press { continue; }
            match key.code {
                KeyCode::Up | KeyCode::Char('k') => {
                    if cursor > 0 { cursor -= 1; }
                }
                KeyCode::Down | KeyCode::Char('j') | KeyCode::Tab => {
                    if cursor + 1 < packs.len() { cursor += 1; }
                }
                KeyCode::Char(' ') => {
                    checked[cursor] = !checked[cursor];
                }
                KeyCode::Char('a') => {
                    let all_on = checked.iter().all(|&c| c);
                    for c in checked.iter_mut() { *c = !all_on; }
                }
                KeyCode::Enter => {
                    let selected: Vec<String> = packs.iter().zip(&checked)
                        .filter(|(_, &on)| on)
                        .map(|(p, _)| p.id.to_string())
                        .collect();
                    return Ok(Some(selected));
                }
                KeyCode::Esc | KeyCode::Char('q') => return Ok(None),
                _ => {}
            }
        }
    }
}

fn draw_pack_selector(
    frame: &mut Frame,
    packs: &[&Pack],
    checked: &[bool],
    installed: &[String],
    cursor: usize,
) {
    let area = frame.area();

    let layout = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3),  // title
            Constraint::Length(3),  // always-installed info
            Constraint::Min(8),    // checkbox list
            Constraint::Length(1), // summary
            Constraint::Length(2), // help
        ])
        .split(area);

    // Title
    let title = Paragraph::new("  FAVOR Setup — Add-on Packs")
        .style(Style::default().fg(Color::Cyan).bold())
        .block(Block::default().padding(Padding::top(1)));
    frame.render_widget(title, layout[0]);

    // Always-installed info
    let always: Vec<String> = Pack::required().iter()
        .map(|p| format!("{} ({})", p.name, p.size_human))
        .collect();
    let info = Paragraph::new(Line::from(vec![
        Span::styled("  Always installed: ", Style::default().fg(Color::DarkGray)),
        Span::styled(always.join(", "), Style::default().fg(Color::DarkGray)),
    ]))
    .block(Block::default().padding(Padding::top(1)));
    frame.render_widget(info, layout[1]);

    // Checkbox list
    let items: Vec<ListItem> = packs.iter().enumerate().map(|(i, pack)| {
        let is_cursor = i == cursor;
        let is_installed = installed.iter().any(|id| id == pack.id);
        let mark = if checked[i] { "[x]" } else { "[ ]" };

        let mark_style = if checked[i] {
            Style::default().fg(Color::Green).bold()
        } else if is_cursor {
            Style::default().fg(Color::White)
        } else {
            Style::default().fg(Color::DarkGray)
        };

        let name_style = if is_cursor {
            Style::default().fg(Color::Cyan).bold()
        } else if checked[i] {
            Style::default().fg(Color::White)
        } else {
            Style::default().fg(Color::DarkGray)
        };

        let size_style = if is_cursor {
            Style::default().fg(Color::Yellow)
        } else {
            Style::default().fg(Color::DarkGray)
        };

        let status_tag = if is_installed { " installed " } else { "" };
        let status_style = Style::default().fg(Color::Green);

        let line = Line::from(vec![
            Span::styled(format!("  {mark} "), mark_style),
            Span::styled(format!("{:<16}", pack.id), name_style),
            Span::styled(format!("{:>6}  ", pack.size_human), size_style),
            Span::styled(status_tag, status_style),
            Span::styled(pack.description, name_style),
        ]);
        ListItem::new(line)
    }).collect();

    let list = List::new(items)
        .block(Block::default()
            .borders(Borders::ALL)
            .title(" Select packs (space toggle) ")
            .border_style(Style::default().fg(Color::Cyan)));
    frame.render_widget(list, layout[2]);

    // Summary
    let selected_count = checked.iter().filter(|&&c| c).count();
    let selected_bytes: u64 = packs.iter().zip(checked)
        .filter(|(_, &on)| on)
        .map(|(p, _)| p.size_bytes)
        .sum();
    let total_gb = selected_bytes as f64 / (1024.0 * 1024.0 * 1024.0);
    let summary_text = if selected_count > 0 {
        format!("  Selected: {selected_count} packs ({total_gb:.0} GB)")
    } else {
        "  No packs selected (you can add them later with `favor data pull --pack <name>`)".to_string()
    };
    let summary = Paragraph::new(summary_text)
        .style(Style::default().fg(Color::Yellow));
    frame.render_widget(summary, layout[3]);

    // Help
    let help = Paragraph::new("  up/down navigate    space toggle    a toggle all    enter done    esc skip")
        .style(Style::default().fg(Color::DarkGray));
    frame.render_widget(help, layout[4]);
}

struct EnvOption {
    env: Environment,
    label: &'static str,
    description: &'static str,
}

const ENV_OPTIONS: &[EnvOption] = &[
    EnvOption {
        env: Environment::Hpc,
        label: "HPC cluster",
        description: "Shared cluster with SLURM (srun/sbatch). Memory budget is your \
                       default; srun allocations automatically override it.",
    },
    EnvOption {
        env: Environment::Workstation,
        label: "Workstation",
        description: "Personal machine or dedicated server. Memory budget is the \
                       hard limit for all operations.",
    },
];

/// Interactive environment selector. Returns None if user skipped/cancelled.
fn select_environment() -> io::Result<Option<Environment>> {
    let _guard = TermGuard::enter()?;

    let backend = CrosstermBackend::new(io::stdout());
    let mut terminal = Terminal::new(backend)?;
    let mut selected: usize = 0;

    loop {
        terminal.draw(|frame| draw_environment(frame, selected))?;

        if let Event::Key(key) = event::read()? {
            if key.kind != KeyEventKind::Press { continue; }
            match key.code {
                KeyCode::Up | KeyCode::Char('k') => {
                    if selected > 0 { selected -= 1; }
                }
                KeyCode::Down | KeyCode::Char('j') | KeyCode::Tab => {
                    if selected < ENV_OPTIONS.len() - 1 { selected += 1; }
                }
                KeyCode::Enter => return Ok(Some(ENV_OPTIONS[selected].env)),
                KeyCode::Esc | KeyCode::Char('q') => return Ok(None),
                _ => {}
            }
        }
    }
}

fn draw_environment(frame: &mut Frame, selected: usize) {
    let area = frame.area();

    let layout = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3),  // title
            Constraint::Length(6),  // selector
            Constraint::Min(6),    // description
            Constraint::Length(2), // help
        ])
        .split(area);

    let title = Paragraph::new("  FAVOR Setup — Environment")
        .style(Style::default().fg(Color::Cyan).bold())
        .block(Block::default().padding(Padding::top(1)));
    frame.render_widget(title, layout[0]);

    let mut selector_lines = Vec::new();
    for (i, opt) in ENV_OPTIONS.iter().enumerate() {
        let marker = if i == selected { " > " } else { "   " };
        let style = if i == selected {
            Style::default().fg(Color::Cyan).bold()
        } else {
            Style::default().fg(Color::DarkGray)
        };
        selector_lines.push(Line::from(vec![
            Span::raw(marker),
            Span::styled(opt.label, style),
        ]));
    }
    let selector = Paragraph::new(selector_lines)
        .block(Block::default()
            .borders(Borders::ALL)
            .title(" Are you on an HPC cluster or a workstation? ")
            .border_style(Style::default().fg(Color::Cyan)));
    frame.render_widget(selector, layout[1]);

    // Description panel
    let desc = Paragraph::new(format!("  {}", ENV_OPTIONS[selected].description))
        .wrap(ratatui::widgets::Wrap { trim: true })
        .block(Block::default()
            .borders(Borders::ALL)
            .title(format!(" {} ", ENV_OPTIONS[selected].label))
            .border_style(Style::default().fg(Color::DarkGray)));
    frame.render_widget(desc, layout[2]);

    let help = Paragraph::new("  up/down navigate    enter select    esc skip")
        .style(Style::default().fg(Color::DarkGray));
    frame.render_widget(help, layout[3]);
}

struct MemoryPreset {
    label: &'static str,
    value: &'static str,
    bytes: u64,
}

const MEMORY_PRESETS: &[MemoryPreset] = &[
    MemoryPreset { label: "8 GB",   value: "8GB",   bytes: 8  * 1024 * 1024 * 1024 },
    MemoryPreset { label: "16 GB",  value: "16GB",  bytes: 16 * 1024 * 1024 * 1024 },
    MemoryPreset { label: "32 GB",  value: "32GB",  bytes: 32 * 1024 * 1024 * 1024 },
    MemoryPreset { label: "64 GB",  value: "64GB",  bytes: 64 * 1024 * 1024 * 1024 },
    MemoryPreset { label: "128 GB", value: "128GB", bytes: 128 * 1024 * 1024 * 1024 },
    MemoryPreset { label: "256 GB", value: "256GB", bytes: 256 * 1024 * 1024 * 1024 },
];

/// Interactive memory budget selector. Shows detected memory as context.
/// Returns None if user skipped, or Some("16GB") etc.
fn select_memory_budget(resources: &Resources) -> io::Result<Option<String>> {
    let _guard = TermGuard::enter()?;

    let backend = CrosstermBackend::new(io::stdout());
    let mut terminal = Terminal::new(backend)?;

    // Default selection: nearest preset <= detected memory
    let detected = resources.memory_bytes;
    let mut selected: usize = MEMORY_PRESETS.iter()
        .rposition(|p| p.bytes <= detected * 100 / 80) // undo 80% factor
        .unwrap_or(1); // default to 16GB

    let custom_idx = MEMORY_PRESETS.len(); // "Custom" is after all presets
    let mut custom_mode = false;
    let mut custom_buf = String::new();

    loop {
        let draw_selected = selected;
        let draw_custom = custom_mode;
        let draw_buf = custom_buf.clone();
        terminal.draw(|frame| {
            draw_memory_budget(frame, resources, draw_selected, custom_idx, draw_custom, &draw_buf);
        })?;

        if let Event::Key(key) = event::read()? {
            if key.kind != KeyEventKind::Press { continue; }

            if custom_mode {
                match key.code {
                    KeyCode::Enter => {
                        if ResourceConfig::parse_memory_bytes(&custom_buf).is_some() {
                            return Ok(Some(custom_buf));
                        }
                        // Invalid input — stay in custom mode
                    }
                    KeyCode::Esc => {
                        custom_mode = false;
                        custom_buf.clear();
                    }
                    KeyCode::Backspace => { custom_buf.pop(); }
                    KeyCode::Char(c) => { custom_buf.push(c); }
                    _ => {}
                }
                continue;
            }

            match key.code {
                KeyCode::Up | KeyCode::Char('k') => {
                    if selected > 0 { selected -= 1; }
                }
                KeyCode::Down | KeyCode::Char('j') | KeyCode::Tab => {
                    if selected < custom_idx { selected += 1; }
                }
                KeyCode::Enter => {
                    if selected == custom_idx {
                        custom_mode = true;
                        custom_buf.clear();
                    } else {
                        return Ok(Some(MEMORY_PRESETS[selected].value.to_string()));
                    }
                }
                KeyCode::Esc | KeyCode::Char('q') => return Ok(None),
                _ => {}
            }
        }
    }
}

fn draw_memory_budget(
    frame: &mut Frame,
    resources: &Resources,
    selected: usize,
    custom_idx: usize,
    custom_mode: bool,
    custom_buf: &str,
) {
    let area = frame.area();

    let layout = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3),  // title
            Constraint::Length(3),  // detected info
            Constraint::Min(10),   // preset list
            Constraint::Length(2), // help
        ])
        .split(area);

    let title = Paragraph::new("  FAVOR Setup — Memory Budget")
        .style(Style::default().fg(Color::Cyan).bold())
        .block(Block::default().padding(Padding::top(1)));
    frame.render_widget(title, layout[0]);

    let detected_text = format!(
        "  Detected: {}    Environment: {}",
        resources.memory_human(), resources.environment(),
    );
    let info = Paragraph::new(Line::from(vec![
        Span::styled(&detected_text, Style::default().fg(Color::Yellow)),
    ]))
    .block(Block::default().padding(Padding::top(1)));
    frame.render_widget(info, layout[1]);

    // Build list items
    let mut items: Vec<ListItem> = MEMORY_PRESETS.iter().enumerate().map(|(i, preset)| {
        let is_sel = i == selected && !custom_mode;
        let marker = if is_sel { " > " } else { "   " };
        let style = if is_sel {
            Style::default().fg(Color::Cyan).bold()
        } else {
            Style::default().fg(Color::White)
        };
        let note = if preset.bytes <= resources.memory_bytes * 100 / 80 {
            Span::styled("  (within detected)", Style::default().fg(Color::Green))
        } else {
            Span::styled("  (exceeds detected)", Style::default().fg(Color::Yellow))
        };
        ListItem::new(Line::from(vec![
            Span::raw(marker),
            Span::styled(format!("{:<10}", preset.label), style),
            note,
        ]))
    }).collect();

    // Custom entry
    let is_custom_sel = selected == custom_idx && !custom_mode;
    if custom_mode {
        let valid = ResourceConfig::parse_memory_bytes(custom_buf).is_some();
        let color = if valid { Color::Green } else { Color::Red };
        items.push(ListItem::new(Line::from(vec![
            Span::styled(" > Custom: ", Style::default().fg(Color::Cyan).bold()),
            Span::styled(custom_buf, Style::default().fg(color)),
            Span::styled("_", Style::default().fg(Color::Cyan)),
            if valid {
                Span::styled("  (valid)", Style::default().fg(Color::Green))
            } else {
                Span::styled("  (e.g. 48GB, 12288MB)", Style::default().fg(Color::DarkGray))
            },
        ])));
    } else {
        let style = if is_custom_sel {
            Style::default().fg(Color::Cyan).bold()
        } else {
            Style::default().fg(Color::White)
        };
        let marker = if is_custom_sel { " > " } else { "   " };
        items.push(ListItem::new(Line::from(vec![
            Span::raw(marker),
            Span::styled("Custom...", style),
        ])));
    }

    let list = List::new(items)
        .block(Block::default()
            .borders(Borders::ALL)
            .title(" Default memory budget (srun allocations override this) ")
            .border_style(Style::default().fg(Color::Cyan)));
    frame.render_widget(list, layout[2]);

    let help_text = if custom_mode {
        "  type amount (e.g. 48GB)    enter confirm    esc cancel"
    } else {
        "  up/down navigate    enter select    esc skip"
    };
    let help = Paragraph::new(help_text)
        .style(Style::default().fg(Color::DarkGray));
    frame.render_widget(help, layout[3]);
}
