use clap::Parser;

mod annotation_db;
mod cli;
mod commands;
mod config;
mod db;
mod error;
mod feedback;
mod ingest;
mod mode;
mod output;
mod packs;
mod resource;
mod staar;
#[cfg(test)]
mod test_fixtures;
mod tui;
mod types;
mod variant_set;

use cli::{Cli, Command};
use error::FavorError;
use mode::OutputMode;

fn main() {
    let cli = Cli::parse();
    let mode = OutputMode::detect(&cli.format);
    let out = output::create(&mode);

    let dry_run = cli.dry_run;
    let result = run(cli.command, &*out, &mode, dry_run);

    if let Err(e) = result {
        out.error(&e);
        std::process::exit(e.exit_code());
    }
}

fn run(
    command: Command,
    out: &dyn output::Output,
    mode: &OutputMode,
    dry_run: bool,
) -> Result<(), FavorError> {
    match command {
        Command::Init { path, force } => commands::init::run(path, force, out, mode),
        Command::Setup { environment, memory_budget } => {
            commands::setup::run(out, mode, environment, memory_budget)
        }
        Command::Data { action } => commands::data::run(action, out),
        Command::Ingest { input, output, emit_sql, build } => {
            commands::ingest::run(input, output, emit_sql, build, out, dry_run)
        }
        Command::Annotate {
            input,
            output: output_path,
            full,
        } => commands::annotate::run(input, output_path, full, out, dry_run),
        Command::Enrich {
            input,
            tissue,
            output: output_path,
        } => commands::enrich::run(input, tissue, output_path, out, dry_run),
        Command::Interpret {
            input,
            tissue,
            disease,
            output: output_path,
        } => commands::interpret::run(input, tissue, disease, output_path, out),
        Command::Staar {
            genotypes,
            phenotype,
            trait_name,
            covariates,
            annotations,
            masks,
            maf_cutoff,
            window_size,
            individual,
            spa,
            ancestry_col,
            scang_lmin,
            scang_lmax,
            scang_step,
            known_loci,
            emit_sumstats,
            adaptive: _adaptive,
            output: output_path,
        } => commands::staar::run(
            genotypes, phenotype, trait_name, covariates, annotations,
            masks, maf_cutoff, window_size, individual, spa, ancestry_col,
            scang_lmin, scang_lmax, scang_step,
            known_loci, emit_sumstats, output_path, out, dry_run,
        ),
        Command::MetaStaar {
            studies,
            masks,
            maf_cutoff,
            window_size,
            output: output_path,
        } => commands::meta_staar::run(
            studies, masks, maf_cutoff, window_size, output_path, out, dry_run,
        ),
        Command::Schema { table } => commands::schema_cmd::run(table, out),
        Command::Manifest => commands::manifest::run(out),
        Command::Uninstall => commands::uninstall::run(out),
    }
}
