use std::path::PathBuf;

use crate::config::{Environment, ResourceConfig};

/// Detected system resources for DuckDB and parallel execution.
///
/// Memory resolution priority:
/// 1. SLURM allocation (SLURM_MEM_PER_NODE) — srun jobs get full allocation
/// 2. Config memory_budget — user's chosen default from setup
/// 3. Auto-detect — cgroup → /proc/meminfo → 4GB fallback
#[derive(Debug, Clone)]
pub struct Resources {
    pub memory_bytes: u64,
    pub threads: usize,
    pub temp_dir: PathBuf,
    /// Configured environment, if set during setup
    environment: Option<Environment>,
}

impl Resources {
    /// Auto-detect resources without config. Used during setup (before config exists).
    pub fn detect() -> Self {
        Self {
            memory_bytes: detect_memory(),
            threads: detect_threads(None),
            temp_dir: detect_temp(None),
            environment: None,
        }
    }

    /// Detect resources respecting user's configured budget.
    ///
    /// - If inside a SLURM job, uses the SLURM allocation (overrides budget).
    /// - Otherwise, uses the configured memory_budget as the ceiling.
    /// - Falls back to auto-detect if no budget is configured.
    pub fn detect_with_config(config: &ResourceConfig) -> Self {
        let memory_bytes = if let Some(bytes) = slurm_memory() {
            bytes
        } else if let Some(budget) = config.memory_budget_bytes() {
            budget * 80 / 100
        } else {
            detect_memory()
        };

        Self {
            memory_bytes,
            threads: detect_threads(config.threads),
            temp_dir: detect_temp(Some(config)),
            environment: config.environment,
        }
    }

    pub fn memory_human(&self) -> String {
        let gb = self.memory_bytes as f64 / (1024.0 * 1024.0 * 1024.0);
        format!("{gb:.1} GiB")
    }

    pub fn environment(&self) -> &'static str {
        if std::env::var("SLURM_JOB_ID").is_ok() {
            return "SLURM job";
        }
        if let Some(env) = self.environment {
            return match env {
                Environment::Hpc => "HPC",
                Environment::Workstation => "workstation",
            };
        }
        if let Ok(hostname) = std::env::var("HOSTNAME") {
            if hostname.contains("login") {
                return "login node";
            }
            return "compute";
        }
        "local"
    }

    /// DuckDB memory_limit string.
    ///
    /// Uses memory_bytes directly — the 80% safety margin is already applied
    /// at the budget level (in detect_with_config). DuckDB also spills to
    /// temp_dir when it exceeds this limit, so this is a soft ceiling.
    pub fn duckdb_memory(&self) -> String {
        let gb = self.memory_bytes / (1024 * 1024 * 1024);
        if gb > 0 {
            format!("{}GB", gb)
        } else {
            let mb = self.memory_bytes / (1024 * 1024);
            format!("{}MB", mb.max(512))
        }
    }

}

/// Check for SLURM memory allocation. Returns bytes if inside an srun job.
fn slurm_memory() -> Option<u64> {
    let val = std::env::var("SLURM_MEM_PER_NODE").ok()?;
    let mb = val.trim_end_matches('M').parse::<u64>().ok()?;
    Some(mb * 1024 * 1024 * 80 / 100)
}

fn detect_memory() -> u64 {
    // SLURM job allocation - most reliable on HPC
    if let Some(bytes) = slurm_memory() {
        return bytes;
    }

    // cgroup memory limit (works on login nodes with resource limits)
    for path in &[
        "/sys/fs/cgroup/memory/memory.limit_in_bytes",
        "/sys/fs/cgroup/memory.max",
    ] {
        if let Ok(content) = std::fs::read_to_string(path) {
            let trimmed = content.trim();
            if trimmed != "max" && trimmed != "9223372036854771712" {
                if let Ok(bytes) = trimmed.parse::<u64>() {
                    if bytes < 1024 * 1024 * 1024 * 1024 { // sanity: < 1TB
                        return bytes * 80 / 100;
                    }
                }
            }
        }
    }

    // /proc/meminfo but cap at reasonable limit on shared login nodes
    if let Ok(content) = std::fs::read_to_string("/proc/meminfo") {
        for line in content.lines() {
            if line.starts_with("MemAvailable:") {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if let Some(kb_str) = parts.get(1) {
                    if let Ok(kb) = kb_str.parse::<u64>() {
                        return kb * 1024 * 75 / 100;
                    }
                }
            }
        }
    }

    // Fallback: 4GB
    4 * 1024 * 1024 * 1024
}

fn detect_threads(config_threads: Option<usize>) -> usize {
    // SLURM_CPUS_PER_TASK overrides everything
    if let Ok(val) = std::env::var("SLURM_CPUS_PER_TASK") {
        if let Ok(n) = val.parse::<usize>() {
            return n;
        }
    }

    // Config override
    if let Some(n) = config_threads {
        return n;
    }

    // cgroup cpu quota
    if let Ok(quota) = std::fs::read_to_string("/sys/fs/cgroup/cpu/cpu.cfs_quota_us") {
        if let Ok(period) = std::fs::read_to_string("/sys/fs/cgroup/cpu/cpu.cfs_period_us") {
            if let (Ok(q), Ok(p)) = (
                quota.trim().parse::<i64>(),
                period.trim().parse::<i64>(),
            ) {
                if q > 0 && p > 0 {
                    return (q / p).max(1) as usize;
                }
            }
        }
    }

    // nproc
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4)
}

fn detect_temp(config: Option<&ResourceConfig>) -> PathBuf {
    // Config override is highest priority
    if let Some(cfg) = config {
        if !cfg.temp_dir.is_empty() {
            let p = PathBuf::from(&cfg.temp_dir);
            if p.exists() {
                return p;
            }
        }
    }

    for var in &["TMPDIR", "SCRATCH", "LOCAL_SCRATCH"] {
        if let Ok(val) = std::env::var(var) {
            let p = PathBuf::from(&val);
            if p.exists() {
                return p;
            }
        }
    }

    // /scratch if exists
    let scratch = PathBuf::from("/scratch");
    if scratch.exists() {
        return scratch;
    }

    std::env::temp_dir()
}
