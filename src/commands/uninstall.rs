use crate::config::Config;
use crate::error::FavorError;
use crate::output::Output;

pub fn run(out: &dyn Output) -> Result<(), FavorError> {
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
