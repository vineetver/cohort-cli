use std::path::PathBuf;
use std::sync::{Arc, Mutex};

use crate::config::Config;
use crate::tui::action::ActionScope;
use crate::tui::screens::help::HelpState;
use crate::tui::screens::stage_view::StageViewState;
use crate::tui::screens::variant::VariantViewState;
use crate::tui::shell::ErrorMessage;
use crate::tui::stages::types::{RunRequest, SetupConfig};
use crate::tui::widgets::palette::Palette;
use crate::tui::widgets::run_overlay::RunOverlay;

use super::session::{SessionId, SessionState, SessionStore};
use super::workspace::WorkspaceState;

pub struct AppState {
    pub workspace: WorkspaceState,
    pub view: View,
    pub modal: Option<Modal>,
    pub error: Option<ErrorMessage>,
    pub session_store: Option<SessionStore>,
    pub setup_sink: Option<Arc<Mutex<Option<SetupConfig>>>>,
    pub list_state: ratatui::widgets::ListState,
    pub pending_focus: Option<PathBuf>,
}

pub enum View {
    Workspace,
    Stage(StageViewState),
    Variant(VariantViewState),
}

pub enum Modal {
    Help(HelpState),
    Palette(Palette),
    Run(RunOverlay),
}

pub enum Outcome {
    Stay,
    Quit,
    Run(RunRequest),
}

impl AppState {
    pub fn new(cwd: PathBuf) -> Self {
        let mut extra_roots: Vec<PathBuf> = Vec::new();
        if let Ok(cfg) = Config::load() {
            let r = cfg.root_dir();
            if !r.as_os_str().is_empty() && r.exists() && r != cwd {
                extra_roots.push(r);
            }
        }
        let workspace = WorkspaceState::new(cwd, extra_roots);
        let mut list_state = ratatui::widgets::ListState::default();
        if !workspace.artifacts.is_empty() {
            list_state.select(Some(0));
        }
        Self {
            workspace,
            view: View::Workspace,
            modal: None,
            error: None,
            session_store: SessionStore::from_home(),
            setup_sink: None,
            list_state,
            pending_focus: None,
        }
    }

    pub fn active_scope(&self) -> ActionScope {
        if let Some(modal) = &self.modal {
            return match modal {
                Modal::Help(_) => ActionScope::Help,
                Modal::Palette(_) => ActionScope::Palette,
                Modal::Run(_) => ActionScope::Run,
            };
        }
        match &self.view {
            View::Workspace => ActionScope::Workspace,
            View::Stage(s) => s.active_scope(),
            View::Variant(v) => v.active_scope(),
        }
    }

    pub fn collect_session(&self) -> SessionState {
        SessionState {
            cwd: self.workspace.cwd.clone(),
            last_artifact: self.workspace.focused().map(|a| a.path.clone()),
        }
    }

    pub fn restore_session(&mut self, state: &SessionState) {
        if let Some(p) = &state.last_artifact {
            self.pending_focus = Some(p.clone());
            self.try_apply_pending_focus();
            self.sync_focus();
        }
    }

    pub fn save_session(&self) {
        let Some(store) = self.session_store.as_ref() else {
            return;
        };
        let state = self.collect_session();
        if state.cwd.as_os_str().is_empty() {
            return;
        }
        let id = SessionId::for_cwd(&state.cwd);
        let _ = store.save(&id, &state);
    }

    pub fn try_apply_pending_focus(&mut self) {
        let Some(target) = self.pending_focus.as_ref() else {
            return;
        };
        if let Some(idx) = self.workspace.artifacts.iter().position(|a| &a.path == target) {
            self.workspace.focus = idx;
            self.pending_focus = None;
        }
    }

    pub fn sync_focus(&mut self) {
        self.list_state.select(if self.workspace.artifacts.is_empty() {
            None
        } else {
            Some(self.workspace.focus)
        });
    }

}
