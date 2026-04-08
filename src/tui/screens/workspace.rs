use std::path::PathBuf;

use crossterm::event::{KeyCode, KeyModifiers};
use ratatui::layout::{Constraint, Direction, Layout, Rect};
use ratatui::style::{Style, Stylize};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Borders, List, ListItem, ListState, Paragraph};
use ratatui::Frame;

use crate::config::Config;
use crate::tui::action::{Action, ActionScope, KeyMap};
use crate::tui::screen::{Screen, Transition};
use crate::tui::screens::setup::SetupScreen;
use crate::tui::screens::transform::TransformScreen;
use crate::tui::screens::variant::VariantScreen;
use crate::tui::state::artifacts::ArtifactKind;
use crate::tui::state::workspace::WorkspaceState;
use crate::tui::theme;
use crate::tui::widgets::log_tail::LogTail;
use crate::tui::widgets::status_bar::StatusBar;

pub struct WorkspaceScreen {
    title: String,
    state: WorkspaceState,
    list_state: ListState,
    notice: Option<String>,
}

impl WorkspaceScreen {
    pub fn new(cwd: PathBuf) -> Self {
        let mut extra_roots: Vec<PathBuf> = Vec::new();
        if let Ok(cfg) = Config::load() {
            let r = cfg.root_dir();
            if !r.as_os_str().is_empty() && r.exists() && r != cwd {
                extra_roots.push(r);
            }
        }
        let state = WorkspaceState::new(cwd, extra_roots);
        let mut list_state = ListState::default();
        if !state.artifacts.is_empty() {
            list_state.select(Some(0));
        }
        Self {
            title: "Workspace".to_string(),
            state,
            list_state,
            notice: None,
        }
    }

    fn sync_focus(&mut self) {
        self.list_state.select(if self.state.artifacts.is_empty() {
            None
        } else {
            Some(self.state.focus)
        });
    }
}

fn actions_for(kind: &ArtifactKind) -> &'static [&'static str] {
    match kind {
        ArtifactKind::RawVcf => &["ingest", "browse"],
        ArtifactKind::PhenotypeTsv => &["use as phenotype"],
        ArtifactKind::KinshipTsv => &["use as kinship"],
        ArtifactKind::ParquetFile => &["browse"],
        ArtifactKind::IngestedSet => &["annotate", "browse"],
        ArtifactKind::AnnotatedSet { .. } => &["staar", "enrich", "browse"],
        ArtifactKind::GenotypeStore => &["inspect"],
        ArtifactKind::StaarResults => &["browse"],
        ArtifactKind::AnnotationRoot => &["status"],
    }
}

fn fmt_size(n: u64) -> String {
    const KB: f64 = 1024.0;
    const MB: f64 = KB * 1024.0;
    const GB: f64 = MB * 1024.0;
    let n = n as f64;
    if n >= GB {
        format!("{:.1} GB", n / GB)
    } else if n >= MB {
        format!("{:.1} MB", n / MB)
    } else if n >= KB {
        format!("{:.1} KB", n / KB)
    } else {
        format!("{n} B")
    }
}

impl Screen for WorkspaceScreen {
    fn title(&self) -> &str {
        &self.title
    }

    fn on_focus(&mut self) {
        self.state.rescan();
        self.sync_focus();
    }

    fn draw(&mut self, frame: &mut Frame, area: Rect, log: &LogTail) {
        if self.state.drain_scan() {
            self.sync_focus();
        }

        let v = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(1),
                Constraint::Min(4),
                Constraint::Length(4),
                Constraint::Length(1),
                Constraint::Length(1),
            ])
            .split(area);

        let header = Layout::default()
            .direction(Direction::Horizontal)
            .constraints([Constraint::Min(0), Constraint::Length(28)])
            .split(v[0]);

        let count = self.state.artifacts.len();
        let status_right = if self.state.scanning {
            format!("{count} artifacts · scanning ")
        } else {
            format!("{count} artifacts ")
        };
        frame.render_widget(
            Paragraph::new(Line::from(Span::styled(
                format!(" {}", self.title),
                Style::default().fg(theme::FG).bold(),
            ))),
            header[0],
        );
        frame.render_widget(
            Paragraph::new(Line::from(Span::styled(
                status_right,
                Style::default().fg(theme::MUTED),
            )))
            .alignment(ratatui::layout::Alignment::Right),
            header[1],
        );

        let h = Layout::default()
            .direction(Direction::Horizontal)
            .constraints([Constraint::Percentage(58), Constraint::Percentage(42)])
            .split(v[1]);

        let focused_idx = self.state.focus;
        let items: Vec<ListItem> = self
            .state
            .artifacts
            .iter()
            .enumerate()
            .map(|(i, a)| {
                let is_focus = i == focused_idx;
                let gutter = if is_focus { "▍ " } else { "  " };
                let name_style = if is_focus {
                    Style::default().fg(theme::FG).bold()
                } else {
                    Style::default().fg(theme::FG)
                };
                let line = Line::from(vec![
                    Span::styled(gutter, Style::default().fg(theme::ACCENT)),
                    Span::styled(format!("{} ", a.kind.glyph()), Style::default().fg(theme::MUTED)),
                    Span::styled(format!("{:<28}", a.display_name), name_style),
                    Span::styled(
                        a.path.to_string_lossy().into_owned(),
                        Style::default().fg(theme::MUTED),
                    ),
                ]);
                ListItem::new(line)
            })
            .collect();

        let list = List::new(items)
            .block(
                Block::default()
                    .borders(Borders::ALL)
                    .border_style(Style::default().fg(theme::ACCENT)),
            )
            .highlight_style(Style::default().fg(theme::FG));
        frame.render_stateful_widget(list, h[0], &mut self.list_state);

        let detail_lines: Vec<Line> = match self.state.focused() {
            Some(a) => {
                let mut lines = vec![Line::from(Span::styled(
                    a.kind.title().to_string(),
                    Style::default().fg(theme::FG).bold(),
                ))];
                lines.push(Line::from(Span::styled(
                    a.path.to_string_lossy().into_owned(),
                    Style::default().fg(theme::MUTED),
                )));
                if a.size_bytes > 0 {
                    lines.push(Line::from(Span::styled(
                        fmt_size(a.size_bytes),
                        Style::default().fg(theme::MUTED),
                    )));
                }
                if let ArtifactKind::AnnotatedSet { tier } = &a.kind {
                    lines.push(Line::from(Span::styled(
                        tier.as_str().to_string(),
                        Style::default().fg(theme::MUTED),
                    )));
                }
                lines.push(Line::from(Span::styled(
                    actions_for(&a.kind).join("  "),
                    Style::default().fg(theme::FG),
                )));
                lines
            }
            None => {
                let msg = if self.state.scanning {
                    "scanning…"
                } else {
                    "no artifacts — r to rescan"
                };
                vec![Line::from(Span::styled(
                    msg,
                    Style::default().fg(theme::MUTED),
                ))]
            }
        };

        frame.render_widget(Paragraph::new(detail_lines), h[1]);

        log.draw(frame, v[2], "");

        let err_line = match &self.notice {
            Some(msg) => Line::from(Span::styled(
                format!(" {msg}"),
                Style::default().fg(theme::BAD),
            )),
            None => Line::from(""),
        };
        frame.render_widget(Paragraph::new(err_line), v[3]);

        StatusBar {
            title: &self.title,
            keys: "j/k move  enter open  r rescan  s setup  q quit",
        }
        .render(frame, v[4]);
    }

    fn scope(&self) -> ActionScope {
        ActionScope::Workspace
    }

    fn keys(&self) -> KeyMap {
        let none = KeyModifiers::NONE;
        KeyMap::new()
            .bind(KeyCode::Char('q'), none, Action::Quit)
            .bind(KeyCode::Esc, none, Action::Quit)
            .bind(KeyCode::Down, none, Action::WorkspaceDown)
            .bind(KeyCode::Char('j'), none, Action::WorkspaceDown)
            .bind(KeyCode::Up, none, Action::WorkspaceUp)
            .bind(KeyCode::Char('k'), none, Action::WorkspaceUp)
            .bind(KeyCode::Char('r'), none, Action::WorkspaceRescan)
            .bind(KeyCode::Enter, none, Action::WorkspaceOpenFocused)
            .bind(KeyCode::Char('s'), none, Action::WorkspaceOpenSetup)
    }

    fn on_action(&mut self, action: Action) -> Transition {
        match action {
            Action::Quit => Transition::Quit,
            Action::WorkspaceDown => {
                self.state.move_focus(1);
                self.sync_focus();
                Transition::Stay
            }
            Action::WorkspaceUp => {
                self.state.move_focus(-1);
                self.sync_focus();
                Transition::Stay
            }
            Action::WorkspaceRescan => {
                self.state.rescan();
                self.sync_focus();
                self.notice = None;
                Transition::Stay
            }
            Action::WorkspaceOpenSetup => Transition::Push(Box::new(SetupScreen::new())),
            Action::WorkspaceOpenFocused => match self.state.focused() {
                Some(a) => match &a.kind {
                    ArtifactKind::RawVcf => {
                        self.notice = None;
                        Transition::Push(Box::new(TransformScreen::new_ingest(Some(a))))
                    }
                    ArtifactKind::IngestedSet => {
                        self.notice = None;
                        Transition::Push(Box::new(TransformScreen::new_annotate(a)))
                    }
                    ArtifactKind::AnnotatedSet { .. } => {
                        match VariantScreen::new_for_annotated_set(a.path.clone()) {
                            Ok(screen) => {
                                self.notice = None;
                                Transition::Push(Box::new(screen))
                            }
                            Err(e) => {
                                self.notice = Some(format!("cannot open: {e}"));
                                Transition::Stay
                            }
                        }
                    }
                    ArtifactKind::ParquetFile => {
                        match VariantScreen::new_for_parquet(a.path.clone()) {
                            Ok(screen) => {
                                self.notice = None;
                                Transition::Push(Box::new(screen))
                            }
                            Err(e) => {
                                self.notice = Some(format!("cannot open: {e}"));
                                Transition::Stay
                            }
                        }
                    }
                    other => {
                        self.notice = Some(format!("no transform for {}", other.title()));
                        Transition::Stay
                    }
                },
                None => Transition::Stay,
            },
            _ => Transition::Stay,
        }
    }
}
