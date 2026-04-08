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
use crate::tui::screens::results::ResultsScreen;
use crate::tui::screens::setup::SetupScreen;
use crate::tui::screens::transform::TransformScreen;
use crate::tui::screens::variant::VariantScreen;
use crate::tui::state::artifacts::ArtifactKind;
use crate::tui::state::workspace::{derived_from_label, WorkspaceState};
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

    fn error_slot(&self) -> Option<&str> {
        self.notice.as_deref()
    }

    fn on_focus(&mut self) {
        self.state.rescan();
        self.sync_focus();
    }

    fn focused_path(&self) -> Option<PathBuf> {
        self.state.focused().map(|a| a.path.clone())
    }

    fn draw(&mut self, frame: &mut Frame, area: Rect, log: &LogTail) {
        if self.state.drain_scan() {
            self.sync_focus();
        }

        let v = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Min(4),
                Constraint::Length(6),
                Constraint::Length(1),
            ])
            .split(area);

        let h = Layout::default()
            .direction(Direction::Horizontal)
            .constraints([Constraint::Percentage(55), Constraint::Percentage(45)])
            .split(v[0]);

        let items: Vec<ListItem> = self
            .state
            .artifacts
            .iter()
            .enumerate()
            .map(|(i, a)| {
                let selected = self.state.selection.contains(&i);
                let mark = if selected { "[x] " } else { "[ ] " };
                let mark_style = Style::default().fg(if selected {
                    theme::ACCENT
                } else {
                    theme::MUTED
                });
                let mut spans = vec![
                    Span::styled(mark, mark_style),
                    Span::styled(
                        format!(" {} ", a.kind.glyph()),
                        Style::default().fg(theme::WARN),
                    ),
                    Span::styled(
                        format!("{:<28}", a.display_name),
                        Style::default().fg(theme::FG),
                    ),
                ];
                if let Some(d) = &a.derived_from {
                    spans.push(Span::styled(
                        derived_from_label(d),
                        Style::default().fg(theme::MUTED),
                    ));
                    if self.state.artifact_is_stale(i) {
                        spans.push(Span::styled(" (stale)", Style::default().fg(theme::BAD)));
                    }
                    spans.push(Span::raw("  "));
                }
                spans.push(Span::styled(
                    a.path.to_string_lossy().into_owned(),
                    Style::default().fg(theme::MUTED),
                ));
                ListItem::new(Line::from(spans))
            })
            .collect();

        let visual_tag = if self.state.visual_anchor.is_some() {
            " [VISUAL]"
        } else {
            ""
        };
        let selected_tag = if self.state.selection.is_empty() {
            String::new()
        } else {
            format!(" · {} selected", self.state.selection.len())
        };
        let list_title = if self.state.scanning {
            format!(
                " Artifacts ({}, scanning…){}{} ",
                self.state.artifacts.len(),
                selected_tag,
                visual_tag,
            )
        } else {
            format!(
                " Artifacts ({}){}{} ",
                self.state.artifacts.len(),
                selected_tag,
                visual_tag,
            )
        };
        let list = List::new(items)
            .block(
                Block::default()
                    .borders(Borders::ALL)
                    .title(list_title)
                    .border_style(Style::default().fg(theme::ACCENT)),
            )
            .highlight_style(Style::default().bg(theme::MUTED).fg(theme::FG))
            .highlight_symbol(" > ");
        frame.render_stateful_widget(list, h[0], &mut self.list_state);

        let detail_lines: Vec<Line> = match self.state.focused() {
            Some(a) => {
                let mut lines = vec![
                    Line::from(Span::styled(
                        format!("  {}", a.kind.title()),
                        Style::default().fg(theme::ACCENT).bold(),
                    )),
                    Line::from(""),
                    Line::from(vec![
                        Span::styled("  path  ", Style::default().fg(theme::MUTED)),
                        Span::styled(
                            a.path.to_string_lossy().into_owned(),
                            Style::default().fg(theme::FG),
                        ),
                    ]),
                ];
                if a.size_bytes > 0 {
                    lines.push(Line::from(vec![
                        Span::styled("  size  ", Style::default().fg(theme::MUTED)),
                        Span::styled(fmt_size(a.size_bytes), Style::default().fg(theme::FG)),
                    ]));
                }
                if let ArtifactKind::AnnotatedSet { tier } = &a.kind {
                    lines.push(Line::from(vec![
                        Span::styled("  tier  ", Style::default().fg(theme::MUTED)),
                        Span::styled(tier.as_str(), Style::default().fg(theme::FG)),
                    ]));
                }
                lines.push(Line::from(""));
                lines.push(Line::from(Span::styled(
                    "  actions",
                    Style::default().fg(theme::ACCENT),
                )));
                for action in actions_for(&a.kind) {
                    lines.push(Line::from(Span::styled(
                        format!("    {action}"),
                        Style::default().fg(theme::FG),
                    )));
                }
                lines
            }
            None => {
                let msg = if self.state.scanning {
                    "  Scanning current directory…"
                } else {
                    "  No artifacts found. Press r to rescan."
                };
                vec![Line::from(Span::styled(
                    msg,
                    Style::default().fg(theme::MUTED),
                ))]
            }
        };

        let detail = Paragraph::new(detail_lines).block(
            Block::default()
                .borders(Borders::ALL)
                .title(" Detail ")
                .border_style(Style::default().fg(theme::MUTED)),
        );
        frame.render_widget(detail, h[1]);

        log.draw(frame, v[1], "Log");

        StatusBar {
            title: &self.title,
            keys: "q quit  j/k move  enter transform  p parent  space select  v visual  a group  s staar  S setup  r rescan  ? help",
        }
        .render(frame, v[2]);
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
            .bind(KeyCode::Char('p'), none, Action::WorkspaceOpenParent)
            .bind(KeyCode::Char(' '), none, Action::WorkspaceToggleSelect)
            .bind(KeyCode::Char('v'), none, Action::WorkspaceToggleVisualMode)
            .bind(KeyCode::Char('a'), none, Action::WorkspaceSelectAllInGroup)
            .bind(KeyCode::Char('s'), none, Action::OpenStaarForm)
            .bind(KeyCode::Char('S'), KeyModifiers::SHIFT, Action::WorkspaceOpenSetup)
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
            Action::WorkspaceToggleSelect => {
                self.state.toggle_selection();
                Transition::Stay
            }
            Action::WorkspaceToggleVisualMode => {
                self.state.toggle_visual_mode();
                Transition::Stay
            }
            Action::WorkspaceSelectAllInGroup => {
                self.state.select_all_in_group();
                Transition::Stay
            }
            Action::WorkspaceOpenParent => {
                if let Some(art) = self.state.focused() {
                    if let Some(idx) = self.state.first_parent_index(art) {
                        self.state.set_focus(idx);
                        self.sync_focus();
                    }
                }
                Transition::Stay
            }
            Action::WorkspaceOpenSetup => Transition::Push(Box::new(SetupScreen::new())),
            Action::OpenStaarForm => match self.state.focused() {
                Some(a) if matches!(a.kind, ArtifactKind::AnnotatedSet { .. }) => {
                    self.notice = None;
                    Transition::Push(Box::new(TransformScreen::new_staar(Some(a.path.clone()))))
                }
                Some(_) => {
                    self.notice = None;
                    Transition::Push(Box::new(TransformScreen::new_staar(None)))
                }
                None => Transition::Stay,
            },
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
                    ArtifactKind::StaarResults => {
                        match ResultsScreen::new(a.path.clone()) {
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
