use crossterm::event::{KeyCode, KeyModifiers};
use ratatui::layout::{Constraint, Direction, Layout, Rect};
use ratatui::style::{Style, Stylize};
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;
use ratatui::Frame;

use crate::tui::action::{format_binding, Action, ActionScope, KeyMap};
use crate::tui::screen::{Screen, Transition};
use crate::tui::theme;
use crate::tui::widgets::log_tail::LogTail;

pub struct HelpScreen {
    active_scope: ActionScope,
}

impl HelpScreen {
    pub fn new(active_scope: ActionScope) -> Self {
        Self { active_scope }
    }
}

impl Screen for HelpScreen {
    fn title(&self) -> &str {
        "Help"
    }

    fn scope(&self) -> ActionScope {
        ActionScope::Help
    }

    fn keys(&self) -> KeyMap {
        let none = KeyModifiers::NONE;
        KeyMap::new()
            .bind(KeyCode::Esc, none, Action::HelpClose)
            .bind(KeyCode::Char('q'), none, Action::HelpClose)
            .bind(KeyCode::Char('?'), none, Action::HelpClose)
    }

    fn on_action(&mut self, action: Action) -> Transition {
        match action {
            Action::HelpClose => Transition::Pop,
            _ => Transition::Stay,
        }
    }

    fn draw(&mut self, frame: &mut Frame, area: Rect, _log: &LogTail) {
        let dim = Style::default().fg(theme::MUTED);
        let key_style = Style::default().fg(theme::ACCENT).bold();
        let scope_style = Style::default().fg(theme::FG).bold();

        let chrome = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(1),
                Constraint::Min(1),
                Constraint::Length(1),
                Constraint::Length(1),
            ])
            .split(area);

        let header = Layout::default()
            .direction(Direction::Horizontal)
            .constraints([Constraint::Min(1), Constraint::Length(20)])
            .split(chrome[0]);
        frame.render_widget(
            Paragraph::new(Span::styled("help", scope_style)),
            header[0],
        );
        frame.render_widget(
            Paragraph::new(Line::from(Span::styled(
                format!("scope {}", self.active_scope.title()),
                dim,
            )))
            .alignment(ratatui::layout::Alignment::Right),
            header[1],
        );

        let mut lines: Vec<Line> = Vec::with_capacity(Action::all().len() + ActionScope::ordered().len());
        let mut first = true;
        for scope in ActionScope::ordered() {
            let mut bound: Vec<&Action> = Action::all()
                .iter()
                .filter(|a| a.scope() == *scope && a.default_key().is_some())
                .collect();
            if bound.is_empty() {
                continue;
            }
            if !first {
                lines.push(Line::from(""));
            }
            first = false;
            let marker = if *scope == self.active_scope { "▌ " } else { "  " };
            lines.push(Line::from(vec![
                Span::styled(marker, key_style),
                Span::styled(scope.title(), scope_style),
            ]));
            bound.sort_by_key(|a| a.title());
            for action in bound {
                let (c, m) = action.default_key().unwrap();
                let key = format_binding(c, m);
                lines.push(Line::from(vec![
                    Span::raw("  "),
                    Span::styled(format!("{key:>12}"), key_style),
                    Span::raw("  "),
                    Span::styled(action.title(), dim),
                ]));
            }
        }
        frame.render_widget(Paragraph::new(lines), chrome[1]);

        frame.render_widget(Paragraph::new(""), chrome[2]);
        frame.render_widget(
            Paragraph::new(Span::styled("esc close   q close   ? close", dim)),
            chrome[3],
        );
    }
}
