use ratatui::buffer::Buffer;
use ratatui::layout::Rect;
use ratatui::style::{Style, Stylize};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Clear, Paragraph, Widget};
use ratatui::Frame;

use crate::tui::action::{format_binding, Action, ActionScope};
use crate::tui::shell::{ScreenChrome, Shell};
use crate::tui::theme;

pub struct HelpState {
    pub scopes: Vec<ActionScope>,
    pub focused: usize,
}

impl HelpState {
    pub fn new(origin: ActionScope) -> Self {
        let origin = if origin == ActionScope::Help {
            ActionScope::Global
        } else {
            origin
        };
        let mut scopes: Vec<ActionScope> = Vec::with_capacity(ActionScope::ordered().len());
        scopes.push(origin);
        for s in ActionScope::ordered() {
            if *s != origin && Action::all().iter().any(|a| a.scope() == *s) {
                scopes.push(*s);
            }
        }
        Self { scopes, focused: 0 }
    }

    pub fn cycle(&mut self) {
        if !self.scopes.is_empty() {
            self.focused = (self.focused + 1) % self.scopes.len();
        }
    }
}

pub fn render(frame: &mut Frame, area: Rect, state: &HelpState) {
    let focused_scope = state
        .scopes
        .get(state.focused)
        .copied()
        .unwrap_or(ActionScope::Global);
    let title = format!("Help — {}", focused_scope.title());
    Clear.render(area, frame.buffer_mut());
    let chrome = ScreenChrome {
        title: &title,
        status: None,
        error: None,
        scope: ActionScope::Help,
        graph: None,
    };
    let scopes = state.scopes.clone();
    let focused = state.focused;
    let body = |inner: Rect, buf: &mut Buffer| {
        let mut lines: Vec<Line> = Vec::new();
        for (i, scope) in scopes.iter().enumerate() {
            let is_focused = i == focused;
            let glyph = if is_focused { "▸ " } else { "  " };
            let header_style = if is_focused {
                Style::default().fg(theme::ACCENT).bold()
            } else {
                Style::default().fg(theme::MUTED).bold()
            };
            lines.push(Line::from(vec![
                Span::styled(glyph, header_style),
                Span::styled(scope.title().to_string(), header_style),
            ]));
            if !is_focused {
                continue;
            }
            for action in Action::all().iter().filter(|a| a.scope() == *scope) {
                let key = action
                    .default_key()
                    .map(|(c, m)| format_binding(c, m))
                    .unwrap_or_else(|| "—".to_string());
                lines.push(Line::from(vec![
                    Span::styled("    ", Style::default()),
                    Span::styled(format!("{key:<14}"), Style::default().fg(theme::WARN)),
                    Span::styled(
                        format!("{:<26}", action.title()),
                        Style::default().fg(theme::FG),
                    ),
                    Span::styled(action.description(), Style::default().fg(theme::MUTED)),
                ]));
            }
        }
        Paragraph::new(lines).render(inner, buf);
    };
    Shell::new(chrome, body).render(area, frame.buffer_mut());
}
