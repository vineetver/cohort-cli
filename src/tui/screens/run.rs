use crossterm::event::{KeyCode, KeyModifiers};
use ratatui::layout::{Constraint, Direction, Layout, Rect};
use ratatui::style::{Style, Stylize};
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;
use ratatui::Frame;

use crate::tui::action::{Action, ActionScope, KeyMap};
use crate::tui::event::AppEvent;
use crate::tui::screen::{Screen, Transition};
use crate::tui::theme;
use crate::tui::widgets::log_tail::LogTail;
use crate::tui::widgets::status_bar::StatusBar;

pub enum RunState {
    Running,
    Done,
    Err(String),
}

pub struct RunScreen {
    title: String,
    state: RunState,
}

impl RunScreen {
    pub fn new(title: String) -> Self {
        Self {
            title,
            state: RunState::Running,
        }
    }
}

impl Screen for RunScreen {
    fn title(&self) -> &str {
        &self.title
    }

    fn scope(&self) -> ActionScope {
        ActionScope::Run
    }

    fn keys(&self) -> KeyMap {
        let none = KeyModifiers::NONE;
        match self.state {
            RunState::Running => KeyMap::new()
                .bind(KeyCode::Esc, none, Action::RunCancelRequest)
                .bind(KeyCode::Enter, none, Action::RunReturn),
            RunState::Done | RunState::Err(_) => KeyMap::new()
                .bind(KeyCode::Enter, none, Action::RunReturn)
                .bind(KeyCode::Esc, none, Action::RunReturn),
        }
    }

    fn on_action(&mut self, action: Action) -> Transition {
        match action {
            Action::RunCancelRequest | Action::RunReturn => Transition::Pop,
            _ => Transition::Stay,
        }
    }

    fn on_other_event(&mut self, event: &AppEvent) -> Transition {
        match event {
            AppEvent::CommandDone(res) => {
                self.state = match res {
                    Ok(()) => RunState::Done,
                    Err(e) => RunState::Err(e.to_string()),
                };
                Transition::Stay
            }
            _ => Transition::Stay,
        }
    }

    fn draw(&mut self, frame: &mut Frame, area: Rect, log: &LogTail) {
        let v = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(1),
                Constraint::Min(1),
                Constraint::Length(1),
                Constraint::Length(1),
            ])
            .split(area);

        let (status_label, status_color) = match &self.state {
            RunState::Running => ("running", theme::ACCENT),
            RunState::Done => ("done", theme::FG),
            RunState::Err(_) => ("failed", theme::BAD),
        };
        let title_w = self.title.chars().count();
        let status_text = format!("[{status_label}]");
        let status_w = status_text.chars().count();
        let total_w = v[0].width as usize;
        let pad = total_w
            .saturating_sub(title_w + status_w + 2)
            .max(1);
        let header = Paragraph::new(Line::from(vec![
            Span::styled(
                format!(" {}", self.title),
                Style::default().fg(theme::FG).bold(),
            ),
            Span::raw(" ".repeat(pad)),
            Span::styled(status_text, Style::default().fg(status_color).bold()),
            Span::raw(" "),
        ]));
        frame.render_widget(header, v[0]);

        log.draw(frame, v[1], "output");

        let err_line = match &self.state {
            RunState::Err(msg) => {
                Line::from(Span::styled(format!(" {msg}"), theme::error_slot_style()))
            }
            _ => Line::from(""),
        };
        frame.render_widget(Paragraph::new(err_line), v[2]);

        let keys = match &self.state {
            RunState::Running => "detach [esc]   quit [ctrl-c]",
            RunState::Done => "return [enter]",
            RunState::Err(_) => "return [enter]",
        };
        StatusBar {
            title: &self.title,
            keys,
        }
        .render(frame, v[3]);
    }
}
