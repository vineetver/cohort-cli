use std::path::PathBuf;

use crossterm::event::KeyCode;
use ratatui::buffer::Buffer;
use ratatui::layout::{Constraint, Direction, Layout, Rect};
use ratatui::style::{Style, Stylize};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Borders, Paragraph, Widget};
use ratatui::Frame;

use crate::tui::action::{Action, ActionScope};
use crate::tui::shell::{ErrorMessage, ScreenChrome, Shell};
use crate::tui::stages::types::{FormError, FormField, PathKind, SessionCtx};
use crate::tui::stages::Stage;
use crate::tui::state::app::{AppState, Outcome, View};
use crate::tui::theme::{self, Tone, FOCUS_GLYPH};
use crate::tui::widgets::file_picker::{self, DirBrowserState};
use crate::tui::widgets::form::{Form, FormOutcome};

pub struct MultiPickerState {
    pub field_id: &'static str,
    pub label: &'static str,
    pub options: &'static [&'static str],
    pub cursor: usize,
    pub checked: Vec<bool>,
}

pub enum Editor {
    None,
    Path(DirBrowserState, &'static str),
    Multi(MultiPickerState),
}

pub struct StageViewState {
    pub stage: &'static dyn Stage,
    pub title: String,
    pub form: Form,
    pub editor: Editor,
}

impl StageViewState {
    pub fn new(stage: &'static dyn Stage, ctx: SessionCtx<'_>) -> Self {
        let schema = stage.form_schema(&ctx);
        let form = Form::new(schema, "Run");
        let title = format!("Stage: {}", stage.label());
        Self {
            stage,
            title,
            form,
            editor: Editor::None,
        }
    }

    pub fn active_scope(&self) -> ActionScope {
        match &self.editor {
            Editor::Path(_, _) => ActionScope::FilePicker,
            _ => ActionScope::Transform,
        }
    }
}

fn open_path_picker(view: &mut StageViewState, id: &'static str) {
    let cwd = std::env::current_dir().unwrap_or_else(|_| PathBuf::from("."));
    let kind = match view.form.field(id) {
        Some(FormField::Path { kind, .. }) => *kind,
        _ => return,
    };
    let start = view
        .form
        .values()
        .path(id)
        .and_then(|p| {
            if p.is_dir() {
                Some(p.clone())
            } else {
                p.parent().map(PathBuf::from)
            }
        })
        .unwrap_or(cwd);
    let show_files = !matches!(kind, PathKind::Dir);
    view.editor = Editor::Path(
        DirBrowserState::with_files("Select path", &start, show_files),
        id,
    );
}

fn open_multi_picker(view: &mut StageViewState, id: &'static str) {
    let (label, options) = match view.form.field(id) {
        Some(FormField::MultiSelect { label, options, .. }) => (*label, *options),
        _ => return,
    };
    let current = view.form.values().multi(id).cloned().unwrap_or_default();
    let checked: Vec<bool> = options.iter().map(|o| current.iter().any(|c| c == o)).collect();
    view.editor = Editor::Multi(MultiPickerState {
        field_id: id,
        label,
        options,
        cursor: 0,
        checked,
    });
}

fn commit_multi(view: &mut StageViewState) {
    if let Editor::Multi(state) = &view.editor {
        let id = state.field_id;
        let picks: Vec<String> = state
            .options
            .iter()
            .zip(state.checked.iter())
            .filter(|(_, &on)| on)
            .map(|(o, _)| (*o).to_string())
            .collect();
        view.form.set_multi(id, picks);
    }
    view.editor = Editor::None;
}

fn try_run(state: &mut AppState) -> Outcome {
    let View::Stage(view) = &mut state.view else {
        return Outcome::Stay;
    };
    match view.stage.build_command(view.form.values()) {
        Ok(req) => Outcome::Run(req),
        Err(err) => {
            state.error = Some(ErrorMessage {
                text: form_error_text(&err),
            });
            Outcome::Stay
        }
    }
}

fn drive_form(state: &mut AppState, code: KeyCode) -> Outcome {
    state.error = None;
    let View::Stage(view) = &mut state.view else {
        return Outcome::Stay;
    };
    match view.form.handle(code) {
        FormOutcome::Continue | FormOutcome::OpenAdvanced => Outcome::Stay,
        FormOutcome::Cancel => {
            state.view = View::Workspace;
            Outcome::Stay
        }
        FormOutcome::Submit => try_run(state),
        FormOutcome::RequestEdit(id) => {
            match view.form.field(id) {
                Some(FormField::Path { .. }) => open_path_picker(view, id),
                Some(FormField::MultiSelect { .. }) => open_multi_picker(view, id),
                _ => {}
            }
            Outcome::Stay
        }
    }
}

fn clear_focused(view: &mut StageViewState) {
    let Some(id) = view.form.focused_field_id() else {
        return;
    };
    if matches!(view.form.field(id), Some(FormField::Path { .. })) {
        view.form.set_path(id, None);
    }
}

fn form_error_text(err: &FormError) -> String {
    match err {
        FormError::Missing(field) => format!("missing field: {field}"),
        FormError::Invalid { field, reason } => format!("{field}: {reason}"),
    }
}

pub fn render(state: &mut AppState, frame: &mut Frame, area: Rect) {
    let View::Stage(view) = &mut state.view else {
        return;
    };
    match &mut view.editor {
        Editor::Path(picker, _) => {
            file_picker::draw(frame, area, picker);
            return;
        }
        Editor::Multi(state) => {
            draw_multi(frame, area, state);
            return;
        }
        Editor::None => {}
    }
    let chrome = ScreenChrome {
        title: &view.title,
        status: None,
        error: state.error.as_ref(),
        scope: ActionScope::Transform,
        graph: None,
    };
    let form = &view.form;
    let body = |inner: Rect, buf: &mut Buffer| {
        form.render(inner, buf);
    };
    Shell::new(chrome, body).render(area, frame.buffer_mut());
}

fn draw_multi(frame: &mut Frame, area: Rect, state: &MultiPickerState) {
    let layout = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(2),
            Constraint::Min(4),
            Constraint::Length(1),
        ])
        .split(area);
    let title = Paragraph::new(Line::from(Span::styled(
        format!("  {}", state.label),
        Style::default().fg(theme::ACCENT).bold(),
    )));
    frame.render_widget(title, layout[0]);

    let lines: Vec<Line> = state
        .options
        .iter()
        .enumerate()
        .map(|(i, opt)| {
            let g = if i == state.cursor { FOCUS_GLYPH } else { " " };
            let mark = if state.checked[i] { "[x]" } else { "[ ]" };
            let tone = if i == state.cursor { Tone::Focus } else { Tone::Normal };
            Line::from(vec![
                Span::styled(format!(" {g} {mark} "), tone.style()),
                Span::styled((*opt).to_string(), tone.style()),
            ])
        })
        .collect();
    let body = Paragraph::new(lines).block(
        Block::default()
            .borders(Borders::ALL)
            .border_style(Style::default().fg(theme::ACCENT)),
    );
    frame.render_widget(body, layout[1]);

    let hint = "  space toggle   enter done   esc cancel";
    frame.render_widget(Paragraph::new(hint).style(theme::hint_bar_style()), layout[2]);
}

pub fn handle_action(state: &mut AppState, action: Action) -> Outcome {
    let View::Stage(view) = &mut state.view else {
        return Outcome::Stay;
    };
    if let Editor::Path(picker, _) = &mut view.editor {
        match action {
            Action::PickerCancel => {
                view.editor = Editor::None;
            }
            Action::PickerUp => picker.select_up(),
            Action::PickerDown => picker.select_down(),
            Action::PickerParent => picker.go_parent(),
            Action::PickerInto => {
                if let Some(chosen) = picker.enter_selected() {
                    if let Editor::Path(_, id) =
                        std::mem::replace(&mut view.editor, Editor::None)
                    {
                        view.form.set_path(id, Some(chosen));
                        state.error = None;
                    }
                }
            }
            _ => {}
        }
        return Outcome::Stay;
    }
    if let Editor::Multi(picker) = &mut view.editor {
        match action {
            Action::TransformCancel => view.editor = Editor::None,
            Action::TransformPrevField => {
                if picker.cursor > 0 {
                    picker.cursor -= 1;
                }
            }
            Action::TransformNextField => {
                if picker.cursor + 1 < picker.options.len() {
                    picker.cursor += 1;
                }
            }
            Action::TransformToggleBool => {
                if let Some(c) = picker.checked.get_mut(picker.cursor) {
                    *c = !*c;
                }
            }
            Action::TransformActivate => commit_multi(view),
            _ => {}
        }
        return Outcome::Stay;
    }

    match action {
        Action::TransformCancel => drive_form(state, KeyCode::Esc),
        Action::TransformNextField => drive_form(state, KeyCode::Down),
        Action::TransformPrevField => drive_form(state, KeyCode::Up),
        Action::TransformActivate => drive_form(state, KeyCode::Enter),
        Action::TransformToggleBool => drive_form(state, KeyCode::Char(' ')),
        Action::TransformClearField => {
            clear_focused(view);
            Outcome::Stay
        }
        _ => Outcome::Stay,
    }
}
