use std::path::{Path, PathBuf};

use crossterm::event::{KeyCode, KeyModifiers};
use ratatui::layout::{Constraint, Direction, Layout, Rect};
use ratatui::style::{Style, Stylize};
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;
use ratatui::Frame;

use crate::cli::GenomeBuild;
use crate::commands::{AnnotateConfig, IngestConfig};
use crate::config::{Config, Tier};
use crate::tui::action::{Action, ActionScope, KeyMap};
use crate::tui::screen::{RunRequest, Screen, Transition};
use crate::tui::state::artifacts::{Artifact, ArtifactKind};
use crate::tui::theme;
use crate::tui::widgets::file_picker::{self, DirBrowserState};
use crate::tui::widgets::log_tail::LogTail;

pub struct IngestForm {
    inputs: Vec<PathBuf>,
    output: Option<PathBuf>,
    output_user_set: bool,
    emit_sql: bool,
    build: Option<GenomeBuild>,
}

impl IngestForm {
    fn from_artifact(art: Option<&Artifact>) -> Self {
        let inputs = match art {
            Some(a) if matches!(a.kind, ArtifactKind::RawVcf) => vec![a.path.clone()],
            _ => Vec::new(),
        };
        let output = inputs.first().map(|p| default_ingest_output(p));
        Self {
            inputs,
            output,
            output_user_set: false,
            emit_sql: false,
            build: None,
        }
    }

    fn refresh_default_output(&mut self) {
        if self.output_user_set {
            return;
        }
        self.output = self.inputs.first().map(|p| default_ingest_output(p));
    }

    fn to_config(&self) -> Result<IngestConfig, String> {
        if self.inputs.is_empty() {
            return Err("at least one input file is required".into());
        }
        for p in &self.inputs {
            if !p.exists() {
                return Err(format!("input not found: {}", p.display()));
            }
        }
        let output = self
            .output
            .clone()
            .unwrap_or_else(|| default_ingest_output(&self.inputs[0]));
        Ok(IngestConfig {
            inputs: self.inputs.clone(),
            output,
            emit_sql: self.emit_sql,
            build_override: self.build.clone(),
        })
    }

    fn command_preview(&self) -> String {
        let mut parts = vec!["cohort ingest".to_string()];
        for p in &self.inputs {
            parts.push(short_name(p));
        }
        if let Some(o) = &self.output {
            parts.push(format!("-o {}", short_name(o)));
        }
        if self.emit_sql {
            parts.push("--emit-sql".into());
        }
        match self.build {
            Some(GenomeBuild::Hg38) => parts.push("--build hg38".into()),
            Some(GenomeBuild::Hg19) => parts.push("--build hg19".into()),
            None => {}
        }
        parts.join(" ")
    }
}

fn default_ingest_output(first: &Path) -> PathBuf {
    let stem = first.file_stem().unwrap_or_default().to_string_lossy().into_owned();
    let stem = stem
        .strip_suffix(".vcf")
        .or_else(|| stem.strip_suffix(".tsv"))
        .or_else(|| stem.strip_suffix(".csv"))
        .unwrap_or(&stem);
    let stem = stem
        .split("_b0_")
        .next()
        .or_else(|| stem.split("_b0.").next())
        .unwrap_or(stem);
    first
        .parent()
        .unwrap_or(first)
        .join(format!("{stem}.ingested"))
}

pub struct AnnotateForm {
    input: PathBuf,
    output: Option<PathBuf>,
    output_user_set: bool,
    tier: Tier,
    data_root: PathBuf,
}

impl AnnotateForm {
    fn from_artifact(art: &Artifact) -> Self {
        let cfg = Config::load().ok();
        let tier = cfg.as_ref().map(|c| c.data.tier).unwrap_or(Tier::Base);
        let data_root = cfg.map(|c| c.root_dir()).unwrap_or_else(PathBuf::new);
        Self {
            input: art.path.clone(),
            output: Some(default_annotate_output(&art.path)),
            output_user_set: false,
            tier,
            data_root,
        }
    }

    fn refresh_default_output(&mut self) {
        if self.output_user_set {
            return;
        }
        self.output = Some(default_annotate_output(&self.input));
    }

    fn to_config(&self) -> Result<AnnotateConfig, String> {
        if !self.input.exists() {
            return Err(format!("input not found: {}", self.input.display()));
        }
        if !self.input.join("meta.json").exists() {
            return Err(format!(
                "'{}' is not an ingested set (missing meta.json)",
                self.input.display()
            ));
        }
        if self.data_root.as_os_str().is_empty() {
            return Err("data root not configured — run setup from workspace".into());
        }
        if !self.data_root.exists() {
            return Err(format!("data root not found: {}", self.data_root.display()));
        }
        let output = self
            .output
            .clone()
            .unwrap_or_else(|| default_annotate_output(&self.input));
        Ok(AnnotateConfig {
            input: self.input.clone(),
            output,
            tier: self.tier,
            data_root: self.data_root.clone(),
        })
    }

    fn command_preview(&self) -> String {
        let mut parts = vec!["cohort annotate".to_string(), short_name(&self.input)];
        if matches!(self.tier, Tier::Full) {
            parts.push("--full".into());
        }
        if let Some(o) = &self.output {
            parts.push(format!("-o {}", short_name(o)));
        }
        parts.join(" ")
    }
}

fn default_annotate_output(input: &Path) -> PathBuf {
    let name = input.file_name().unwrap_or_default().to_string_lossy().into_owned();
    let stem = name.strip_suffix(".ingested").unwrap_or(&name);
    input
        .parent()
        .unwrap_or(input)
        .join(format!("{stem}.annotated"))
}

fn short_name(p: &Path) -> String {
    p.file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_else(|| p.display().to_string())
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum FieldKind {
    Path,
    Toggle,
    Cycle,
    Expander,
    Run,
}

struct FieldSpec {
    kind: FieldKind,
    label: &'static str,
    advanced: bool,
}

const ING_FIELDS: &[FieldSpec] = &[
    FieldSpec { kind: FieldKind::Path, label: "input", advanced: false },
    FieldSpec { kind: FieldKind::Path, label: "output", advanced: false },
    FieldSpec { kind: FieldKind::Expander, label: "advanced", advanced: false },
    FieldSpec { kind: FieldKind::Toggle, label: "emit SQL", advanced: true },
    FieldSpec { kind: FieldKind::Cycle, label: "build", advanced: true },
    FieldSpec { kind: FieldKind::Run, label: "Run", advanced: false },
];

const ANN_FIELDS: &[FieldSpec] = &[
    FieldSpec { kind: FieldKind::Path, label: "input", advanced: false },
    FieldSpec { kind: FieldKind::Path, label: "output", advanced: false },
    FieldSpec { kind: FieldKind::Cycle, label: "tier", advanced: false },
    FieldSpec { kind: FieldKind::Expander, label: "advanced", advanced: false },
    FieldSpec { kind: FieldKind::Path, label: "data root", advanced: true },
    FieldSpec { kind: FieldKind::Run, label: "Run", advanced: false },
];

enum FormState {
    Ingest(IngestForm),
    Annotate(AnnotateForm),
}

impl FormState {
    fn fields(&self) -> &'static [FieldSpec] {
        match self {
            FormState::Ingest(_) => ING_FIELDS,
            FormState::Annotate(_) => ANN_FIELDS,
        }
    }
}

enum PickerTarget {
    IngestAddInput,
    IngestOutput,
    AnnotateOutput,
    AnnotateDataRoot,
}

pub struct TransformScreen {
    title: String,
    form: FormState,
    visible: Vec<usize>,
    focus: usize,
    advanced: bool,
    picker: Option<(DirBrowserState, PickerTarget)>,
    error: Option<String>,
}

impl TransformScreen {
    pub fn new_ingest(focused: Option<&Artifact>) -> Self {
        let form = FormState::Ingest(IngestForm::from_artifact(focused));
        let mut s = Self {
            title: "ingest".into(),
            form,
            visible: Vec::new(),
            focus: 0,
            advanced: false,
            picker: None,
            error: None,
        };
        s.rebuild_visible();
        s
    }

    pub fn new_annotate(art: &Artifact) -> Self {
        let form = FormState::Annotate(AnnotateForm::from_artifact(art));
        let mut s = Self {
            title: "annotate".into(),
            form,
            visible: Vec::new(),
            focus: 0,
            advanced: false,
            picker: None,
            error: None,
        };
        s.rebuild_visible();
        s
    }

    fn rebuild_visible(&mut self) {
        let prev_real = self.visible.get(self.focus).copied();
        self.visible = self
            .form
            .fields()
            .iter()
            .enumerate()
            .filter(|(_, f)| self.advanced || !f.advanced)
            .map(|(i, _)| i)
            .collect();
        self.focus = match prev_real.and_then(|r| self.visible.iter().position(|&i| i == r)) {
            Some(p) => p,
            None => 0,
        };
    }

    fn current_field(&self) -> &'static FieldSpec {
        let real = self.visible[self.focus];
        &self.form.fields()[real]
    }

    fn cycle_focus(&mut self, delta: isize) {
        let n = self.visible.len() as isize;
        if n == 0 {
            return;
        }
        self.focus = ((self.focus as isize + delta).rem_euclid(n)) as usize;
    }

    fn field_value(&self, real_idx: usize) -> String {
        let spec = &self.form.fields()[real_idx];
        match (&self.form, spec.label) {
            (FormState::Ingest(f), "input") => match f.inputs.len() {
                0 => "(empty — Enter to add)".into(),
                1 => f.inputs[0].display().to_string(),
                n => format!("{n} files"),
            },
            (FormState::Ingest(f), "output") => match &f.output {
                Some(p) => p.display().to_string(),
                None => "(none)".into(),
            },
            (FormState::Ingest(f), "emit SQL") => if f.emit_sql { "yes" } else { "no" }.into(),
            (FormState::Ingest(f), "build") => match f.build {
                Some(GenomeBuild::Hg38) => "hg38".into(),
                Some(GenomeBuild::Hg19) => "hg19".into(),
                None => "auto".into(),
            },
            (FormState::Annotate(f), "input") => f.input.display().to_string(),
            (FormState::Annotate(f), "output") => match &f.output {
                Some(p) => p.display().to_string(),
                None => "(none)".into(),
            },
            (FormState::Annotate(f), "tier") => f.tier.as_str().into(),
            (FormState::Annotate(f), "data root") => {
                if f.data_root.as_os_str().is_empty() {
                    "(unset — run setup)".into()
                } else {
                    f.data_root.display().to_string()
                }
            }
            (_, "advanced") => {
                if self.advanced {
                    "[-] advanced (Enter to hide)".into()
                } else {
                    "[+] advanced (Enter to show)".into()
                }
            }
            _ => String::new(),
        }
    }

    fn try_run(&mut self) -> Transition {
        let result: Result<RunRequest, String> = match &self.form {
            FormState::Ingest(f) => f.to_config().map(RunRequest::Ingest),
            FormState::Annotate(f) => f.to_config().map(RunRequest::Annotate),
        };
        match result {
            Ok(req) => Transition::Run(req),
            Err(msg) => {
                self.error = Some(msg);
                Transition::Stay
            }
        }
    }

    fn open_picker_for_focus(&mut self) {
        let real = self.visible[self.focus];
        let spec = &self.form.fields()[real];
        if !matches!(spec.kind, FieldKind::Path) {
            return;
        }
        let (target, start, prompt, show_files) = match (&self.form, spec.label) {
            (FormState::Ingest(_), "input") => (
                PickerTarget::IngestAddInput,
                std::env::current_dir().unwrap_or_else(|_| PathBuf::from(".")),
                "select input file",
                true,
            ),
            (FormState::Ingest(f), "output") => (
                PickerTarget::IngestOutput,
                f.output
                    .as_deref()
                    .and_then(|p| p.parent())
                    .map(PathBuf::from)
                    .or_else(|| std::env::current_dir().ok())
                    .unwrap_or_else(|| PathBuf::from(".")),
                "select output directory",
                false,
            ),
            (FormState::Annotate(f), "output") => (
                PickerTarget::AnnotateOutput,
                f.output
                    .as_deref()
                    .and_then(|p| p.parent())
                    .map(PathBuf::from)
                    .or_else(|| std::env::current_dir().ok())
                    .unwrap_or_else(|| PathBuf::from(".")),
                "select output directory",
                false,
            ),
            (FormState::Annotate(f), "data root") => (
                PickerTarget::AnnotateDataRoot,
                if f.data_root.as_os_str().is_empty() || !f.data_root.exists() {
                    std::env::current_dir().unwrap_or_else(|_| PathBuf::from("/"))
                } else {
                    f.data_root.clone()
                },
                "select data root",
                false,
            ),
            _ => return,
        };
        self.picker = Some((DirBrowserState::with_files(prompt, &start, show_files), target));
    }

    fn apply_picker_choice(&mut self, target: PickerTarget, chosen: PathBuf) {
        match (&mut self.form, target) {
            (FormState::Ingest(f), PickerTarget::IngestAddInput) => {
                f.inputs.push(chosen);
                f.refresh_default_output();
            }
            (FormState::Ingest(f), PickerTarget::IngestOutput) => {
                f.output = Some(chosen);
                f.output_user_set = true;
            }
            (FormState::Annotate(f), PickerTarget::AnnotateOutput) => {
                f.output = Some(chosen);
                f.output_user_set = true;
            }
            (FormState::Annotate(f), PickerTarget::AnnotateDataRoot) => {
                f.data_root = chosen;
            }
            _ => {}
        }
        self.error = None;
    }

    fn activate_field(&mut self) -> Transition {
        let spec = self.current_field();
        match spec.kind {
            FieldKind::Run => self.try_run(),
            FieldKind::Expander => {
                self.advanced = !self.advanced;
                self.rebuild_visible();
                Transition::Stay
            }
            FieldKind::Toggle => {
                if let FormState::Ingest(f) = &mut self.form {
                    if spec.label == "emit SQL" {
                        f.emit_sql = !f.emit_sql;
                    }
                }
                Transition::Stay
            }
            FieldKind::Cycle => {
                match (&mut self.form, spec.label) {
                    (FormState::Ingest(f), "build") => {
                        f.build = match f.build {
                            None => Some(GenomeBuild::Hg38),
                            Some(GenomeBuild::Hg38) => Some(GenomeBuild::Hg19),
                            Some(GenomeBuild::Hg19) => None,
                        };
                    }
                    (FormState::Annotate(f), "tier") => {
                        f.tier = match f.tier {
                            Tier::Base => Tier::Full,
                            Tier::Full => Tier::Base,
                        };
                    }
                    _ => {}
                }
                Transition::Stay
            }
            FieldKind::Path => {
                self.open_picker_for_focus();
                Transition::Stay
            }
        }
    }

    fn clear_input_at_focus(&mut self) {
        let spec = self.current_field();
        match (&mut self.form, spec.label) {
            (FormState::Ingest(f), "input") => {
                f.inputs.clear();
                f.refresh_default_output();
            }
            (FormState::Ingest(f), "output") => {
                f.output_user_set = false;
                f.refresh_default_output();
            }
            (FormState::Annotate(f), "output") => {
                f.output_user_set = false;
                f.refresh_default_output();
            }
            _ => {}
        }
    }

    fn command_preview(&self) -> String {
        match &self.form {
            FormState::Ingest(f) => f.command_preview(),
            FormState::Annotate(f) => f.command_preview(),
        }
    }
}

impl Screen for TransformScreen {
    fn title(&self) -> &str {
        &self.title
    }

    fn draw(&mut self, frame: &mut Frame, area: Rect, log: &LogTail) {
        if let Some((picker, _)) = self.picker.as_mut() {
            file_picker::draw(frame, area, picker);
            return;
        }

        let v = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(1),
                Constraint::Min(4),
                Constraint::Length(1),
                Constraint::Length(1),
                Constraint::Length(1),
                Constraint::Length(if area.height >= 18 { 5 } else { 0 }),
                Constraint::Length(1),
                Constraint::Length(1),
            ])
            .split(area);

        let header_split = Layout::default()
            .direction(Direction::Horizontal)
            .constraints([Constraint::Min(10), Constraint::Length(24)])
            .split(v[0]);
        let title_left = Paragraph::new(Line::from(Span::styled(
            format!(" {} ", self.title),
            Style::default().fg(theme::FG).bold(),
        )));
        let mode = if self.advanced { "advanced" } else { "basic" };
        let status_right = Paragraph::new(Line::from(Span::styled(
            format!("{mode} · {}/{} ", self.visible.len(), self.form.fields().len()),
            Style::default().fg(theme::MUTED),
        )))
        .alignment(ratatui::layout::Alignment::Right);
        frame.render_widget(title_left, header_split[0]);
        frame.render_widget(status_right, header_split[1]);

        let mut rows: Vec<Line> = Vec::with_capacity(self.visible.len());
        for (vi, &real) in self.visible.iter().enumerate() {
            let spec = &self.form.fields()[real];
            let is_focus = vi == self.focus;
            let value = self.field_value(real);
            let glyph = if is_focus { "▌ " } else { "  " };
            let label_style = if is_focus {
                Style::default().fg(theme::ACCENT).bold()
            } else {
                Style::default().fg(theme::MUTED)
            };
            let value_style = match spec.kind {
                FieldKind::Run => Style::default().fg(theme::OK).bold(),
                _ if is_focus => Style::default().fg(theme::FG).bold(),
                _ => Style::default().fg(theme::FG),
            };
            let show_label = !matches!(spec.kind, FieldKind::Run | FieldKind::Expander)
                && !((spec.label == "input" || spec.label == "output") && !is_focus);
            let label_text = if show_label {
                format!("{:<10}", spec.label)
            } else {
                String::new()
            };
            let display_value = if matches!(spec.kind, FieldKind::Run) {
                "Run [Enter]".to_string()
            } else {
                value
            };
            rows.push(Line::from(vec![
                Span::styled(glyph, Style::default().fg(theme::ACCENT)),
                Span::styled(label_text, label_style),
                Span::styled(display_value, value_style),
            ]));
        }
        frame.render_widget(Paragraph::new(rows), v[1]);

        let preview = Paragraph::new(Line::from(Span::styled(
            format!("  $ {}", self.command_preview()),
            Style::default().fg(theme::MUTED),
        )));
        frame.render_widget(preview, v[2]);

        let error_line = match &self.error {
            Some(msg) => Line::from(Span::styled(
                format!("  {msg}"),
                Style::default().fg(theme::BAD).bold(),
            )),
            None => Line::from(""),
        };
        frame.render_widget(Paragraph::new(error_line), v[3]);

        let primary_hint = Line::from(vec![
            Span::raw("  "),
            Span::styled("Run [Enter]", Style::default().fg(theme::OK).bold()),
            Span::styled("   tab fields · space toggle · esc back", Style::default().fg(theme::MUTED)),
        ]);
        frame.render_widget(Paragraph::new(primary_hint), v[4]);

        if v[5].height > 0 {
            log.draw(frame, v[5], "log");
        }

        let hint = Paragraph::new(Line::from(Span::styled(
            " backspace clears · enter on [advanced] toggles flags ",
            Style::default().fg(theme::MUTED),
        )));
        frame.render_widget(hint, v[7]);
    }

    fn scope(&self) -> ActionScope {
        if self.picker.is_some() {
            ActionScope::FilePicker
        } else {
            ActionScope::Transform
        }
    }

    fn keys(&self) -> KeyMap {
        let none = KeyModifiers::NONE;
        if self.picker.is_some() {
            KeyMap::new()
                .bind(KeyCode::Esc, none, Action::PickerCancel)
                .bind(KeyCode::Up, none, Action::PickerUp)
                .bind(KeyCode::Char('k'), none, Action::PickerUp)
                .bind(KeyCode::Down, none, Action::PickerDown)
                .bind(KeyCode::Char('j'), none, Action::PickerDown)
                .bind(KeyCode::Left, none, Action::PickerParent)
                .bind(KeyCode::Backspace, none, Action::PickerParent)
                .bind(KeyCode::Right, none, Action::PickerInto)
                .bind(KeyCode::Enter, none, Action::PickerInto)
        } else {
            KeyMap::new()
                .bind(KeyCode::Esc, none, Action::TransformCancel)
                .bind(KeyCode::Tab, none, Action::TransformNextField)
                .bind(KeyCode::Down, none, Action::TransformNextField)
                .bind(KeyCode::BackTab, KeyModifiers::SHIFT, Action::TransformPrevField)
                .bind(KeyCode::Up, none, Action::TransformPrevField)
                .bind(KeyCode::Enter, none, Action::TransformActivate)
                .bind(KeyCode::Char(' '), none, Action::TransformToggleBool)
                .bind(KeyCode::Backspace, none, Action::TransformClearField)
        }
    }

    fn on_action(&mut self, action: Action) -> Transition {
        if let Some((picker, _)) = self.picker.as_mut() {
            match action {
                Action::PickerCancel => {
                    self.picker = None;
                }
                Action::PickerUp => picker.select_up(),
                Action::PickerDown => picker.select_down(),
                Action::PickerParent => picker.go_parent(),
                Action::PickerInto => {
                    if let Some(chosen) = picker.enter_selected() {
                        let (_, target) = self.picker.take().unwrap();
                        self.apply_picker_choice(target, chosen);
                    }
                }
                _ => {}
            }
            return Transition::Stay;
        }

        match action {
            Action::TransformCancel => Transition::Pop,
            Action::TransformNextField => {
                self.cycle_focus(1);
                Transition::Stay
            }
            Action::TransformPrevField => {
                self.cycle_focus(-1);
                Transition::Stay
            }
            Action::TransformActivate | Action::TransformToggleBool => self.activate_field(),
            Action::TransformClearField => {
                self.clear_input_at_focus();
                Transition::Stay
            }
            _ => Transition::Stay,
        }
    }
}
