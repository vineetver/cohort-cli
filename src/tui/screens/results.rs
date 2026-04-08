use std::fs::File;
use std::path::{Path, PathBuf};

use arrow::array::{Array, Float64Array, StringArray, UInt32Array};
use crossterm::event::{KeyCode, KeyModifiers};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use ratatui::layout::{Constraint, Direction, Layout, Rect};
use ratatui::style::{Style, Stylize};
use ratatui::symbols::Marker;
use ratatui::text::{Line, Span};
use ratatui::widgets::canvas::{Canvas, Line as CanvasLine, Points};
use ratatui::widgets::{Block, Borders, Cell, Paragraph, Row, Table};
use ratatui::Frame;

use crate::error::CohortError;
use crate::tui::action::{Action, ActionScope, KeyMap};
use crate::tui::screen::{Screen, Transition};
use crate::tui::screens::variant::VariantScreen;
use crate::tui::state::artifacts::{classify, ArtifactKind};
use crate::tui::theme;
use crate::tui::widgets::log_tail::LogTail;
use crate::tui::widgets::status_bar::StatusBar;

const SIGNIFICANCE: f64 = 5e-8;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SortKey {
    PValue,
    Gene,
    NVariants,
}

impl SortKey {
    fn next(self) -> Self {
        match self {
            Self::PValue => Self::Gene,
            Self::Gene => Self::NVariants,
            Self::NVariants => Self::PValue,
        }
    }

    fn label(self) -> &'static str {
        match self {
            Self::PValue => "p_value",
            Self::Gene => "gene",
            Self::NVariants => "n_variants",
        }
    }
}

#[derive(Debug, Clone)]
struct ResultRow {
    gene: String,
    chromosome: String,
    n_variants: u32,
    p_value: f64,
}

pub struct ResultsScreen {
    title: String,
    results_dir: PathBuf,
    mask_files: Vec<PathBuf>,
    mask_idx: usize,
    rows: Vec<ResultRow>,
    sort: SortKey,
    focus: usize,
    error: Option<String>,
    notice: Option<String>,
}

impl ResultsScreen {
    pub fn new(results_dir: PathBuf) -> Result<Self, CohortError> {
        let mask_files = list_mask_parquets(&results_dir)?;
        if mask_files.is_empty() {
            return Err(CohortError::Input(format!(
                "no mask result parquet files in {}",
                results_dir.display()
            )));
        }
        let title = format!(
            "STAAR results: {}",
            results_dir
                .file_name()
                .map(|n| n.to_string_lossy().into_owned())
                .unwrap_or_else(|| results_dir.display().to_string())
        );
        let mut s = Self {
            title,
            results_dir,
            mask_files,
            mask_idx: 0,
            rows: Vec::new(),
            sort: SortKey::PValue,
            focus: 0,
            error: None,
            notice: None,
        };
        s.load_current_mask();
        Ok(s)
    }

    fn current_mask_name(&self) -> String {
        self.mask_files
            .get(self.mask_idx)
            .and_then(|p| p.file_stem())
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_default()
    }

    fn load_current_mask(&mut self) {
        let Some(path) = self.mask_files.get(self.mask_idx).cloned() else {
            self.rows.clear();
            return;
        };
        match read_result_rows(&path) {
            Ok(rows) => {
                self.rows = rows;
                self.error = None;
            }
            Err(e) => {
                self.rows.clear();
                self.error = Some(format!("{e}"));
            }
        }
        self.focus = 0;
        self.apply_sort();
    }

    fn apply_sort(&mut self) {
        match self.sort {
            SortKey::PValue => self.rows.sort_by(|a, b| {
                a.p_value
                    .partial_cmp(&b.p_value)
                    .unwrap_or(std::cmp::Ordering::Equal)
            }),
            SortKey::Gene => self.rows.sort_by(|a, b| a.gene.cmp(&b.gene)),
            SortKey::NVariants => self.rows.sort_by(|a, b| b.n_variants.cmp(&a.n_variants)),
        }
        self.focus = self.focus.min(self.rows.len().saturating_sub(1));
    }

    fn move_focus(&mut self, delta: isize) {
        if self.rows.is_empty() {
            self.focus = 0;
            return;
        }
        let len = self.rows.len() as isize;
        let next = (self.focus as isize + delta).clamp(0, len - 1);
        self.focus = next as usize;
    }

    fn cycle_sort(&mut self) {
        self.sort = self.sort.next();
        self.apply_sort();
    }

    fn cycle_mask(&mut self) {
        if self.mask_files.is_empty() {
            return;
        }
        self.mask_idx = (self.mask_idx + 1) % self.mask_files.len();
        self.load_current_mask();
    }

    fn open_variants_for_focused(&mut self) -> Transition {
        let Some(row) = self.rows.get(self.focus) else {
            return Transition::Stay;
        };
        let Some(annotated) = find_sibling_annotated_set(&self.results_dir) else {
            self.notice = Some(format!(
                "no sibling annotated set found near {}",
                self.results_dir.display()
            ));
            return Transition::Stay;
        };
        let filter = format!("gene_name = \"{}\"", row.gene);
        match VariantScreen::new_for_annotated_set(annotated) {
            Ok(screen) => {
                self.notice = None;
                Transition::Push(Box::new(screen.with_initial_filter(filter)))
            }
            Err(e) => {
                self.notice = Some(format!("cannot open variants: {e}"));
                Transition::Stay
            }
        }
    }

    fn draw_table(&self, frame: &mut Frame, area: Rect) {
        let block = Block::default()
            .borders(Borders::ALL)
            .title(format!(
                " {} · sort:{} · {}/{} mask · {} rows ",
                self.title,
                self.sort.label(),
                self.mask_idx + 1,
                self.mask_files.len(),
                self.rows.len(),
            ))
            .border_style(Style::default().fg(theme::ACCENT));
        let inner = block.inner(area);
        frame.render_widget(block, area);

        if self.rows.is_empty() {
            let msg = self
                .error
                .as_deref()
                .unwrap_or("  no rows in this mask file");
            frame.render_widget(
                Paragraph::new(Line::from(Span::styled(
                    msg,
                    Style::default().fg(theme::MUTED),
                ))),
                inner,
            );
            return;
        }

        let mask = self.current_mask_name();
        let header = Row::new([
            Cell::from("gene"),
            Cell::from("mask"),
            Cell::from("chrom"),
            Cell::from("p_value"),
            Cell::from("n_variants"),
        ])
        .style(Style::default().fg(theme::ACCENT).bold());

        let visible = inner.height.saturating_sub(2) as usize;
        let total = self.rows.len();
        let start = self
            .focus
            .saturating_sub(visible / 2)
            .min(total.saturating_sub(visible.max(1)));
        let end = (start + visible).min(total);

        let rows: Vec<Row> = (start..end)
            .map(|i| {
                let r = &self.rows[i];
                let hot = r.p_value < SIGNIFICANCE && !r.p_value.is_nan();
                let base = if hot {
                    Style::default().fg(theme::BAD).bold()
                } else {
                    Style::default().fg(theme::FG)
                };
                let style = if i == self.focus {
                    base.bg(theme::MUTED)
                } else {
                    base
                };
                Row::new([
                    Cell::from(r.gene.clone()),
                    Cell::from(mask.clone()),
                    Cell::from(r.chromosome.clone()),
                    Cell::from(format_pvalue(r.p_value)),
                    Cell::from(r.n_variants.to_string()),
                ])
                .style(style)
            })
            .collect();

        let widths = [
            Constraint::Length(18),
            Constraint::Length(28),
            Constraint::Length(6),
            Constraint::Length(12),
            Constraint::Length(12),
        ];
        let table = Table::new(rows, widths).header(header);
        frame.render_widget(table, inner);
    }

    fn draw_manhattan(&self, frame: &mut Frame, area: Rect) {
        let block = Block::default()
            .borders(Borders::ALL)
            .title(" Manhattan -log10(p) ")
            .border_style(Style::default().fg(theme::MUTED));
        let inner = block.inner(area);
        frame.render_widget(block, area);

        let mut points: Vec<(f64, f64)> = Vec::with_capacity(self.rows.len());
        let mut max_y: f64 = 1.0;
        for (i, r) in self.rows.iter().enumerate() {
            let y = neg_log10(r.p_value);
            if y > max_y {
                max_y = y;
            }
            points.push((i as f64, y));
        }
        let n = self.rows.len().max(1) as f64;

        let canvas = Canvas::default()
            .marker(Marker::Braille)
            .x_bounds([0.0, n])
            .y_bounds([0.0, max_y * 1.1])
            .paint(move |ctx| {
                ctx.draw(&Points {
                    coords: &points,
                    color: theme::ACCENT,
                });
                ctx.draw(&CanvasLine {
                    x1: 0.0,
                    y1: -SIGNIFICANCE.log10(),
                    x2: n,
                    y2: -SIGNIFICANCE.log10(),
                    color: theme::BAD,
                });
            });
        frame.render_widget(canvas, inner);
    }

    fn draw_qq(&self, frame: &mut Frame, area: Rect) {
        let block = Block::default()
            .borders(Borders::ALL)
            .title(" QQ obs vs exp ")
            .border_style(Style::default().fg(theme::MUTED));
        let inner = block.inner(area);
        frame.render_widget(block, area);

        let mut sorted: Vec<f64> = self
            .rows
            .iter()
            .map(|r| r.p_value)
            .filter(|p| p.is_finite() && *p > 0.0)
            .collect();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let n = sorted.len();
        let mut points: Vec<(f64, f64)> = Vec::with_capacity(n);
        let mut max_axis: f64 = 1.0;
        for (i, p) in sorted.iter().enumerate() {
            let expected = -((i + 1) as f64 / (n as f64 + 1.0)).log10();
            let observed = -p.log10();
            if expected > max_axis {
                max_axis = expected;
            }
            if observed > max_axis {
                max_axis = observed;
            }
            points.push((expected, observed));
        }

        let bound = max_axis * 1.1;
        let canvas = Canvas::default()
            .marker(Marker::Braille)
            .x_bounds([0.0, bound])
            .y_bounds([0.0, bound])
            .paint(move |ctx| {
                ctx.draw(&CanvasLine {
                    x1: 0.0,
                    y1: 0.0,
                    x2: bound,
                    y2: bound,
                    color: theme::MUTED,
                });
                ctx.draw(&Points {
                    coords: &points,
                    color: theme::ACCENT,
                });
            });
        frame.render_widget(canvas, inner);
    }
}

impl Screen for ResultsScreen {
    fn title(&self) -> &str {
        &self.title
    }

    fn error_slot(&self) -> Option<&str> {
        self.notice.as_deref().or(self.error.as_deref())
    }

    fn draw(&mut self, frame: &mut Frame, area: Rect, _log: &LogTail) {
        let v = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Min(8),
                Constraint::Length(8),
                Constraint::Length(8),
                Constraint::Length(1),
            ])
            .split(area);
        self.draw_table(frame, v[0]);
        self.draw_manhattan(frame, v[1]);
        self.draw_qq(frame, v[2]);

        StatusBar {
            title: &self.title,
            keys: "↑↓ navigate · enter variants · s sort · m mask · q back",
        }
        .render(frame, v[3]);
    }

    fn scope(&self) -> ActionScope {
        ActionScope::Results
    }

    fn keys(&self) -> KeyMap {
        let none = KeyModifiers::NONE;
        KeyMap::new()
            .bind(KeyCode::Char('q'), none, Action::ResultsClose)
            .bind(KeyCode::Esc, none, Action::ResultsClose)
            .bind(KeyCode::Char('j'), none, Action::ResultsScrollDown)
            .bind(KeyCode::Down, none, Action::ResultsScrollDown)
            .bind(KeyCode::Char('k'), none, Action::ResultsScrollUp)
            .bind(KeyCode::Up, none, Action::ResultsScrollUp)
            .bind(KeyCode::Char('s'), none, Action::ResultsCycleSort)
            .bind(KeyCode::Char('m'), none, Action::ResultsCycleMask)
            .bind(KeyCode::Enter, none, Action::ResultsOpenVariants)
    }

    fn on_action(&mut self, action: Action) -> Transition {
        match action {
            Action::ResultsClose => Transition::Pop,
            Action::ResultsScrollDown => {
                self.notice = None;
                self.move_focus(1);
                Transition::Stay
            }
            Action::ResultsScrollUp => {
                self.notice = None;
                self.move_focus(-1);
                Transition::Stay
            }
            Action::ResultsCycleSort => {
                self.notice = None;
                self.cycle_sort();
                Transition::Stay
            }
            Action::ResultsCycleMask => {
                self.notice = None;
                self.cycle_mask();
                Transition::Stay
            }
            Action::ResultsOpenVariants => self.open_variants_for_focused(),
            _ => Transition::Stay,
        }
    }
}

fn list_mask_parquets(dir: &Path) -> Result<Vec<PathBuf>, CohortError> {
    let entries = std::fs::read_dir(dir)
        .map_err(|e| CohortError::Resource(format!("read {}: {e}", dir.display())))?;
    let mut out: Vec<PathBuf> = Vec::new();
    for entry in entries {
        let Ok(entry) = entry else { continue };
        let path = entry.path();
        if path.is_file() && path.extension().map(|x| x == "parquet").unwrap_or(false) {
            out.push(path);
        }
    }
    out.sort();
    Ok(out)
}

fn read_result_rows(path: &Path) -> Result<Vec<ResultRow>, CohortError> {
    let file = File::open(path)
        .map_err(|e| CohortError::Resource(format!("open {}: {e}", path.display())))?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)
        .map_err(|e| CohortError::Input(format!("parquet open {}: {e}", path.display())))?;
    let reader = builder
        .build()
        .map_err(|e| CohortError::Analysis(format!("parquet build: {e}")))?;
    let mut rows: Vec<ResultRow> = Vec::new();
    for batch in reader {
        let batch = batch.map_err(|e| CohortError::Analysis(format!("parquet read: {e}")))?;
        let schema = batch.schema();
        let gene_idx = schema
            .index_of("gene_symbol")
            .map_err(|_| CohortError::Input("missing gene_symbol".into()))?;
        let chrom_idx = schema
            .index_of("chromosome")
            .map_err(|_| CohortError::Input("missing chromosome".into()))?;
        let nv_idx = schema
            .index_of("n_variants")
            .map_err(|_| CohortError::Input("missing n_variants".into()))?;
        let p_idx = schema
            .index_of("STAAR-O")
            .map_err(|_| CohortError::Input("missing STAAR-O".into()))?;
        let gene = batch
            .column(gene_idx)
            .as_any()
            .downcast_ref::<StringArray>()
            .ok_or_else(|| CohortError::Input("gene_symbol not Utf8".into()))?;
        let chrom = batch
            .column(chrom_idx)
            .as_any()
            .downcast_ref::<StringArray>()
            .ok_or_else(|| CohortError::Input("chromosome not Utf8".into()))?;
        let nv = batch
            .column(nv_idx)
            .as_any()
            .downcast_ref::<UInt32Array>()
            .ok_or_else(|| CohortError::Input("n_variants not UInt32".into()))?;
        let pv = batch
            .column(p_idx)
            .as_any()
            .downcast_ref::<Float64Array>()
            .ok_or_else(|| CohortError::Input("STAAR-O not Float64".into()))?;
        for i in 0..batch.num_rows() {
            rows.push(ResultRow {
                gene: gene.value(i).to_string(),
                chromosome: chrom.value(i).to_string(),
                n_variants: nv.value(i),
                p_value: if pv.is_null(i) { f64::NAN } else { pv.value(i) },
            });
        }
    }
    Ok(rows)
}

fn find_sibling_annotated_set(start: &Path) -> Option<PathBuf> {
    let mut cursor = start.parent()?;
    for _ in 0..6 {
        if let Ok(entries) = std::fs::read_dir(cursor) {
            for entry in entries.flatten() {
                let p = entry.path();
                if let Some(art) = classify(&p) {
                    if matches!(art.kind, ArtifactKind::AnnotatedSet { .. }) {
                        return Some(p);
                    }
                }
            }
        }
        cursor = cursor.parent()?;
    }
    None
}

fn neg_log10(p: f64) -> f64 {
    if p.is_nan() || p <= 0.0 {
        0.0
    } else {
        -p.log10()
    }
}

fn format_pvalue(p: f64) -> String {
    if p.is_nan() {
        "NaN".to_string()
    } else if p <= 0.0 {
        "0".to_string()
    } else if p < 1e-3 {
        format!("{p:.2e}")
    } else {
        format!("{p:.4}")
    }
}
