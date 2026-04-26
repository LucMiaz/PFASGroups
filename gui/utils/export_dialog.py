"""
Export Format Dialog
====================
A single reusable dialog opened by "Export…" buttons in the Results,
Prioritisation, and Modelling tabs.

The *context* argument drives which panels are shown:
  "results"        — Flat / Wide / Long layout options + full metric/group selection
  "prioritisation" — exports the ranked table as-is (Rank, Name, SMILES, Score)
  "modelling"      — exports Scores and/or Bayesian t-test tables
"""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QDialog, QDialogButtonBox, QFileDialog,
    QFormLayout, QGridLayout, QGroupBox, QHBoxLayout, QLabel, QMessageBox,
    QPushButton, QRadioButton, QScrollArea, QSizePolicy,
    QVBoxLayout, QWidget, QButtonGroup, QFrame,
)

from gui import style

if TYPE_CHECKING:
    from PFASGroups.PFASEmbeddings import PFASEmbeddingSet

# ── Metric definitions ─────────────────────────────────────────────────────

_COMPONENT_METRICS: list[tuple[str, list[tuple[str, str]]]] = [
    ("Presence", [
        ("binary",          "Binary (0/1 — matched or not)"),
        ("count",           "Count (number of matched components)"),
        ("max_component",   "Max component size per group"),
        ("total_component", "Total component size per group"),
    ]),
    ("Graph topology", [
        ("size",                          "Size (C-atom count)"),
        ("branching",                     "Branching (1 = linear, 0 = fully branched)"),
        ("effective_graph_resistance",     "Effective graph resistance (EGR)"),
        ("effective_graph_resistance_BDE", "EGR – BDE weighted"),
        ("component_fraction",             "Component fraction (% of mol C-atoms)"),
    ]),
    ("Eccentricity", [
        ("mean_eccentricity",   "Mean eccentricity"),
        ("median_eccentricity", "Median eccentricity"),
        ("diameter",            "Diameter (longest shortest path)"),
        ("radius",              "Radius (minimum eccentricity)"),
    ]),
    ("Distance", [
        ("min_dist_to_centre",              "Min distance to centre"),
        ("max_dist_to_periphery",           "Max distance to periphery"),
        ("min_dist_to_barycentre",          "Min distance to barycentre"),
        ("min_resistance_dist_to_centre",   "Min resistance distance to centre"),
        ("max_resistance_dist_to_periphery","Max resistance distance to periphery"),
        ("min_resistance_dist_to_barycentre","Min resistance distance to barycentre"),
    ]),
    ("Structure", [
        ("n_spacer",  "N spacer (telomeric CH₂ linker length)"),
        ("ring_size", "Ring size (smallest overlapping ring, 0 = acyclic)"),
    ]),
]

# Flat set of keys for the graph-topology / eccentricity / distance / structure
# sections — used to decide whether to show the Aggregation widget.
_GRAPH_METRICS: frozenset[str] = frozenset(
    key
    for section, items in _COMPONENT_METRICS[1:]  # skip "Presence"
    for key, _ in items
)

_MOLECULE_METRICS: list[tuple[str, str]] = [
    ("n_components",           "N components (total matched)"),
    ("total_size",             "Total size (sum of component sizes)"),
    ("mean_size",              "Mean component size"),
    ("max_size",               "Max component size"),
    ("mean_branching",         "Mean branching"),
    ("max_branching",          "Max branching"),
    ("mean_eccentricity",      "Mean eccentricity (across components)"),
    ("max_diameter",           "Max diameter"),
    ("mean_component_fraction","Mean component fraction"),
    ("max_component_fraction", "Max component fraction"),
]

_GROUP_SELECTIONS: list[tuple[str, str]] = [
    ("all",              "All groups"),
    ("oecd",             "OECD groups (IDs 1–28)"),
    ("generic",          "Generic groups (IDs 29–76)"),
    ("telomers",         "Telomers (IDs 77–118)"),
    ("generic+telomers", "Generic + Telomers (IDs 29–118)"),
]

_FORMAT_EXTS = {
    "CSV": ".csv",
    "Excel (.xlsx)": ".xlsx",
    "TSV": ".tsv",
    "JSON": ".json",
}


# ── Helper: horizontal rule ────────────────────────────────────────────────

def _hline() -> QFrame:
    line = QFrame()
    line.setFrameShape(QFrame.Shape.HLine)
    line.setFrameShadow(QFrame.Shadow.Sunken)
    line.setStyleSheet(f"color: {style.C_BORDER};")
    return line


# ── Data building ──────────────────────────────────────────────────────────

def _get_groups(emb) -> list[str]:
    """Return list of matched group names for a single embedding."""
    groups: list[str] = []
    try:
        matched = emb.get("matches") or emb.get("groups") or {}
        if isinstance(matched, dict):
            for gid, match_data in matched.items():
                if match_data:
                    name = match_data.get("name") or str(gid)
                    groups.append(name)
        elif isinstance(matched, list):
            for m in matched:
                if isinstance(m, dict):
                    groups.append(m.get("name", str(m.get("id", ""))))
                else:
                    groups.append(str(m))
    except Exception:
        pass
    return groups


def _get_definitions(emb) -> list[str]:
    """Return list of matched PFAS definition IDs (as strings)."""
    defs: list[str] = []
    try:
        for d in (emb.get("pfas_definitions") or emb.get("definitions") or []):
            did = d.get("id") if isinstance(d, dict) else getattr(d, "id", None)
            if did is not None:
                defs.append(str(did))
    except Exception:
        pass
    return defs


def _build_flat(emb_set, cols: set[str]) -> pd.DataFrame:
    rows = []
    for emb in emb_set:
        row: dict = {}
        if "name" in cols:
            row["name"] = emb.get("name", "")
        if "smiles" in cols:
            row["smiles"] = emb.get("smiles", "")
        if "matched_groups" in cols:
            row["matched_groups"] = "; ".join(_get_groups(emb))
        if "n_groups" in cols:
            row["n_groups"] = len(_get_groups(emb))
        if "matched_definitions" in cols:
            row["matched_definitions"] = "; ".join(_get_definitions(emb))
        rows.append(row)
    return pd.DataFrame(rows)


def _build_wide(
    emb_set,
    group_selection: str,
    component_metrics: list[str],
    molecule_metrics: list[str],
    aggregation: str,
) -> pd.DataFrame:
    """Build a wide DataFrame using to_array() + column_names()."""
    kwargs: dict = dict(
        group_selection=group_selection,
        aggregation=aggregation,
        progress=False,
    )
    if component_metrics:
        kwargs["component_metrics"] = component_metrics
    if molecule_metrics:
        kwargs["molecule_metrics"] = molecule_metrics

    col_names: list[str] = emb_set.column_names(
        component_metrics=component_metrics if component_metrics else None,
        molecule_metrics=molecule_metrics if molecule_metrics else None,
        group_selection=group_selection,
    )
    arr = emb_set.to_array(**kwargs)
    # arr may be an EmbeddingArray (subclass of ndarray) — convert to plain array
    import numpy as np
    mat = np.asarray(arr)

    # Prepend name + smiles columns
    names = [emb.get("name", "") for emb in emb_set]
    smiles = [emb.get("smiles", "") for emb in emb_set]
    df = pd.DataFrame(mat, columns=col_names)
    df.insert(0, "smiles", smiles)
    df.insert(0, "name", names)
    return df


def _build_long(emb_set, group_selection: str, cols: set[str]) -> pd.DataFrame:
    """One row per (molecule, matched group)."""
    # Determine allowed group IDs for filtering
    group_id_filter: set[int] | None = None
    if group_selection == "oecd":
        group_id_filter = set(range(1, 29))
    elif group_selection == "generic":
        group_id_filter = set(range(29, 77))
    elif group_selection == "telomers":
        group_id_filter = set(range(77, 119))
    elif group_selection == "generic+telomers":
        group_id_filter = set(range(29, 119))

    rows = []
    for emb in emb_set:
        name = emb.get("name", "")
        smiles = emb.get("smiles", "")
        matched = emb.get("matches") or emb.get("groups") or {}

        if isinstance(matched, dict):
            items = list(matched.items())
        elif isinstance(matched, list):
            items = [(i, m) for i, m in enumerate(matched)]
        else:
            items = []

        for raw_id, match_data in items:
            if not match_data:
                continue
            if isinstance(match_data, dict):
                group_name = match_data.get("name", str(raw_id))
                group_id = match_data.get("group_id") or match_data.get("id") or raw_id
                match_count = match_data.get("count") or match_data.get("match_count") or 1
            elif hasattr(match_data, "group_name"):
                group_name = match_data.group_name or str(raw_id)
                group_id = getattr(match_data, "group_id", raw_id)
                match_count = getattr(match_data, "match_count", 1)
            else:
                group_name = str(raw_id)
                group_id = raw_id
                match_count = 1

            try:
                gid_int = int(group_id)
            except (TypeError, ValueError):
                gid_int = None

            if group_id_filter is not None and gid_int not in group_id_filter:
                continue

            row: dict = {}
            if "name" in cols:
                row["name"] = name
            if "smiles" in cols:
                row["smiles"] = smiles
            if "group_name" in cols:
                row["group_name"] = group_name
            if "group_id" in cols:
                row["group_id"] = group_id
            if "match_count" in cols:
                row["match_count"] = match_count
            rows.append(row)
    return pd.DataFrame(rows)


def _build_prioritisation(ranked: list) -> pd.DataFrame:
    rows = [
        {"rank": i + 1, "name": name, "smiles": smiles, "score": score}
        for i, (name, smiles, score) in enumerate(ranked)
    ]
    return pd.DataFrame(rows)


def _write_df(df: pd.DataFrame, path: str, fmt: str) -> None:
    if fmt == "CSV":
        df.to_csv(path, index=False)
    elif fmt == "Excel (.xlsx)":
        df.to_excel(path, index=False)
    elif fmt == "TSV":
        df.to_csv(path, index=False, sep="\t")
    elif fmt == "JSON":
        df.to_json(path, orient="records", indent=2)


# ── Main dialog ────────────────────────────────────────────────────────────

class ExportDialog(QDialog):
    """
    Pop-up dialog for configuring and performing data exports.

    Parameters
    ----------
    data : PFASEmbeddingSet | list[tuple] | dict
        The data to export.  Interpretation depends on *context*.
    context : str
        One of ``"results"``, ``"prioritisation"``, or ``"modelling"``.
    """

    def __init__(self, data, context: str, parent=None):
        super().__init__(parent)
        self._data = data
        self._context = context
        self.setWindowTitle("Export…")
        self.setMinimumWidth(560)
        self.setMinimumHeight(480)
        self._build_ui()

    # ── UI construction ────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setSpacing(6)
        root.setContentsMargins(14, 14, 14, 10)

        # ── Scroll area wraps format + content ─────────────────────────────
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        scroll_inner = QWidget()
        inner_lay = QVBoxLayout(scroll_inner)
        inner_lay.setSpacing(10)
        inner_lay.setContentsMargins(0, 0, 6, 0)
        scroll.setWidget(scroll_inner)
        root.addWidget(scroll, stretch=1)

        # ── Format ────────────────────────────────────────────────────────
        fmt_grp = QGroupBox("Export format")
        fmt_lay = QHBoxLayout(fmt_grp)
        fmt_lay.setSpacing(16)
        self._fmt_group = QButtonGroup(self)
        for i, label in enumerate(_FORMAT_EXTS):
            rb = QRadioButton(label)
            rb.setChecked(i == 0)
            self._fmt_group.addButton(rb, i)
            fmt_lay.addWidget(rb)
        fmt_lay.addStretch()
        inner_lay.addWidget(fmt_grp)

        # ── Context-specific content ───────────────────────────────────────
        if self._context == "results":
            self._build_results_section(inner_lay)
        elif self._context == "prioritisation":
            note = QLabel("Exports: Rank, Name, SMILES, Score — all rows.")
            note.setObjectName("muted")
            inner_lay.addWidget(note)
        elif self._context == "modelling":
            self._build_modelling_section(inner_lay)

        inner_lay.addStretch()

        # ── Buttons (outside the scroll area, always visible) ──────────────
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Cancel
        )
        export_btn = buttons.addButton(
            "Export…", QDialogButtonBox.ButtonRole.AcceptRole
        )
        export_btn.setObjectName("btn_primary")
        export_btn.clicked.connect(self._do_export)
        buttons.rejected.connect(self.reject)
        root.addWidget(buttons)

    def _build_results_section(self, root: QVBoxLayout) -> None:
        # ── Layout radios ──────────────────────────────────────────────────
        layout_grp = QGroupBox("Layout")
        layout_lay = QHBoxLayout(layout_grp)
        layout_lay.setSpacing(20)
        self._layout_group = QButtonGroup(self)
        for i, label in enumerate(("Flat", "Wide", "Long")):
            rb = QRadioButton(label)
            rb.setChecked(i == 0)
            self._layout_group.addButton(rb, i)
            layout_lay.addWidget(rb)
        layout_lay.addStretch()
        root.addWidget(layout_grp)

        # Container that switches between panels
        self._panels_container = QWidget()
        panels_lay = QVBoxLayout(self._panels_container)
        panels_lay.setContentsMargins(0, 0, 0, 0)
        panels_lay.setSpacing(6)
        root.addWidget(self._panels_container)

        self._flat_panel = self._build_flat_panel()
        self._wide_panel = self._build_wide_panel()
        self._long_panel = self._build_long_panel()

        panels_lay.addWidget(self._flat_panel)
        panels_lay.addWidget(self._wide_panel)
        panels_lay.addWidget(self._long_panel)

        self._wide_panel.setVisible(False)
        self._long_panel.setVisible(False)

        self._layout_group.idToggled.connect(self._on_layout_changed)

    def _build_flat_panel(self) -> QGroupBox:
        grp = QGroupBox("Columns to include")
        lay = QVBoxLayout(grp)
        self._flat_cbs: dict[str, QCheckBox] = {}
        for key, label in [
            ("name",               "Name / identifier"),
            ("smiles",             "SMILES"),
            ("matched_groups",     "Matched groups (semicolon-separated list)"),
            ("n_groups",           "Number of matched groups"),
            ("matched_definitions","Matched PFAS definitions"),
        ]:
            cb = QCheckBox(label)
            cb.setChecked(True)
            self._flat_cbs[key] = cb
            lay.addWidget(cb)
        return grp

    def _build_wide_panel(self) -> QWidget:
        outer = QWidget()
        lay = QVBoxLayout(outer)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(8)

        # Group filter
        grp_filter = QGroupBox("Group filter")
        gf_lay = QFormLayout(grp_filter)
        self._wide_group_sel = QComboBox()
        for key, label in _GROUP_SELECTIONS:
            self._wide_group_sel.addItem(label, key)
        gf_lay.addRow("Groups:", self._wide_group_sel)
        lay.addWidget(grp_filter)

        # ── Component metrics in a scrollable two-column grid ──────────────
        comp_grp = QGroupBox("Component metrics")
        comp_outer_lay = QVBoxLayout(comp_grp)
        comp_outer_lay.setSpacing(4)

        # Select All / Clear All toolbar
        btn_row = QHBoxLayout()
        sel_all = QPushButton("Select All")
        sel_all.setObjectName("btn_secondary")
        sel_all.setFixedWidth(90)
        clr_all = QPushButton("Clear All")
        clr_all.setObjectName("btn_secondary")
        clr_all.setFixedWidth(90)
        sel_all.clicked.connect(lambda: self._set_all_comp_metrics(True))
        clr_all.clicked.connect(lambda: self._set_all_comp_metrics(False))
        btn_row.addWidget(sel_all)
        btn_row.addWidget(clr_all)
        btn_row.addStretch()
        comp_outer_lay.addLayout(btn_row)

        # Scroll area wrapping a two-column grid
        scroll_comp = QScrollArea()
        scroll_comp.setWidgetResizable(True)
        scroll_comp.setFrameShape(QScrollArea.Shape.NoFrame)
        scroll_comp.setMinimumHeight(180)
        scroll_comp.setMaximumHeight(320)

        inner = QWidget()
        # Two-column grid: sections alternate left/right columns
        grid = QGridLayout(inner)
        grid.setSpacing(4)
        grid.setContentsMargins(4, 4, 4, 4)
        grid.setColumnStretch(0, 1)
        grid.setColumnStretch(1, 1)

        self._comp_metric_cbs: dict[str, QCheckBox] = {}
        # Lay out sections: first half in column 0, second half in column 1
        sections = _COMPONENT_METRICS
        mid = (len(sections) + 1) // 2  # split point
        col_widgets: list[list[QWidget]] = [[], []]  # widgets per column

        for s_idx, (section_name, items) in enumerate(sections):
            col = 0 if s_idx < mid else 1
            sec_lbl = QLabel(f"<b>{section_name}</b>")
            sec_lbl.setStyleSheet(f"color: {style.C_DARK_PURPLE}; margin-top: 4px;")
            col_widgets[col].append(sec_lbl)
            for key, desc in items:
                cb = QCheckBox(desc)
                cb.setChecked(key == "binary")
                cb.stateChanged.connect(self._on_comp_metric_changed)
                self._comp_metric_cbs[key] = cb
                col_widgets[col].append(cb)
            sep = _hline()
            col_widgets[col].append(sep)

        # Place widgets into the grid, one column at a time
        for col_idx, widgets in enumerate(col_widgets):
            for row_idx, widget in enumerate(widgets):
                grid.addWidget(widget, row_idx, col_idx)
            # Push remaining rows down with a stretch row
            grid.setRowStretch(len(widgets), 1)

        scroll_comp.setWidget(inner)
        comp_outer_lay.addWidget(scroll_comp)
        lay.addWidget(comp_grp)

        # ── Molecule-level metrics (two-column grid) ───────────────────────
        mol_grp = QGroupBox("Molecule-level metrics")
        mol_grid_lay = QGridLayout(mol_grp)
        mol_grid_lay.setSpacing(4)
        mol_grid_lay.setColumnStretch(0, 1)
        mol_grid_lay.setColumnStretch(1, 1)
        self._mol_metric_cbs: dict[str, QCheckBox] = {}
        for i, (key, desc) in enumerate(_MOLECULE_METRICS):
            cb = QCheckBox(desc)
            cb.setChecked(False)
            self._mol_metric_cbs[key] = cb
            mol_grid_lay.addWidget(cb, i // 2, i % 2)
        lay.addWidget(mol_grp)

        # Aggregation (shown only when graph metrics selected)
        self._agg_grp = QGroupBox("Aggregation (for multiple components per group)")
        agg_lay = QHBoxLayout(self._agg_grp)
        self._agg_bg = QButtonGroup(self)
        rb_mean = QRadioButton("Mean")
        rb_median = QRadioButton("Median")
        rb_mean.setChecked(True)
        self._agg_bg.addButton(rb_mean, 0)
        self._agg_bg.addButton(rb_median, 1)
        agg_lay.addWidget(rb_mean)
        agg_lay.addWidget(rb_median)
        agg_lay.addStretch()
        self._agg_grp.setVisible(False)
        lay.addWidget(self._agg_grp)

        return outer

    def _build_long_panel(self) -> QWidget:
        outer = QWidget()
        lay = QVBoxLayout(outer)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(8)

        # Group filter
        grp_filter = QGroupBox("Group filter")
        gf_lay = QFormLayout(grp_filter)
        self._long_group_sel = QComboBox()
        for key, label in _GROUP_SELECTIONS:
            self._long_group_sel.addItem(label, key)
        gf_lay.addRow("Groups:", self._long_group_sel)
        lay.addWidget(grp_filter)

        # Columns
        cols_grp = QGroupBox("Columns to include per row")
        cols_lay = QVBoxLayout(cols_grp)
        self._long_cbs: dict[str, QCheckBox] = {}
        for key, label in [
            ("name",        "Name / identifier"),
            ("smiles",      "SMILES"),
            ("group_name",  "Group name"),
            ("group_id",    "Group ID"),
            ("match_count", "Match count (number of matched components)"),
        ]:
            cb = QCheckBox(label)
            cb.setChecked(True)
            self._long_cbs[key] = cb
            cols_lay.addWidget(cb)
        lay.addWidget(cols_grp)
        return outer

    def _build_modelling_section(self, root: QVBoxLayout) -> None:
        grp = QGroupBox("Tables to export")
        lay = QVBoxLayout(grp)
        self._cb_scores = QCheckBox("Scores (fingerprint set vs ROC-AUC, MCC, Bal. Acc., Avg. Prec.)")
        self._cb_scores.setChecked(True)
        self._cb_bayes = QCheckBox("Bayesian t-test (pairwise comparisons)")
        self._cb_bayes.setChecked(True)
        lay.addWidget(self._cb_scores)
        lay.addWidget(self._cb_bayes)
        note = QLabel("Excel: both tables as separate sheets.  "
                      "Other formats: separate files (second file gets '_bayes' suffix).")
        note.setObjectName("muted")
        note.setWordWrap(True)
        lay.addWidget(note)
        root.addWidget(grp)

    # ── Slots ──────────────────────────────────────────────────────────────

    def _on_layout_changed(self, btn_id: int, checked: bool) -> None:
        if not checked:
            return
        self._flat_panel.setVisible(btn_id == 0)
        self._wide_panel.setVisible(btn_id == 1)
        self._long_panel.setVisible(btn_id == 2)

    def _on_comp_metric_changed(self) -> None:
        any_graph = any(
            self._comp_metric_cbs[k].isChecked()
            for k in _GRAPH_METRICS
            if k in self._comp_metric_cbs
        )
        self._agg_grp.setVisible(any_graph)

    def _set_all_comp_metrics(self, checked: bool) -> None:
        for cb in self._comp_metric_cbs.values():
            cb.setChecked(checked)

    # ── Export ─────────────────────────────────────────────────────────────

    def _selected_format(self) -> str:
        labels = list(_FORMAT_EXTS.keys())
        return labels[self._fmt_group.checkedId()]

    def _do_export(self) -> None:
        fmt = self._selected_format()
        ext = _FORMAT_EXTS[fmt]

        try:
            if self._context == "results":
                self._export_results(fmt, ext)
            elif self._context == "prioritisation":
                self._export_prioritisation(fmt, ext)
            elif self._context == "modelling":
                self._export_modelling(fmt, ext)
        except Exception as exc:
            QMessageBox.critical(self, "Export error", str(exc))

    def _pick_path(self, default_name: str, ext: str, fmt: str) -> str | None:
        if fmt == "CSV":
            filt = "CSV files (*.csv);;All files (*)"
        elif fmt == "Excel (.xlsx)":
            filt = "Excel files (*.xlsx);;All files (*)"
        elif fmt == "TSV":
            filt = "TSV files (*.tsv);;All files (*)"
        else:
            filt = "JSON files (*.json);;All files (*)"
        path, _ = QFileDialog.getSaveFileName(
            self, "Save export", default_name + ext, filt
        )
        return path or None

    # Results export ──────────────────────────────────────────────────────

    def _export_results(self, fmt: str, ext: str) -> None:
        layout_id = self._layout_group.checkedId()

        if layout_id == 0:  # Flat
            cols = {k for k, cb in self._flat_cbs.items() if cb.isChecked()}
            if not cols:
                QMessageBox.warning(self, "Nothing selected",
                                    "Please select at least one column to export.")
                return
            df = _build_flat(self._data, cols)

        elif layout_id == 1:  # Wide
            comp_metrics = [k for k, cb in self._comp_metric_cbs.items() if cb.isChecked()]
            mol_metrics = [k for k, cb in self._mol_metric_cbs.items() if cb.isChecked()]
            if not comp_metrics and not mol_metrics:
                QMessageBox.warning(self, "Nothing selected",
                                    "Please select at least one metric to include.")
                return
            group_sel = self._wide_group_sel.currentData()
            agg = "median" if self._agg_bg.checkedId() == 1 else "mean"
            df = _build_wide(self._data, group_sel, comp_metrics, mol_metrics, agg)

        else:  # Long
            cols = {k for k, cb in self._long_cbs.items() if cb.isChecked()}
            if not cols:
                QMessageBox.warning(self, "Nothing selected",
                                    "Please select at least one column to export.")
                return
            group_sel = self._long_group_sel.currentData()
            df = _build_long(self._data, group_sel, cols)
            if df.empty:
                QMessageBox.information(self, "No rows",
                                        "No group matches found for the selected group filter.")
                return

        path = self._pick_path("pfas_results", ext, fmt)
        if not path:
            return
        _write_df(df, path, fmt)
        QMessageBox.information(self, "Export complete",
                                f"Saved {len(df):,} rows to:\n{path}")
        self.accept()

    # Prioritisation export ───────────────────────────────────────────────

    def _export_prioritisation(self, fmt: str, ext: str) -> None:
        if not self._data:
            QMessageBox.information(self, "No data", "No prioritisation results to export.")
            return
        df = _build_prioritisation(self._data)
        path = self._pick_path("pfas_prioritisation", ext, fmt)
        if not path:
            return
        _write_df(df, path, fmt)
        QMessageBox.information(self, "Export complete",
                                f"Saved {len(df):,} rows to:\n{path}")
        self.accept()

    # Modelling export ────────────────────────────────────────────────────

    def _export_modelling(self, fmt: str, ext: str) -> None:
        export_scores = self._cb_scores.isChecked()
        export_bayes = self._cb_bayes.isChecked()
        if not export_scores and not export_bayes:
            QMessageBox.warning(self, "Nothing selected",
                                "Please select at least one table to export.")
            return

        scores_df: pd.DataFrame | None = None
        bayes_df: pd.DataFrame | None = None

        if export_scores and self._data.get("scores") is not None:
            scores_df = self._data["scores"].reset_index()
            scores_df.columns = ["fingerprint_set"] + list(scores_df.columns[1:])

        if export_bayes and self._data.get("bayes") is not None:
            bayes_df = pd.DataFrame(self._data["bayes"])

        path = self._pick_path("pfas_modelling", ext, fmt)
        if not path:
            return

        if fmt == "Excel (.xlsx)":
            with pd.ExcelWriter(path, engine="openpyxl") as writer:
                if scores_df is not None:
                    scores_df.to_excel(writer, sheet_name="Scores", index=False)
                if bayes_df is not None:
                    bayes_df.to_excel(writer, sheet_name="Bayesian t-test", index=False)
        else:
            p = Path(path)
            n_written = 0
            if scores_df is not None:
                _write_df(scores_df, str(p), fmt)
                n_written += 1
            if bayes_df is not None:
                bayes_path = str(p.with_stem(p.stem + "_bayes"))
                _write_df(bayes_df, bayes_path, fmt)
                n_written += 1

        msg = "Export complete."
        if fmt != "Excel (.xlsx)" and export_scores and export_bayes and scores_df is not None and bayes_df is not None:
            msg += (f"\nScores → {path}"
                    f"\nBayesian t-test → {Path(path).with_stem(Path(path).stem + '_bayes')}")
        else:
            msg += f"\nSaved to: {path}"
        QMessageBox.information(self, "Export complete", msg)
        self.accept()
