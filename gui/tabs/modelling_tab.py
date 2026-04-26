"""
Tab 6 — ML Modelling

Benchmark multiple fingerprint sets (PFASGroups presets, Morgan, ToxPrint,
TxP_PFAS, and optionally a custom TSV) on a user-specified binary target
column.  Uses HistGradientBoostingClassifier + RepeatedStratifiedKFold with
Bayesian correlated t-test comparison.

pyCSRML (ToxPrint + TxP_PFAS) is used only when installed; a warning banner
is shown otherwise.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QSplitter,
    QLabel, QPushButton, QComboBox, QSpinBox, QCheckBox,
    QGroupBox, QFormLayout, QLineEdit, QFileDialog,
    QTableWidget, QTableWidgetItem, QHeaderView,
    QProgressBar, QMessageBox, QTabWidget,
    QDialog, QDialogButtonBox, QScrollArea,
)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg

from gui import style
from gui.workers import ModelWorker
from gui.utils.fingerprints import is_pycsrml_available
from gui.utils.export_dialog import ExportDialog
from PFASGroups.embeddings import FINGERPRINT_PRESETS


class ModellingTab(QWidget):
    """Tab 6: ML modelling benchmark."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._embedding_set = None
        self._df = None
        self._worker = None
        self._results_data: dict | None = None
        self._build_ui()

    # ── Build ─────────────────────────────────────────────────────────────

    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(14, 14, 14, 10)
        root.setSpacing(10)

        title = QLabel("Step 6 · ML Modelling Benchmark")
        title.setObjectName("section_title")
        root.addWidget(title)

        note = QLabel(
            "Compare fingerprint descriptors for binary property prediction using "
            "HistGradientBoosting with repeated stratified k-fold CV. "
            "Pairwise Bayesian correlated t-tests (ROPE = 0.01 ROC-AUC) are reported."
        )
        note.setObjectName("muted")
        note.setWordWrap(True)
        root.addWidget(note)

        if not is_pycsrml_available():
            warn = QLabel(
                "pyCSRML is not installed.  ToxPrint and TxP_PFAS fingerprints are "
                "unavailable.  Install pyCSRML to enable them."
            )
            warn.setStyleSheet("color: #7a4300; background: #fff3cd; padding: 8px; "
                               "border-radius: 4px;")
            warn.setWordWrap(True)
            root.addWidget(warn)

        splitter = QSplitter(Qt.Orientation.Horizontal)

        # ── Left: configuration ───────────────────────────────────────────
        left = QWidget()
        left.setMaximumWidth(360)
        left_lay = QVBoxLayout(left)
        left_lay.setContentsMargins(0, 0, 8, 0)
        left_lay.setSpacing(10)

        # Target column
        target_grp = QGroupBox("Target (binary label)")
        target_form = QFormLayout(target_grp)
        self._target_combo = QComboBox()
        self._target_combo.addItem("— choose column —")
        target_form.addRow("Label column:", self._target_combo)
        lbl = QLabel("Load data in Tab 1 first; the column must be 0/1 or True/False.")
        lbl.setObjectName("muted")
        lbl.setWordWrap(True)
        target_form.addRow(lbl)
        left_lay.addWidget(target_grp)

        # Fingerprint sets
        fp_grp = QGroupBox("Fingerprint sets")
        fp_lay = QVBoxLayout(fp_grp)

        self._fp_summary_lbl = QLabel()
        self._fp_summary_lbl.setWordWrap(True)
        self._fp_summary_lbl.setObjectName("muted")
        fp_lay.addWidget(self._fp_summary_lbl)

        cfg_btn = QPushButton("Configure fingerprints…")
        cfg_btn.setObjectName("btn_secondary")
        cfg_btn.clicked.connect(self._open_fp_dialog)
        fp_lay.addWidget(cfg_btn)
        left_lay.addWidget(fp_grp)

        # Default selection state (mirrors old defaults)
        self._fp_state: dict = {
            "presets": {k: (k == "best") for k in FINGERPRINT_PRESETS},
            "morgan": True,
            "toxprint": False,
            "txp_pfas": False,
            "custom": False,
            "custom_path": "",
        }
        self._update_fp_summary()

        # CV settings
        cv_grp = QGroupBox("Cross-validation settings")
        cv_form = QFormLayout(cv_grp)
        self._cv_splits = QSpinBox()
        self._cv_splits.setRange(2, 10)
        self._cv_splits.setValue(3)
        cv_form.addRow("k-folds:", self._cv_splits)
        self._cv_repeats = QSpinBox()
        self._cv_repeats.setRange(1, 20)
        self._cv_repeats.setValue(5)
        cv_form.addRow("Repeats:", self._cv_repeats)
        left_lay.addWidget(cv_grp)

        self._progress = QProgressBar()
        self._progress.setValue(0)
        self._progress.setFixedHeight(20)
        left_lay.addWidget(self._progress)

        self._run_btn = QPushButton("Run Benchmark")
        self._run_btn.setObjectName("btn_primary")
        self._run_btn.setEnabled(False)
        self._run_btn.clicked.connect(self._run)
        left_lay.addWidget(self._run_btn)
        left_lay.addStretch()

        # ── Right: results ────────────────────────────────────────────────
        right = QWidget()
        right_lay = QVBoxLayout(right)
        right_lay.setContentsMargins(8, 0, 0, 0)
        right_lay.setSpacing(8)

        results_tabs = QTabWidget()
        results_tabs.setDocumentMode(True)

        # Scores table
        self._scores_table = QTableWidget(0, 5)
        self._scores_table.setHorizontalHeaderLabels(
            ["Fingerprint set", "ROC-AUC", "MCC", "Bal. Acc.", "Avg. Prec."]
        )
        self._scores_table.horizontalHeader().setSectionResizeMode(
            0, QHeaderView.ResizeMode.Stretch
        )
        self._scores_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._scores_table.setAlternatingRowColors(True)
        results_tabs.addTab(self._scores_table, "Scores")

        # Bayes table
        self._bayes_table = QTableWidget(0, 6)
        self._bayes_table.setHorizontalHeaderLabels(
            ["A", "B", "P(A>B)", "P(ROPE)", "P(B>A)", "Mean Δ"]
        )
        for col in range(6):
            self._bayes_table.horizontalHeader().setSectionResizeMode(
                col, QHeaderView.ResizeMode.ResizeToContents
            )
        self._bayes_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._bayes_table.setAlternatingRowColors(True)
        results_tabs.addTab(self._bayes_table, "Bayesian t-test")

        # ROC-AUC chart
        chart_widget = QWidget()
        chart_lay = QVBoxLayout(chart_widget)
        self._fig, self._ax = plt.subplots(figsize=(7, 3))
        self._canvas = FigureCanvasQTAgg(self._fig)
        chart_lay.addWidget(self._canvas)
        results_tabs.addTab(chart_widget, "ROC-AUC Chart")

        right_lay.addWidget(results_tabs, stretch=1)

        self._export_btn = QPushButton("Export\u2026")
        self._export_btn.setObjectName("btn_secondary")
        self._export_btn.setEnabled(False)
        self._export_btn.clicked.connect(self._export)
        right_lay.addWidget(self._export_btn)

        splitter.addWidget(left)
        splitter.addWidget(right)
        splitter.setSizes([350, 530])
        root.addWidget(splitter, stretch=1)

    # ── Public API ────────────────────────────────────────────────────────

    def set_results(self, embedding_set, df=None):
        """Called by MainWindow when classification finishes."""
        self._embedding_set = embedding_set
        self._df = df
        self._target_combo.blockSignals(True)
        self._target_combo.clear()
        self._target_combo.addItem("— choose column —")
        if df is not None:
            for col in df.columns:
                self._target_combo.addItem(col)
        self._target_combo.blockSignals(False)
        self._run_btn.setEnabled(True)

    # ── Slots ─────────────────────────────────────────────────────────────

    def _open_fp_dialog(self):
        dlg = _FpSelectDialog(self._fp_state, parent=self)
        if dlg.exec() == QDialog.DialogCode.Accepted:
            self._fp_state = dlg.get_state()
            self._update_fp_summary()

    def _update_fp_summary(self):
        parts = []
        for key, checked in self._fp_state["presets"].items():
            if checked:
                parts.append(f"PFASGroups:{key}")
        if self._fp_state.get("morgan"):
            parts.append("Morgan")
        if self._fp_state.get("toxprint"):
            parts.append("ToxPrint")
        if self._fp_state.get("txp_pfas"):
            parts.append("TxP_PFAS")
        if self._fp_state.get("custom") and self._fp_state.get("custom_path"):
            parts.append("Custom")
        if parts:
            self._fp_summary_lbl.setText(", ".join(parts))
        else:
            self._fp_summary_lbl.setText("None selected — click 'Configure' to choose.")

    def _run(self):
        if self._embedding_set is None:
            return

        target_idx = self._target_combo.currentIndex()
        if target_idx == 0:
            QMessageBox.warning(self, "No target",
                                "Please choose a target (label) column.")
            return

        target_col = self._target_combo.currentText()
        if self._df is None or target_col not in self._df.columns:
            QMessageBox.warning(self, "Column missing",
                                f"Column '{target_col}' not found in the data.")
            return

        try:
            y = self._df[target_col].astype(int).values
        except Exception as exc:
            QMessageBox.critical(self, "Target error",
                                 f"Cannot convert target column to int:\n{exc}")
            return

        smiles_list = [str(m.get("smiles") or "") for m in self._embedding_set]

        # Build feature sets
        from gui.utils.fingerprints import (
            get_pfasgroups_fingerprints, get_morgan_fingerprints,
            get_toxprint_fingerprints, get_txp_pfas_fingerprints,
            load_custom_fingerprints,
        )

        feature_sets = {}
        for key, checked in self._fp_state["presets"].items():
            if checked:
                try:
                    feature_sets[f"PFASGroups:{key}"] = get_pfasgroups_fingerprints(
                        self._embedding_set, preset=key
                    )
                except Exception as exc:
                    QMessageBox.warning(self, f"Fingerprint error ({key})", str(exc))

        if self._fp_state.get("morgan"):
            feature_sets["Morgan (r=2, 512b)"] = get_morgan_fingerprints(smiles_list)

        if self._fp_state.get("toxprint") and is_pycsrml_available():
            arr = get_toxprint_fingerprints(smiles_list)
            if arr is not None:
                feature_sets["ToxPrint (729b)"] = arr

        if self._fp_state.get("txp_pfas") and is_pycsrml_available():
            arr = get_txp_pfas_fingerprints(smiles_list)
            if arr is not None:
                feature_sets["TxP_PFAS (129b)"] = arr

        if self._fp_state.get("custom"):
            path = self._fp_state.get("custom_path", "").strip()
            if path:
                arr = load_custom_fingerprints(path, smiles_list)
                if arr is not None:
                    feature_sets["Custom"] = arr
                else:
                    QMessageBox.warning(self, "Custom fingerprint error",
                                        "Failed to load custom fingerprint file.")

        if not feature_sets:
            QMessageBox.warning(self, "No fingerprints",
                                "No fingerprint sets selected.")
            return

        self._progress.setValue(0)
        self._run_btn.setEnabled(False)

        worker = ModelWorker(
            feature_sets=feature_sets,
            y=y,
            cv_splits=self._cv_splits.value(),
            cv_repeats=self._cv_repeats.value(),
        )
        worker.progress.connect(self._progress.setValue)
        worker.result.connect(self._show_results)
        worker.error.connect(self._on_error)
        worker.finished.connect(lambda: self._run_btn.setEnabled(True))
        self._worker = worker
        worker.start()

    def _show_results(self, payload: dict):
        self._results_data = payload
        self._export_btn.setEnabled(True)
        summary_df = payload["scores"]
        bayes = payload["bayes"]

        # ── Scores table
        self._scores_table.setRowCount(len(summary_df))
        for row_i, (idx, row) in enumerate(summary_df.iterrows()):
            self._scores_table.setItem(row_i, 0, QTableWidgetItem(str(idx)))
            for col_i, metric in enumerate(["roc_auc", "mcc", "bal_acc", "avg_prec"]):
                val = row.get(metric, float("nan"))
                item = QTableWidgetItem(f"{val:.4f}")
                item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                self._scores_table.setItem(row_i, col_i + 1, item)

        # ── Bayes table
        self._bayes_table.setRowCount(len(bayes))
        for row_i, bt in enumerate(bayes):
            self._bayes_table.setItem(row_i, 0, QTableWidgetItem(bt.get("A", "")))
            self._bayes_table.setItem(row_i, 1, QTableWidgetItem(bt.get("B", "")))
            for col_i, key in enumerate(["p_A_wins", "p_rope", "p_B_wins", "mean_diff"]):
                val = bt.get(key, float("nan"))
                item = QTableWidgetItem(f"{val:.4f}" if isinstance(val, float) else str(val))
                item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                self._bayes_table.setItem(row_i, col_i + 2, item)

        # ── ROC-AUC chart
        names = list(summary_df.index)
        vals = [summary_df.loc[n, "roc_auc"] for n in names]
        self._ax.clear()
        colours = [style.C_BLUE if "PFASGroups" in n else style.C_ORANGE for n in names]
        self._ax.barh(names, vals, color=colours, alpha=0.85)
        self._ax.set_xlabel("Mean ROC-AUC")
        self._ax.set_title("Fingerprint benchmark – ROC-AUC", fontsize=9)
        self._ax.set_xlim(0, 1)
        self._ax.axvline(0.5, color="grey", lw=0.8, ls="--")
        self._fig.tight_layout(pad=0.8)
        self._canvas.draw()

    def _on_error(self, msg: str):
        QMessageBox.critical(self, "Benchmark error", msg)
        self._run_btn.setEnabled(True)

    def _export(self):
        if not self._results_data:
            QMessageBox.information(self, "No Data", "Run the benchmark first.")
            return
        dlg = ExportDialog(self._results_data, "modelling", parent=self)
        dlg.exec()


def _hline() -> QWidget:
    line = QWidget()
    line.setFixedHeight(1)
    line.setStyleSheet("background: #ddd;")
    return line


# ── Fingerprint selection dialog ───────────────────────────────────────────────

class _FpSelectDialog(QDialog):
    """Pop-up for selecting which fingerprint sets to benchmark."""

    def __init__(self, state: dict, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Select Fingerprint Sets")
        self.setMinimumWidth(400)
        self._build_ui(state)

    def _build_ui(self, state: dict):
        lay = QVBoxLayout(self)
        lay.setSpacing(10)

        # Scroll area for presets (can be long)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        inner = QWidget()
        inner_lay = QVBoxLayout(inner)
        inner_lay.setSpacing(6)

        grp_presets = QGroupBox("PFASGroups presets")
        preset_lay = QVBoxLayout(grp_presets)
        self._preset_cbs: dict[str, QCheckBox] = {}
        for key in FINGERPRINT_PRESETS:
            cb = QCheckBox(key)
            cb.setChecked(state["presets"].get(key, False))
            self._preset_cbs[key] = cb
            preset_lay.addWidget(cb)
        inner_lay.addWidget(grp_presets)

        grp_other = QGroupBox("Other fingerprints")
        other_lay = QVBoxLayout(grp_other)

        _has_pycsrml = is_pycsrml_available()

        self._cb_morgan = QCheckBox("Morgan (RDKit, r=2, 512 bits)")
        self._cb_morgan.setChecked(state.get("morgan", True))
        other_lay.addWidget(self._cb_morgan)

        self._cb_toxprint = QCheckBox("ToxPrint (729 bits, requires pyCSRML)")
        self._cb_toxprint.setChecked(state.get("toxprint", False))
        self._cb_toxprint.setEnabled(_has_pycsrml)
        if not _has_pycsrml:
            self._cb_toxprint.setToolTip("Install pyCSRML to enable")
        other_lay.addWidget(self._cb_toxprint)

        self._cb_txp_pfas = QCheckBox("TxP_PFAS (129 bits, requires pyCSRML)")
        self._cb_txp_pfas.setChecked(state.get("txp_pfas", False))
        self._cb_txp_pfas.setEnabled(_has_pycsrml)
        if not _has_pycsrml:
            self._cb_txp_pfas.setToolTip("Install pyCSRML to enable")
        other_lay.addWidget(self._cb_txp_pfas)

        other_lay.addWidget(_hline())

        custom_row = QHBoxLayout()
        self._cb_custom = QCheckBox("Custom TSV/CSV:")
        self._cb_custom.setChecked(state.get("custom", False))
        custom_row.addWidget(self._cb_custom)
        self._custom_path = QLineEdit(state.get("custom_path", ""))
        self._custom_path.setPlaceholderText("path/to/fingerprints.tsv")
        self._custom_path.setReadOnly(True)
        custom_row.addWidget(self._custom_path, stretch=1)
        browse_btn = QPushButton("Browse…")
        browse_btn.clicked.connect(self._browse_custom)
        custom_row.addWidget(browse_btn)
        other_lay.addLayout(custom_row)

        inner_lay.addWidget(grp_other)
        scroll.setWidget(inner)
        lay.addWidget(scroll)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        lay.addWidget(buttons)

    def _browse_custom(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select fingerprint file",
            filter="Fingerprint files (*.tsv *.csv);;All (*.*)"
        )
        if path:
            self._custom_path.setText(path)
            self._cb_custom.setChecked(True)

    def get_state(self) -> dict:
        return {
            "presets": {k: cb.isChecked() for k, cb in self._preset_cbs.items()},
            "morgan": self._cb_morgan.isChecked(),
            "toxprint": self._cb_toxprint.isChecked(),
            "txp_pfas": self._cb_txp_pfas.isChecked(),
            "custom": self._cb_custom.isChecked(),
            "custom_path": self._custom_path.text(),
        }
