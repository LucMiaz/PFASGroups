"""
Tab 4 — Prioritise Molecules

Rank the classified embedding set by either:
  • Reference-based KL divergence (supply a reference CSV/Excel/SQLite)
  • Intrinsic metrics: total component size, max component size, or their sum

Controls mirror prioritise_molecules() + get_priority_statistics() API.
"""
from __future__ import annotations

import io

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QSplitter,
    QLabel, QPushButton, QComboBox, QDoubleSpinBox, QSpinBox,
    QGroupBox, QFormLayout, QLineEdit, QFileDialog,
    QTableWidget, QTableWidgetItem, QHeaderView,
    QProgressBar, QMessageBox, QScrollArea,
)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg

from gui import style
from gui.workers import PrioritiseWorker


class PrioritiseTab(QWidget):
    """Tab 4: molecule prioritisation."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._embedding_set = None
        self._worker = None
        self._build_ui()

    # ── Build ─────────────────────────────────────────────────────────────

    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(14, 14, 14, 10)
        root.setSpacing(10)

        title = QLabel("Step 4 · Prioritise Molecules")
        title.setObjectName("section_title")
        root.addWidget(title)

        note = QLabel(
            "Rank molecules by structural diversity against a reference set (KL divergence) "
            "or by intrinsic metrics (total/max component size)."
        )
        note.setObjectName("muted")
        note.setWordWrap(True)
        root.addWidget(note)

        splitter = QSplitter(Qt.Orientation.Horizontal)

        # ── Left: options ─────────────────────────────────────────────────
        left = QWidget()
        left.setMaximumWidth(380)
        left_lay = QVBoxLayout(left)
        left_lay.setContentsMargins(0, 0, 8, 0)
        left_lay.setSpacing(10)

        # Mode
        mode_grp = QGroupBox("Prioritisation mode")
        mode_form = QFormLayout(mode_grp)

        self._mode_combo = QComboBox()
        self._mode_combo.addItems([
            "Intrinsic: total component size",
            "Intrinsic: max component size",
            "Intrinsic: total + max",
            "Reference-based KL divergence",
        ])
        self._mode_combo.currentIndexChanged.connect(self._on_mode_changed)
        mode_form.addRow("Mode:", self._mode_combo)

        # Reference file (hidden unless reference mode)
        self._ref_row = QWidget()
        ref_row_lay = QHBoxLayout(self._ref_row)
        ref_row_lay.setContentsMargins(0, 0, 0, 0)
        self._ref_path = QLineEdit()
        self._ref_path.setPlaceholderText("Path to reference file…")
        self._ref_path.setReadOnly(True)
        browse_btn = QPushButton("Browse…")
        browse_btn.setFixedWidth(80)
        browse_btn.clicked.connect(self._browse_reference)
        ref_row_lay.addWidget(self._ref_path, stretch=1)
        ref_row_lay.addWidget(browse_btn)
        mode_form.addRow("Reference file:", self._ref_row)
        self._ref_row.setVisible(False)

        left_lay.addWidget(mode_grp)

        # Weighting
        weight_grp = QGroupBox("Weighting (intrinsic modes)")
        weight_form = QFormLayout(weight_grp)

        self._a_spin = QDoubleSpinBox()
        self._a_spin.setRange(0.0, 1.0)
        self._a_spin.setSingleStep(0.1)
        self._a_spin.setValue(1.0)
        self._a_spin.setDecimals(2)
        weight_form.addRow("α (total weight):", self._a_spin)

        self._b_spin = QDoubleSpinBox()
        self._b_spin.setRange(0.0, 1.0)
        self._b_spin.setSingleStep(0.1)
        self._b_spin.setValue(1.0)
        self._b_spin.setDecimals(2)
        weight_form.addRow("β (max weight):", self._b_spin)

        self._percentile_spin = QSpinBox()
        self._percentile_spin.setRange(1, 99)
        self._percentile_spin.setValue(90)
        weight_form.addRow("Percentile cutoff:", self._percentile_spin)

        lbl_w = QLabel("Only used when mode is 'Intrinsic: total + max'.")
        lbl_w.setObjectName("muted")
        lbl_w.setWordWrap(True)
        weight_form.addRow(lbl_w)
        left_lay.addWidget(weight_grp)

        # Run
        self._progress = QProgressBar()
        self._progress.setValue(0)
        self._progress.setFixedHeight(20)
        left_lay.addWidget(self._progress)

        self._run_btn = QPushButton("Run Prioritisation")
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

        results_lbl = QLabel("Ranking")
        results_lbl.setObjectName("section_title")
        right_lay.addWidget(results_lbl)

        self._table = QTableWidget(0, 4)
        self._table.setHorizontalHeaderLabels(["Rank", "Name", "SMILES", "Score"])
        self._table.horizontalHeader().setSectionResizeMode(
            2, QHeaderView.ResizeMode.Stretch
        )
        self._table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._table.setAlternatingRowColors(True)
        self._table.setSortingEnabled(True)
        right_lay.addWidget(self._table, stretch=2)

        # Chart area
        chart_lbl = QLabel("Priority score distribution")
        chart_lbl.setObjectName("section_title")
        right_lay.addWidget(chart_lbl)

        self._fig, self._ax = plt.subplots(figsize=(7, 2.4))
        self._fig.tight_layout(pad=1)
        self._canvas = FigureCanvasQTAgg(self._fig)
        self._canvas.setMinimumHeight(150)
        right_lay.addWidget(self._canvas, stretch=1)

        splitter.addWidget(left)
        splitter.addWidget(right)
        splitter.setSizes([360, 520])

        root.addWidget(splitter, stretch=1)

    # ── Public API ────────────────────────────────────────────────────────

    def set_results(self, embedding_set):
        """Called by MainWindow when classification finishes."""
        self._embedding_set = embedding_set
        self._run_btn.setEnabled(True)

    # ── Slots ─────────────────────────────────────────────────────────────

    def _on_mode_changed(self, idx: int):
        self._ref_row.setVisible(idx == 3)

    def _browse_reference(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select reference file",
            filter="Data files (*.csv *.xlsx *.xls *.db *.sqlite);;All (*.*)"
        )
        if path:
            self._ref_path.setText(path)

    def _run(self):
        if self._embedding_set is None:
            return

        mode_idx = self._mode_combo.currentIndex()
        mode_map = {
            0: "total",
            1: "max",
            2: "combined",
            3: "reference",
        }
        mode = mode_map[mode_idx]

        reference = None
        if mode == "reference":
            ref_path = self._ref_path.text().strip()
            if not ref_path:
                QMessageBox.warning(self, "No reference",
                                    "Please choose a reference file.")
                return
            try:
                from gui.utils.io_readers import read_file
                reference = read_file(ref_path)
            except Exception as exc:
                QMessageBox.critical(self, "Reference file error", str(exc))
                return

        self._progress.setValue(0)
        self._run_btn.setEnabled(False)

        worker = PrioritiseWorker(
            embedding_set=self._embedding_set,
            reference=reference,
            mode=mode,
            a=self._a_spin.value(),
            b=self._b_spin.value(),
            percentile=self._percentile_spin.value(),
        )
        worker.progress.connect(self._progress.setValue)
        worker.result.connect(self._show_results)
        worker.error.connect(self._on_error)
        worker.finished.connect(lambda: self._run_btn.setEnabled(True))
        self._worker = worker
        worker.start()

    def _show_results(self, payload):
        """payload = {"ranked": list[tuple(name, smiles, score)], "statistics": dict}"""
        ranked = payload.get("ranked", [])
        stats = payload.get("statistics", {})

        # Populate table
        self._table.setSortingEnabled(False)
        self._table.setRowCount(len(ranked))
        for i, (name, smiles, score) in enumerate(ranked):
            self._table.setItem(i, 0, _num_item(i + 1))
            self._table.setItem(i, 1, QTableWidgetItem(str(name)))
            self._table.setItem(i, 2, QTableWidgetItem(str(smiles)))
            self._table.setItem(i, 3, _num_item(float(score)))
        self._table.setSortingEnabled(True)
        self._table.resizeColumnToContents(0)
        self._table.resizeColumnToContents(1)
        self._table.resizeColumnToContents(3)

        # Chart
        scores = [float(s) for _, _, s in ranked]
        self._ax.clear()
        if scores:
            self._ax.barh(
                range(min(len(scores), 30)),
                scores[:30],
                color=style.C_BLUE,
                alpha=0.8,
            )
            self._ax.set_xlabel("Priority score")
            self._ax.set_title("Top-30 priority scores", fontsize=9)
            self._ax.invert_yaxis()
        self._fig.tight_layout(pad=0.8)
        self._canvas.draw()

    def _on_error(self, msg: str):
        QMessageBox.critical(self, "Prioritisation error", msg)
        self._run_btn.setEnabled(True)


def _num_item(val) -> QTableWidgetItem:
    item = QTableWidgetItem()
    item.setData(Qt.ItemDataRole.DisplayRole, val)
    return item
