"""
Tab 5 — Chemical Space

Plots an interactive 2-D projection of the classified embedding set using
UMAP, PCA, or t-SNE.  Each point is coloured by PFAS group or a user-supplied
label column.  The chart is rendered as a self-contained Plotly HTML page in
a QWebEngineView.
"""
from __future__ import annotations

import tempfile
import webbrowser

from PySide6.QtCore import Qt, QUrl
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QSplitter,
    QLabel, QPushButton, QComboBox, QDoubleSpinBox, QSpinBox,
    QGroupBox, QFormLayout, QProgressBar, QMessageBox,
)

try:
    from PySide6.QtWebEngineWidgets import QWebEngineView
    _HAS_WEBENGINE = True
except ImportError:
    _HAS_WEBENGINE = False

from gui import style
from gui.workers import ChemSpaceWorker
from PFASGroups.embeddings import FINGERPRINT_PRESETS


class ChemSpaceTab(QWidget):
    """Tab 5: Chemical space visualisation."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._embedding_set = None
        self._df_labels = None
        self._worker = None
        self._build_ui()

    # ── Build ─────────────────────────────────────────────────────────────

    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(14, 14, 14, 10)
        root.setSpacing(10)

        title = QLabel("Step 5 · Chemical Space")
        title.setObjectName("section_title")
        root.addWidget(title)

        note = QLabel(
            "Reduce the PFASGroups fingerprint using UMAP, PCA, or t-SNE and "
            "visualise the 2-D chemical space as an interactive Plotly chart."
        )
        note.setObjectName("muted")
        note.setWordWrap(True)
        root.addWidget(note)

        if not _HAS_WEBENGINE:
            warn = QLabel(
                "PyQt6-WebEngine is not installed. Interactive charts are unavailable.\n"
                "Install it with: pip install PyQt6-WebEngine"
            )
            warn.setStyleSheet("color: #b00; background: #fff3cd; padding: 8px; "
                               "border-radius: 4px;")
            warn.setWordWrap(True)
            root.addWidget(warn)

        splitter = QSplitter(Qt.Orientation.Horizontal)

        # ── Left: controls ────────────────────────────────────────────────
        left = QWidget()
        left.setMaximumWidth(320)
        left_lay = QVBoxLayout(left)
        left_lay.setContentsMargins(0, 0, 8, 0)
        left_lay.setSpacing(10)

        proj_grp = QGroupBox("Projection")
        proj_form = QFormLayout(proj_grp)

        self._method_combo = QComboBox()
        self._method_combo.addItems(["UMAP", "PCA", "t-SNE"])
        self._method_combo.currentTextChanged.connect(self._on_method_changed)
        proj_form.addRow("Method:", self._method_combo)

        self._preset_combo = QComboBox()
        for key, info in FINGERPRINT_PRESETS.items():
            self._preset_combo.addItem(key, userData=info)
        proj_form.addRow("Fingerprint preset:", self._preset_combo)

        self._label_combo = QComboBox()
        self._label_combo.addItem("Top PFAS group (auto)")
        proj_form.addRow("Colour by:", self._label_combo)

        left_lay.addWidget(proj_grp)

        # UMAP params
        self._umap_grp = QGroupBox("UMAP parameters")
        umap_form = QFormLayout(self._umap_grp)

        self._n_neighbours = QSpinBox()
        self._n_neighbours.setRange(2, 200)
        self._n_neighbours.setValue(15)
        umap_form.addRow("n_neighbours:", self._n_neighbours)

        self._min_dist = QDoubleSpinBox()
        self._min_dist.setRange(0.0, 1.0)
        self._min_dist.setSingleStep(0.05)
        self._min_dist.setValue(0.1)
        self._min_dist.setDecimals(2)
        umap_form.addRow("min_dist:", self._min_dist)
        left_lay.addWidget(self._umap_grp)

        # t-SNE params
        self._tsne_grp = QGroupBox("t-SNE parameters")
        tsne_form = QFormLayout(self._tsne_grp)
        self._perplexity = QDoubleSpinBox()
        self._perplexity.setRange(5.0, 200.0)
        self._perplexity.setValue(30.0)
        self._perplexity.setDecimals(1)
        tsne_form.addRow("perplexity:", self._perplexity)
        self._tsne_grp.setVisible(False)
        left_lay.addWidget(self._tsne_grp)

        self._progress = QProgressBar()
        self._progress.setValue(0)
        self._progress.setFixedHeight(20)
        left_lay.addWidget(self._progress)

        self._plot_btn = QPushButton("Plot")
        self._plot_btn.setObjectName("btn_primary")
        self._plot_btn.setEnabled(False)
        self._plot_btn.clicked.connect(self._run)
        left_lay.addWidget(self._plot_btn)
        left_lay.addStretch()

        # ── Right: chart view ─────────────────────────────────────────────
        right = QWidget()
        right_lay = QVBoxLayout(right)
        right_lay.setContentsMargins(8, 0, 0, 0)

        if _HAS_WEBENGINE:
            self._web_view = QWebEngineView()
            self._web_view.setHtml("<body style='background:#f5f5f5'>"
                                   "<p style='padding:20px;color:#888'>Run the "
                                   "projection to see the chemical space chart.</p>"
                                   "</body>")
            right_lay.addWidget(self._web_view)
        else:
            no_web = QWidget()
            no_web_lay = QVBoxLayout(no_web)
            no_web_lay.setAlignment(Qt.AlignmentFlag.AlignCenter)
            info = QLabel(
                "PySide6-WebEngineWidgets is not installed.\n"
                "After plotting, the chart will open in your default browser.\n\n"
                "To enable in-app display, install: pip install PySide6-WebEngineWidgets"
            )
            info.setAlignment(Qt.AlignmentFlag.AlignCenter)
            info.setWordWrap(True)
            info.setStyleSheet("color: #555; padding: 12px;")
            self._open_browser_btn = QPushButton("Open last chart in browser")
            self._open_browser_btn.setEnabled(False)
            self._open_browser_btn.clicked.connect(self._open_in_browser)
            no_web_lay.addWidget(info)
            no_web_lay.addWidget(self._open_browser_btn)
            right_lay.addWidget(no_web)

        splitter.addWidget(left)
        splitter.addWidget(right)
        splitter.setSizes([300, 600])
        root.addWidget(splitter, stretch=1)

    # ── Public API ────────────────────────────────────────────────────────

    def set_results(self, embedding_set, df_labels=None):
        """Called by MainWindow after classification."""
        self._embedding_set = embedding_set
        self._df_labels = df_labels

        # Populate label column combo from df_labels
        self._label_combo.blockSignals(True)
        self._label_combo.clear()
        self._label_combo.addItem("Top PFAS group (auto)")
        if df_labels is not None:
            for col in df_labels.columns:
                self._label_combo.addItem(col)
        self._label_combo.blockSignals(False)

        self._plot_btn.setEnabled(True)

    # ── Slots ─────────────────────────────────────────────────────────────

    def _on_method_changed(self, method: str):
        is_umap = method == "UMAP"
        is_tsne = method == "t-SNE"
        self._umap_grp.setVisible(is_umap)
        self._tsne_grp.setVisible(is_tsne)

    def _run(self):
        if self._embedding_set is None:
            return

        method = self._method_combo.currentText()
        preset = self._preset_combo.currentText()
        label_idx = self._label_combo.currentIndex()
        label_col = None if label_idx == 0 else self._label_combo.currentText()

        self._progress.setValue(0)
        self._plot_btn.setEnabled(False)

        worker = ChemSpaceWorker(
            embedding_set=self._embedding_set,
            method=method,
            preset=preset,
            label_col=label_col,
            df_labels=self._df_labels,
            n_neighbours=self._n_neighbours.value(),
            min_dist=self._min_dist.value(),
        )
        worker.progress.connect(self._progress.setValue)
        worker.result.connect(self._show_html)
        worker.error.connect(self._on_error)
        worker.finished.connect(lambda: self._plot_btn.setEnabled(True))
        self._worker = worker
        worker.start()

    def _show_html(self, html_string: str):
        self._progress.setValue(100)
        if _HAS_WEBENGINE:
            self._web_view.setHtml(html_string)
        else:
            # Save to a temp file and open in browser
            tmp = tempfile.NamedTemporaryFile(
                mode="w", suffix=".html", delete=False, encoding="utf-8"
            )
            tmp.write(html_string)
            tmp.close()
            self._last_html_path = tmp.name
            self._open_browser_btn.setEnabled(True)
            webbrowser.open(f"file:///{tmp.name}")

    def _open_in_browser(self):
        if hasattr(self, "_last_html_path"):
            webbrowser.open(f"file:///{self._last_html_path}")

    def _on_error(self, msg: str):
        QMessageBox.critical(self, "Chemical space error", msg)
        self._plot_btn.setEnabled(True)
