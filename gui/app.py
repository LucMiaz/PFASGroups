"""
PFASGroups GUI — Main window.
"""
from __future__ import annotations

import sys
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtGui import QIcon, QPixmap, QFont
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QTabWidget, QWidget,
    QHBoxLayout, QVBoxLayout, QLabel, QStatusBar,
    QSizePolicy,
)

from gui import style

_LOGO_PATH = Path(__file__).parent.parent / "logo" / "PFASGroups_logo.png"
_TASKBAR_LOGO_PATH = Path(__file__).parent.parent / "logo" / "taskbar_logo.svg"


class HeaderWidget(QWidget):
    """Dark purple header bar with logo and app title."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("header_widget")
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.setFixedHeight(70)

        layout = QHBoxLayout(self)
        layout.setContentsMargins(16, 8, 16, 8)
        layout.setSpacing(14)

        # Logo
        if _LOGO_PATH.exists():
            pix = QPixmap(str(_LOGO_PATH)).scaledToHeight(
                48, Qt.TransformationMode.SmoothTransformation
            )
            logo_lbl = QLabel()
            logo_lbl.setPixmap(pix)
            logo_lbl.setFixedSize(pix.size())
            layout.addWidget(logo_lbl)
        else:
            placeholder = QLabel("PFASGroups")
            placeholder.setStyleSheet(
                f"color: {style.C_ORANGE}; font-size: 22px; font-weight: 800;"
            )
            layout.addWidget(placeholder)

        # Title block
        title_block = QVBoxLayout()
        title_block.setSpacing(2)

        title_lbl = QLabel("PFASGroups")
        title_lbl.setObjectName("header_title")
        title_block.addWidget(title_lbl)

        sub_lbl = QLabel(
            "Classification · Screening · Chemical Space · Modelling"
        )
        sub_lbl.setObjectName("header_subtitle")
        title_block.addWidget(sub_lbl)

        layout.addLayout(title_block)
        layout.addStretch()

        version_lbl = QLabel("v3.2")
        version_lbl.setStyleSheet(
            f"color: {style.C_PURPLE_LIGHT}; font-size: 11px;"
        )
        layout.addWidget(version_lbl)


class MainWindow(QMainWindow):
    """Application main window."""

    def __init__(self):
        super().__init__()
        self.setWindowTitle("PFASGroups")
        self.resize(1280, 900)
        self.setMinimumSize(900, 650)

        if _TASKBAR_LOGO_PATH.exists():
            self.setWindowIcon(QIcon(str(_TASKBAR_LOGO_PATH)))
        elif _LOGO_PATH.exists():
            self.setWindowIcon(QIcon(str(_LOGO_PATH)))

        # ── Central widget ─────────────────────────────────────────────────
        central = QWidget()
        self.setCentralWidget(central)
        root_layout = QVBoxLayout(central)
        root_layout.setContentsMargins(0, 0, 0, 0)
        root_layout.setSpacing(0)

        # Header
        root_layout.addWidget(HeaderWidget())

        # Tab widget
        self._tabs = QTabWidget()
        self._tabs.setDocumentMode(True)
        root_layout.addWidget(self._tabs)

        # ── Status bar ─────────────────────────────────────────────────────
        self._status_bar = QStatusBar()
        self.setStatusBar(self._status_bar)
        self._status_bar.showMessage("Ready.")

        # ── Build tabs ─────────────────────────────────────────────────────
        self._build_tabs()

    # ------------------------------------------------------------------
    def _build_tabs(self):
        from gui.tabs.classify_tab import ClassifyTab
        from gui.tabs.results_tab import ResultsTab
        from gui.tabs.definition_tester_tab import DefinitionTesterTab
        from gui.tabs.prioritise_tab import PrioritiseTab
        from gui.tabs.chemspace_tab import ChemSpaceTab
        from gui.tabs.modelling_tab import ModellingTab
        from gui.tabs.help_tab import HelpTab

        self._classify_tab = ClassifyTab()
        self._results_tab = ResultsTab()
        self._def_tab = DefinitionTesterTab()
        self._prior_tab = PrioritiseTab()
        self._chem_tab = ChemSpaceTab()
        self._model_tab = ModellingTab()
        self._help_tab = HelpTab()

        # Wire classification → results
        self._classify_tab.classification_done.connect(self._on_classification_done)

        self._tabs.addTab(self._classify_tab,  "1 · Classification")
        self._tabs.addTab(self._results_tab,   "2 · Results")
        self._tabs.addTab(self._def_tab,       "3 · Definition Tester")
        self._tabs.addTab(self._prior_tab,     "4 · Prioritisation")
        self._tabs.addTab(self._chem_tab,      "5 · Chemical Space")
        self._tabs.addTab(self._model_tab,     "6 · Modelling")
        self._tabs.addTab(self._help_tab,      "7 · Help")

    # ------------------------------------------------------------------
    def _on_classification_done(self, embedding_set):
        """Propagate classification results to all dependent tabs."""
        self._results_tab.set_results(embedding_set)
        self._prior_tab.set_results(embedding_set)
        self._chem_tab.set_results(embedding_set)
        self._model_tab.set_results(embedding_set)
        self._tabs.setCurrentIndex(1)  # jump to Results tab

    # ------------------------------------------------------------------
    def show_status(self, msg: str):
        self._status_bar.showMessage(msg)


def create_app(argv=None) -> tuple[QApplication, MainWindow]:
    app = QApplication(argv or sys.argv)
    app.setApplicationName("PFASGroups")
    app.setOrganizationName("Stockholm University")
    style.apply(app)
    window = MainWindow()
    return app, window
