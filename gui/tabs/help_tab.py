"""
Tab 7 — Help

Displays the bundled help.html file in a QTextBrowser (no WebEngine
dependency required).  Falls back gracefully if the file is missing.
"""
from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt, QUrl
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel,
    QPushButton, QTextBrowser, QSizePolicy,
)

from gui import style

_HELP_HTML = Path(__file__).parent.parent / "data" / "help.html"


class HelpTab(QWidget):
    """Tab 7: embedded documentation."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._build_ui()

    # ── Build ─────────────────────────────────────────────────────────────

    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(14, 14, 14, 10)
        root.setSpacing(10)

        # ── Header row ────────────────────────────────────────────────────
        header = QHBoxLayout()

        title = QLabel("Step 7 · Documentation")
        title.setObjectName("section_title")
        header.addWidget(title)
        header.addStretch()

        # "Open in browser" button (handy for printing / zooming)
        open_btn = QPushButton("Open in Browser")
        open_btn.setObjectName("secondary_button")
        open_btn.setFixedWidth(140)
        open_btn.clicked.connect(self._open_in_browser)
        header.addWidget(open_btn)

        root.addLayout(header)

        # ── Browser ───────────────────────────────────────────────────────
        self._browser = QTextBrowser()
        self._browser.setOpenExternalLinks(True)
        self._browser.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding
        )
        # Make the background consistent with the help html's --bg colour
        self._browser.setStyleSheet(
            "QTextBrowser {"
            f"  background: #f7f7fa;"
            f"  border: 1px solid {style.C_BORDER};"
            "  border-radius: 6px;"
            "  padding: 8px;"
            "}"
        )
        root.addWidget(self._browser)

        self._load_html()

    # ── Helpers ───────────────────────────────────────────────────────────

    def _load_html(self):
        if _HELP_HTML.exists():
            # Use setSource so relative resource paths inside the HTML
            # (images, etc.) can be resolved correctly.
            self._browser.setSource(QUrl.fromLocalFile(str(_HELP_HTML)))
        else:
            self._browser.setHtml(
                "<h2 style='color:#51127C'>Help file not found</h2>"
                f"<p>Expected: <code>{_HELP_HTML}</code></p>"
                "<p>Please re-install or re-clone the PFASGroups repository.</p>"
            )

    def _open_in_browser(self):
        import webbrowser
        if _HELP_HTML.exists():
            webbrowser.open(_HELP_HTML.as_uri())
        else:
            from PySide6.QtWidgets import QMessageBox
            QMessageBox.warning(
                self,
                "Help File Not Found",
                f"Cannot locate {_HELP_HTML}.\n"
                "Please re-install the PFASGroups package.",
            )
