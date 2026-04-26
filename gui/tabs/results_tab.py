"""
Tab 2 — Results

Scrollable grid of compound cards. Each card shows:
  - Compound name / SMILES
  - SVG structure drawing with highlighted matched atoms
  - List of matched PFAS groups
  - Show / hide toggle

Toolbar: Show All, Hide All, Export CSV.
"""
from __future__ import annotations

from PySide6.QtCore import Qt, QByteArray, QTimer
from PySide6.QtGui import QPixmap, QImage
from PySide6.QtSvgWidgets import QSvgWidget
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QScrollArea, QGroupBox, QCheckBox, QFileDialog, QMessageBox,
    QFrame, QSizePolicy, QGridLayout, QLineEdit,
)

from gui import style
from gui.utils.export_dialog import ExportDialog

_CARD_W = 340        # fixed card width in px
_CARD_SPACING = 12   # grid spacing in px
_MAX_COLS = 3        # maximum number of columns


class CompoundCard(QWidget):
    """Visual card for a single classified molecule."""

    def __init__(self, index: int, embedding, parent=None):
        super().__init__(parent)
        self._index = index
        self._embedding = embedding
        self._visible = True
        self._build(embedding)

    # ── Build ────────────────────────────────────────────────────────────

    def _build(self, emb):
        self.setObjectName("compound_card")
        self.setStyleSheet(
            f"QWidget#compound_card {{ background: {style.C_SURFACE}; "
            f"border: 1px solid {style.C_BORDER}; border-radius: 6px; }}"
        )
        self.setFixedWidth(_CARD_W)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(6)

        # Header row: name + hide button
        hdr = QHBoxLayout()
        name = emb.get("name") or emb.get("smiles", "")[:24] or f"Mol {self._index + 1}"
        name_lbl = QLabel(f"<b>{name}</b>")
        name_lbl.setToolTip(emb.get("smiles", ""))
        name_lbl.setWordWrap(True)
        hdr.addWidget(name_lbl, stretch=1)

        self._hide_btn = QPushButton("▼")  # ▼ down arrow
        self._hide_btn.setFixedSize(28, 28)
        self._hide_btn.setStyleSheet(
            "QPushButton { border: none; background: transparent; padding: 0;"
            f" color: {style.C_DARK_PURPLE}; font-size: 16px; }}"
            f"QPushButton:hover {{ color: {style.C_ORANGE}; }}"
        )
        self._hide_btn.clicked.connect(self._toggle_visibility)
        hdr.addWidget(self._hide_btn)
        layout.addLayout(hdr)

        # SVG structure
        self._mol_widget = QSvgWidget()
        self._mol_widget.setFixedSize(280, 220)
        self._mol_widget.setStyleSheet("background: white; border-radius: 4px;")
        self._load_structure(emb)
        layout.addWidget(self._mol_widget)

        # Matched groups
        groups = self._get_groups(emb)
        if groups:
            hdr_lbl = QLabel("<b>Matched groups:</b>")
            hdr_lbl.setStyleSheet(f"color: {style.C_DARK_PURPLE}; font-size: 11px;")
            layout.addWidget(hdr_lbl)

            groups_box = QWidget()
            groups_box.setStyleSheet("background: transparent;")
            box_lay = QVBoxLayout(groups_box)
            box_lay.setContentsMargins(0, 0, 0, 0)
            box_lay.setSpacing(2)
            for g in groups:
                item_lbl = QLabel(f"\u2022 {g}")
                item_lbl.setStyleSheet(
                    f"color: {style.C_DARK_PURPLE}; font-size: 11px; background: transparent;"
                )
                item_lbl.setWordWrap(True)
                box_lay.addWidget(item_lbl)
            box_lay.addStretch()

            groups_scroll = QScrollArea()
            groups_scroll.setWidget(groups_box)
            groups_scroll.setWidgetResizable(True)
            groups_scroll.setMaximumHeight(90)
            groups_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
            groups_scroll.setFrameShape(QFrame.Shape.NoFrame)
            groups_scroll.setStyleSheet("background: transparent;")
            layout.addWidget(groups_scroll)
        else:
            no_match = QLabel("<i>No PFAS groups matched.</i>")
            no_match.setStyleSheet(f"color: {style.C_TEXT_MUTED}; font-size: 11px;")
            layout.addWidget(no_match)

        # PFAS definitions
        self._add_definitions_row(emb, layout)

    def _load_structure(self, emb):
        try:
            from gui.utils.mol_renderer import embedding_to_svg
            svg_bytes = embedding_to_svg(emb, width=280, height=220)
            if svg_bytes:
                self._mol_widget.load(QByteArray(svg_bytes))
                return
        except Exception:
            pass
        # Fallback: plain SMILES rendering
        try:
            from gui.utils.mol_renderer import smiles_to_svg
            svg_bytes = smiles_to_svg(emb.get("smiles", ""), 280, 220)
            if svg_bytes:
                self._mol_widget.load(QByteArray(svg_bytes))
        except Exception:
            pass

    def _get_groups(self, emb) -> list[str]:
        """Return unique matched PFAS group names in match order."""
        seen: dict[str, int] = {}
        try:
            for m in emb.get("matches", []):
                if m.get("type") != "HalogenGroup":
                    continue
                name = m.get("group_name") or str(m.get("id") or "?")
                seen[name] = seen.get(name, 0) + 1
        except Exception:
            pass
        return list(seen.keys())

    def _add_definitions_row(self, emb, layout):
        def_ids: list = []
        try:
            for m in emb.get("matches", []):
                if m.get("type") == "PFASdefinition":
                    did = m.get("id")
                    if did is not None and did not in def_ids:
                        def_ids.append(did)
        except Exception:
            pass

        if def_ids:
            def_row = QHBoxLayout()
            def_row.setSpacing(4)
            for did in sorted(def_ids):
                badge = QLabel(f"{did}")
                badge.setObjectName("badge_pass")
                badge.setToolTip(f"PFAS Definition {did}")
                def_row.addWidget(badge)
            def_row.addStretch()
            layout.addLayout(def_row)

    # ── Toggle ────────────────────────────────────────────────────────────

    def _toggle_visibility(self):
        self.set_content_visible(not self._visible)

    def set_content_visible(self, visible: bool):
        self._visible = visible
        self._mol_widget.setVisible(visible)
        # Find and toggle group labels
        for i in range(self.layout().count()):
            item = self.layout().itemAt(i)
            if item and item.widget() and isinstance(item.widget(), QLabel) and item.widget() != self.layout().itemAt(0).layout().itemAt(0).widget() if self.layout().itemAt(0) else False:
                item.widget().setVisible(visible)
        self._hide_btn.setText("▶" if not visible else "▼")  # ▶ right arrow if hidden, ▼ down arrow if visible

    def toggle(self):
        self._toggle_visibility()


class ResultsTab(QWidget):
    """Tab 2: show classified molecule cards."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._embedding_set = None
        self._cards: list[CompoundCard] = []
        self._current_n_cols: int = 0
        self._build_ui()

    # ── Build ────────────────────────────────────────────────────────────

    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(12, 12, 12, 8)
        root.setSpacing(10)

        # Title row
        title_row = QHBoxLayout()
        title = QLabel("Step 2 · Classification Results")
        title.setObjectName("section_title")
        title_row.addWidget(title)
        title_row.addStretch()
        root.addLayout(title_row)

        # Toolbar
        toolbar = QHBoxLayout()
        toolbar.setSpacing(8)

        self._count_lbl = QLabel("No results loaded.")
        self._count_lbl.setObjectName("muted")
        toolbar.addWidget(self._count_lbl)
        toolbar.addStretch()

        self._search = QLineEdit()
        self._search.setPlaceholderText("Filter by name or SMILES…")
        self._search.setMaximumWidth(260)
        self._search.textChanged.connect(self._filter_cards)
        toolbar.addWidget(self._search)

        show_all = QPushButton("▶ Show All")
        show_all.setObjectName("btn_secondary")
        show_all.clicked.connect(lambda: self._set_all_visible(True))
        toolbar.addWidget(show_all)

        hide_all = QPushButton("▼ Hide All")
        hide_all.setObjectName("btn_secondary")
        hide_all.clicked.connect(lambda: self._set_all_visible(False))
        toolbar.addWidget(hide_all)

        export_btn = QPushButton("Export\u2026")
        export_btn.setObjectName("btn_secondary")
        export_btn.clicked.connect(self._export)
        toolbar.addWidget(export_btn)

        root.addLayout(toolbar)

        # Scroll area with grid
        self._scroll = QScrollArea()
        self._scroll.setWidgetResizable(True)
        self._scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)

        self._grid_container = QWidget()
        self._grid = QGridLayout(self._grid_container)
        self._grid.setSpacing(12)
        self._grid.setContentsMargins(4, 4, 4, 4)
        self._scroll.setWidget(self._grid_container)
        root.addWidget(self._scroll, stretch=1)

        self._placeholder = QLabel("Run classification to see results here.")
        self._placeholder.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._placeholder.setObjectName("muted")
        self._grid.addWidget(self._placeholder, 0, 0)

    # ── API ───────────────────────────────────────────────────────────────

    def set_results(self, embedding_set):
        self._embedding_set = embedding_set
        # Destroy old cards before replacing them.
        for card in self._cards:
            card.deleteLater()
        self._clear_grid()
        self._cards = []

        items = list(embedding_set) if embedding_set else []
        for i, emb in enumerate(items):
            card = CompoundCard(i, emb)
            self._cards.append(card)

        self._populate_grid(self._cards)
        # Defer a re-layout so the viewport has its real width after Qt layout.
        QTimer.singleShot(0, lambda: self._populate_grid(self._cards))
        self._count_lbl.setText(f"{len(self._cards):,} molecule(s) classified.")

    # ── Internal ─────────────────────────────────────────────────────────

    def _clear_grid(self):
        """Remove all items from the grid layout.
        CompoundCards stay parented to _grid_container (takeAt only removes
        the layout item; setParent(None) would make them top-level windows)."""
        while self._grid.count():
            item = self._grid.takeAt(0)
            if item and item.widget():
                w = item.widget()
                if not isinstance(w, CompoundCard):
                    w.deleteLater()

    def _n_cols(self) -> int:
        vp_w = self._scroll.viewport().width()
        return max(1, min(_MAX_COLS, vp_w // (_CARD_W + _CARD_SPACING)))

    def _populate_grid(self, cards: list[CompoundCard]):
        n = self._n_cols()
        self._current_n_cols = n
        self._clear_grid()
        if not cards:
            self._grid.addWidget(self._placeholder, 0, 0)
            return
        for i, card in enumerate(cards):
            row, col = divmod(i, n)
            self._grid.addWidget(card, row, col)

    def _filter_cards(self, text: str):
        text = text.lower()
        visible_cards = []
        for card in self._cards:
            emb = card._embedding
            name = (emb.get("name") or "").lower()
            smi = (emb.get("smiles") or "").lower()
            if not text or text in name or text in smi:
                visible_cards.append(card)
            card.setVisible(not text or text in name or text in smi)
        # Re-layout only visible cards
        visible = [c for c in self._cards if c.isVisible() and c != self._placeholder]
        self._populate_grid(visible)

    def _set_all_visible(self, visible: bool):
        for card in self._cards:
            card.set_content_visible(visible)

    def showEvent(self, event):
        super().showEvent(event)
        # Reflow when the tab is first made visible (viewport width now known).
        if self._cards:
            QTimer.singleShot(0, lambda: self._populate_grid(
                [c for c in self._cards if c.isVisible() and c != self._placeholder]
            ))

    def resizeEvent(self, event):
        super().resizeEvent(event)
        new_cols = self._n_cols()
        if new_cols != self._current_n_cols and self._cards:
            visible = [c for c in self._cards if c.isVisible() and c != self._placeholder]
            self._populate_grid(visible)

    def _export(self):
        if not self._embedding_set:
            QMessageBox.information(self, "No Data", "No results to export.")
            return
        dlg = ExportDialog(self._embedding_set, "results", parent=self)
        dlg.exec()
