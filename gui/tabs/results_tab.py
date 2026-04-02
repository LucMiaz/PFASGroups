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

from PySide6.QtCore import Qt, QByteArray
from PySide6.QtGui import QPixmap, QImage
from PySide6.QtSvgWidgets import QSvgWidget
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QScrollArea, QGroupBox, QCheckBox, QFileDialog, QMessageBox,
    QFrame, QSizePolicy, QGridLayout, QLineEdit,
)

from gui import style


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
        self.setFixedWidth(300)

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

        self._hide_btn = QPushButton("Hide")
        self._hide_btn.setObjectName("btn_secondary")
        self._hide_btn.setFixedSize(46, 22)
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
            groups_lbl = QLabel()
            groups_lbl.setWordWrap(True)
            groups_lbl.setStyleSheet(f"color: {style.C_DARK_PURPLE}; font-size: 11px;")
            text = "<b>Matched groups:</b><br/>" + "<br/>".join(f"• {g}" for g in groups[:8])
            if len(groups) > 8:
                text += f"<br/><i>…and {len(groups) - 8} more</i>"
            groups_lbl.setText(text)
            layout.addWidget(groups_lbl)
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
        groups = []
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

    def _add_definitions_row(self, emb, layout):
        matched_def_ids = set()
        try:
            defs = emb.get("pfas_definitions") or emb.get("definitions") or []
            for d in defs:
                if isinstance(d, dict):
                    matched_def_ids.add(d.get("id"))
                else:
                    matched_def_ids.add(getattr(d, "id", None))
        except Exception:
            pass

        if matched_def_ids:
            def_row = QHBoxLayout()
            def_row.setSpacing(4)
            for did in sorted(matched_def_ids):
                if did is None:
                    continue
                badge = QLabel(f"DEF {did}")
                badge.setObjectName("badge_pass")
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
        self._hide_btn.setText("Show" if not visible else "Hide")

    def toggle(self):
        self._toggle_visibility()


class ResultsTab(QWidget):
    """Tab 2: show classified molecule cards."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._embedding_set = None
        self._cards: list[CompoundCard] = []
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

        show_all = QPushButton("Show All")
        show_all.setObjectName("btn_secondary")
        show_all.clicked.connect(lambda: self._set_all_visible(True))
        toolbar.addWidget(show_all)

        hide_all = QPushButton("Hide All")
        hide_all.setObjectName("btn_secondary")
        hide_all.clicked.connect(lambda: self._set_all_visible(False))
        toolbar.addWidget(hide_all)

        export_btn = QPushButton("Export CSV")
        export_btn.setObjectName("btn_secondary")
        export_btn.clicked.connect(self._export_csv)
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
        self._clear_grid()
        self._cards = []

        items = list(embedding_set) if embedding_set else []
        for i, emb in enumerate(items):
            card = CompoundCard(i, emb)
            self._cards.append(card)

        self._populate_grid(self._cards)
        self._count_lbl.setText(f"{len(self._cards):,} molecule(s) classified.")

    # ── Internal ─────────────────────────────────────────────────────────

    def _clear_grid(self):
        while self._grid.count():
            item = self._grid.takeAt(0)
            if item and item.widget():
                item.widget().deleteLater()

    def _populate_grid(self, cards: list[CompoundCard], n_cols: int = 4):
        self._clear_grid()
        if not cards:
            self._grid.addWidget(self._placeholder, 0, 0)
            return
        for i, card in enumerate(cards):
            row, col = divmod(i, n_cols)
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
        self._clear_grid()
        n_cols = 4
        visible = [c for c in self._cards if c.isVisible() and c != self._placeholder]
        for i, card in enumerate(visible):
            row, col = divmod(i, n_cols)
            self._grid.addWidget(card, row, col)

    def _set_all_visible(self, visible: bool):
        for card in self._cards:
            card.set_content_visible(visible)

    def _export_csv(self):
        if not self._embedding_set:
            QMessageBox.information(self, "No Data", "No results to export.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export Results", "pfas_results.csv",
            "CSV files (*.csv);;All files (*)"
        )
        if not path:
            return
        try:
            df = self._embedding_set.to_dataframe()
            df.to_csv(path, index=False)
            QMessageBox.information(self, "Export Complete",
                                    f"Results saved to:\n{path}")
        except AttributeError:
            # Fallback if to_dataframe() isn't available
            import pandas as pd
            rows = []
            for emb in self._embedding_set:
                rows.append({
                    "name": emb.get("name", ""),
                    "smiles": emb.get("smiles", ""),
                })
            pd.DataFrame(rows).to_csv(path, index=False)
            QMessageBox.information(self, "Export Complete",
                                    f"Results saved to:\n{path}")
        except Exception as exc:
            QMessageBox.critical(self, "Export Error", str(exc))
