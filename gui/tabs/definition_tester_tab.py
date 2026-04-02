"""
Tab 3 — Definition Tester

Live testing of:
  (a) built-in or custom PFASDefinitions (SMARTS / fluorineRatio)
  (b) custom HalogenGroup definitions (ComponentSmarts + SMARTS JSON)

Mirrors the functionality of PFASgroupsJS smarts_tester.html and
pfas_definition_checker.html.
"""
from __future__ import annotations

import json

from PySide6.QtCore import Qt, QByteArray
from PySide6.QtSvgWidgets import QSvgWidget
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QSplitter,
    QLabel, QPushButton, QTextEdit, QLineEdit,
    QGroupBox, QFormLayout, QComboBox, QCheckBox,
    QScrollArea, QFrame, QMessageBox, QTabWidget,
    QProgressBar,
)

from gui import style
from gui.workers import DefinitionTestWorker, SmartsTestWorker


class DefinitionTesterTab(QWidget):
    """Tab 3: test PFAS definitions and custom SMARTS."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._worker = None
        self._build_ui()

    # ── Build ────────────────────────────────────────────────────────────

    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(14, 14, 14, 10)
        root.setSpacing(10)

        title = QLabel("Step 3 · Definition & SMARTS Tester")
        title.setObjectName("section_title")
        root.addWidget(title)

        note = QLabel(
            "Test SMILES against built-in PFAS definitions or define custom "
            "PFASDefinitions and ComponentSmarts to validate them."
        )
        note.setObjectName("muted")
        note.setWordWrap(True)
        root.addWidget(note)

        tabs = QTabWidget()
        tabs.setDocumentMode(True)
        tabs.addTab(self._build_definitions_subtab(), "PFAS Definitions")
        tabs.addTab(self._build_group_subtab(), "Custom Group (SMARTS)")
        root.addWidget(tabs, stretch=1)

    # ── Sub-tab A: PFAS Definitions ──────────────────────────────────────

    def _build_definitions_subtab(self) -> QWidget:
        w = QWidget()
        layout = QVBoxLayout(w)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(10)

        splitter = QSplitter(Qt.Orientation.Horizontal)

        # Left: SMILES input + custom definition JSON
        left = QWidget()
        left_lay = QVBoxLayout(left)
        left_lay.setContentsMargins(0, 0, 6, 0)

        smiles_grp = QGroupBox("Test SMILES")
        smiles_form = QFormLayout(smiles_grp)
        self._def_smiles = QLineEdit()
        self._def_smiles.setPlaceholderText("e.g. OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F")
        smiles_form.addRow("SMILES:", self._def_smiles)

        mol_row = QHBoxLayout()
        self._def_mol_svg = QSvgWidget()
        self._def_mol_svg.setFixedSize(250, 200)
        self._def_mol_svg.setStyleSheet("background: white; border: 1px solid #ddd; border-radius: 4px;")
        mol_row.addWidget(self._def_mol_svg)
        mol_row.addStretch()
        smiles_form.addRow("Structure:", mol_row)
        left_lay.addWidget(smiles_grp)

        # Custom definitions JSON
        custom_grp = QGroupBox("Custom PFASDefinitions JSON (optional)")
        custom_lay = QVBoxLayout(custom_grp)
        self._custom_def_json = QTextEdit()
        self._custom_def_json.setPlaceholderText(
            '[\n  {\n    "id": 99,\n    "name": "My Custom Definition",\n'
            '    "smarts": ["[CX4](F)(F)(F)"],\n    "fluorineRatio": null,\n'
            '    "description": "CF3 group"\n  }\n]'
        )
        self._custom_def_json.setMaximumHeight(160)
        custom_lay.addWidget(self._custom_def_json)
        lbl = QLabel("Leave blank to test against all built-in definitions.")
        lbl.setObjectName("muted")
        custom_lay.addWidget(lbl)
        left_lay.addWidget(custom_grp)

        run_row = QHBoxLayout()
        self._def_progress = QProgressBar()
        self._def_progress.setValue(0)
        self._def_progress.setFixedHeight(18)
        run_row.addWidget(self._def_progress, stretch=1)

        test_btn = QPushButton("Test")
        test_btn.setObjectName("btn_primary")
        test_btn.clicked.connect(self._run_definition_test)
        run_row.addWidget(test_btn)
        left_lay.addLayout(run_row)
        left_lay.addStretch()

        # Right: results
        right = QWidget()
        right_lay = QVBoxLayout(right)
        right_lay.setContentsMargins(6, 0, 0, 0)

        res_lbl = QLabel("Results")
        res_lbl.setObjectName("section_title")
        right_lay.addWidget(res_lbl)

        self._def_results_scroll = QScrollArea()
        self._def_results_scroll.setWidgetResizable(True)
        self._def_results_widget = QWidget()
        self._def_results_layout = QVBoxLayout(self._def_results_widget)
        self._def_results_layout.setSpacing(6)
        self._def_results_layout.addStretch()
        self._def_results_scroll.setWidget(self._def_results_widget)
        right_lay.addWidget(self._def_results_scroll, stretch=1)

        splitter.addWidget(left)
        splitter.addWidget(right)
        splitter.setSizes([440, 400])
        layout.addWidget(splitter, stretch=1)

        # Connect SMILES → preview
        self._def_smiles.textChanged.connect(self._update_def_preview)

        return w

    # ── Sub-tab B: Custom Group ───────────────────────────────────────────

    def _build_group_subtab(self) -> QWidget:
        w = QWidget()
        layout = QVBoxLayout(w)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(10)

        splitter = QSplitter(Qt.Orientation.Horizontal)

        # Left: group definition
        left = QWidget()
        left_lay = QVBoxLayout(left)
        left_lay.setContentsMargins(0, 0, 6, 0)

        grp_grp = QGroupBox("HalogenGroup Definition")
        grp_form = QFormLayout(grp_grp)
        grp_form.setSpacing(8)

        self._grp_name = QLineEdit()
        self._grp_name.setPlaceholderText("Custom PFCA")
        grp_form.addRow("Group name:", self._grp_name)

        self._grp_component_smarts = QComboBox()
        self._grp_component_smarts.setEditable(True)
        self._grp_component_smarts.addItems(["Perfluoroalkyl", "Polyfluoroalkyl", ""])
        grp_form.addRow("componentSmarts:", self._grp_component_smarts)

        self._grp_smarts_json = QTextEdit()
        self._grp_smarts_json.setPlaceholderText(
            '{"[CX3](=O)[OX2H1]": 1}'
        )
        self._grp_smarts_json.setMaximumHeight(90)
        grp_form.addRow("SMARTS dict (JSON):", self._grp_smarts_json)

        self._grp_max_dist = QLineEdit("0")
        grp_form.addRow("max_dist_from_comp:", self._grp_max_dist)

        self._grp_constraints_json = QTextEdit()
        self._grp_constraints_json.setPlaceholderText('{"gte": {"F": 2}}')
        self._grp_constraints_json.setMaximumHeight(70)
        grp_form.addRow("Constraints (JSON):", self._grp_constraints_json)

        left_lay.addWidget(grp_grp)

        test_grp = QGroupBox("Test Molecule")
        test_form = QFormLayout(test_grp)
        self._grp_smiles = QLineEdit()
        self._grp_smiles.setPlaceholderText("SMILES to test")
        test_form.addRow("SMILES:", self._grp_smiles)
        left_lay.addWidget(test_grp)

        run_row2 = QHBoxLayout()
        self._grp_progress = QProgressBar()
        self._grp_progress.setValue(0)
        self._grp_progress.setFixedHeight(18)
        run_row2.addWidget(self._grp_progress, stretch=1)

        grp_test_btn = QPushButton("Test Group")
        grp_test_btn.setObjectName("btn_primary")
        grp_test_btn.clicked.connect(self._run_group_test)
        run_row2.addWidget(grp_test_btn)
        left_lay.addLayout(run_row2)
        left_lay.addStretch()

        # Right: diagnostics
        right = QWidget()
        right_lay = QVBoxLayout(right)
        right_lay.setContentsMargins(6, 0, 0, 0)

        diag_lbl = QLabel("Diagnostics")
        diag_lbl.setObjectName("section_title")
        right_lay.addWidget(diag_lbl)

        self._grp_result_text = QTextEdit()
        self._grp_result_text.setReadOnly(True)
        self._grp_result_text.setPlaceholderText("Test result will appear here…")
        right_lay.addWidget(self._grp_result_text, stretch=1)

        splitter.addWidget(left)
        splitter.addWidget(right)
        splitter.setSizes([440, 400])
        layout.addWidget(splitter, stretch=1)

        return w

    # ── Slots ─────────────────────────────────────────────────────────────

    def _update_def_preview(self, smiles: str):
        if not smiles.strip():
            return
        try:
            from gui.utils.mol_renderer import smiles_to_svg
            svg_bytes = smiles_to_svg(smiles.strip(), 250, 200)
            if svg_bytes:
                self._def_mol_svg.load(QByteArray(svg_bytes))
        except Exception:
            pass

    def _run_definition_test(self):
        smiles = self._def_smiles.text().strip()
        if not smiles:
            QMessageBox.warning(self, "No SMILES", "Please enter a SMILES string.")
            return

        custom_json = self._custom_def_json.toPlainText().strip()
        custom_defs = None
        if custom_json:
            try:
                raw = json.loads(custom_json)
                from PFASGroups import PFASDefinition
                custom_defs = [PFASDefinition(**d) for d in raw]
            except Exception as exc:
                QMessageBox.critical(self, "JSON Error",
                                     f"Invalid custom definition JSON:\n{exc}")
                return

        self._def_progress.setValue(0)
        self._worker = DefinitionTestWorker(smiles, custom_definitions=custom_defs)
        self._worker.progress.connect(self._def_progress.setValue)
        self._worker.result.connect(self._show_definition_results)
        self._worker.error.connect(lambda e: QMessageBox.critical(self, "Error", e))
        self._worker.start()

    def _show_definition_results(self, results: list[dict]):
        # Clear old results
        while self._def_results_layout.count() > 1:
            item = self._def_results_layout.takeAt(0)
            if item and item.widget():
                item.widget().deleteLater()

        for r in results:
            row = QHBoxLayout()
            badge = QLabel("PASS" if r["passed"] else "FAIL")
            badge.setObjectName("badge_pass" if r["passed"] else "badge_fail")
            badge.setFixedWidth(40)
            row.addWidget(badge)

            name_lbl = QLabel(f"<b>{r['id']}. {r['name']}</b>")
            name_lbl.setToolTip(r.get("description", ""))
            row.addWidget(name_lbl, stretch=1)

            container = QWidget()
            container.setLayout(row)
            self._def_results_layout.insertWidget(
                self._def_results_layout.count() - 1, container
            )

    def _run_group_test(self):
        smiles = self._grp_smiles.text().strip()
        if not smiles:
            QMessageBox.warning(self, "No SMILES", "Please enter a SMILES string.")
            return

        try:
            smarts_dict = json.loads(self._grp_smarts_json.toPlainText() or "{}")
            constraints = json.loads(self._grp_constraints_json.toPlainText() or "{}")
        except Exception as exc:
            QMessageBox.critical(self, "JSON Error", f"Invalid JSON:\n{exc}")
            return

        group_dict = {
            "id": 9999,
            "name": self._grp_name.text() or "Custom Group",
            "smarts": smarts_dict or None,
            "componentSmarts": self._grp_component_smarts.currentText() or None,
            "constraints": constraints,
            "max_dist_from_comp": int(self._grp_max_dist.text() or "0"),
        }

        self._grp_progress.setValue(0)
        worker = SmartsTestWorker(smiles, group_dict)
        worker.progress.connect(self._grp_progress.setValue)
        worker.result.connect(self._show_group_result)
        worker.error.connect(lambda e: (
            self._grp_result_text.setPlainText(f"ERROR: {e}"),
        ))
        worker.start()
        self._worker = worker

    def _show_group_result(self, result):
        if isinstance(result, dict):
            text = json.dumps(result, indent=2, default=str)
        else:
            text = str(result)
        self._grp_result_text.setPlainText(text)
