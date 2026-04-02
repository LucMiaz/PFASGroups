"""
Tab 1 — Classification

Lets the user load a compound dataset (CSV / Excel / SQLite / MariaDB / PostgreSQL),
map the SMILES and name columns, configure classification options, and run.
Emits classification_done(PFASEmbeddingSet) when the run succeeds.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QFormLayout,
    QLabel, QPushButton, QComboBox, QCheckBox, QGroupBox,
    QFileDialog, QProgressBar, QSizePolicy, QSplitter,
    QTabWidget, QLineEdit, QSpinBox, QFrame, QScrollArea,
    QMessageBox, QTextEdit, QDialog, QDialogButtonBox,
    QRadioButton, QButtonGroup,
)

from gui import style
from gui.workers import ClassifyWorker

# Path where DB connection presets are stored
_DB_PRESETS_FILE = Path.home() / ".pfasgroups" / "db_presets.json"


def _simple_obfuscate(s: str) -> str:
    """Very lightweight XOR obfuscation to avoid storing plaintext passwords."""
    key = b"PFASGroupsGUI2025"
    b = s.encode()
    return "".join(format(b[i] ^ key[i % len(key)], "02x") for i in range(len(b)))


def _simple_deobfuscate(h: str) -> str:
    key = b"PFASGroupsGUI2025"
    try:
        b = bytes(int(h[i:i+2], 16) for i in range(0, len(h), 2))
        return "".join(chr(b[i] ^ key[i % len(key)]) for i in range(len(b)))
    except Exception:
        return ""


def _load_presets() -> list[dict]:
    if _DB_PRESETS_FILE.exists():
        try:
            return json.loads(_DB_PRESETS_FILE.read_text(encoding="utf-8"))
        except Exception:
            pass
    return []


def _save_presets(presets: list[dict]) -> None:
    _DB_PRESETS_FILE.parent.mkdir(parents=True, exist_ok=True)
    _DB_PRESETS_FILE.write_text(json.dumps(presets, indent=2), encoding="utf-8")


class ClassifyTab(QWidget):
    """Tab 1: load data + run classification."""

    classification_done = Signal(object)  # PFASEmbeddingSet

    def __init__(self, parent=None):
        super().__init__(parent)
        self._df = None
        self._worker = None
        self._last_embedding = None
        self._build_ui()

    # ── Build ────────────────────────────────────────────────────────────

    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(16, 16, 16, 12)
        root.setSpacing(12)

        # Section title
        title = QLabel("Step 1 · Load & Classify")
        title.setObjectName("section_title")
        root.addWidget(title)

        # Splitter: left = data source, right = options
        splitter = QSplitter(Qt.Orientation.Horizontal)
        root.addWidget(splitter, stretch=1)

        splitter.addWidget(self._build_source_panel())
        splitter.addWidget(self._build_options_panel())
        splitter.setSizes([520, 340])

        # Progress + Run row
        bottom = QHBoxLayout()
        bottom.setSpacing(10)

        self._progress = QProgressBar()
        self._progress.setValue(0)
        self._progress.setTextVisible(True)
        self._progress.setFixedHeight(22)
        bottom.addWidget(self._progress, stretch=1)

        self._status_lbl = QLabel("Ready.")
        self._status_lbl.setObjectName("muted")
        self._status_lbl.setMaximumWidth(300)
        bottom.addWidget(self._status_lbl)

        self._run_btn = QPushButton("Run Classification")
        self._run_btn.setObjectName("btn_primary")
        self._run_btn.setMinimumWidth(160)
        self._run_btn.clicked.connect(self._run)
        bottom.addWidget(self._run_btn)

        root.addLayout(bottom)

    # ── Source panel ─────────────────────────────────────────────────────

    def _build_source_panel(self) -> QWidget:
        panel = QGroupBox("Data Source")
        layout = QVBoxLayout(panel)
        layout.setSpacing(10)

        # Source type tabs
        self._src_tabs = QTabWidget()
        self._src_tabs.setDocumentMode(True)
        self._src_tabs.addTab(self._build_file_tab(), "File")
        self._src_tabs.addTab(self._build_db_tab(), "Database")
        self._src_tabs.addTab(self._build_manual_tab(), "Manual / Paste")
        layout.addWidget(self._src_tabs)

        # Column mapping
        col_group = QGroupBox("Column Mapping")
        col_form = QFormLayout(col_group)
        col_form.setSpacing(8)

        self._smiles_col = QComboBox()
        self._smiles_col.setEditable(True)
        col_form.addRow("SMILES column:", self._smiles_col)

        self._name_col = QComboBox()
        self._name_col.setEditable(True)
        self._name_col.addItem("(none)")
        col_form.addRow("Name / ID column:", self._name_col)

        layout.addWidget(col_group)

        # Preview label
        self._preview_lbl = QLabel("No data loaded.")
        self._preview_lbl.setObjectName("muted")
        layout.addWidget(self._preview_lbl)

        return panel

    def _build_file_tab(self) -> QWidget:
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setContentsMargins(8, 8, 8, 8)
        lay.setSpacing(8)

        row = QHBoxLayout()
        self._file_path = QLineEdit()
        self._file_path.setPlaceholderText("No file selected…")
        self._file_path.setReadOnly(True)
        row.addWidget(self._file_path, stretch=1)

        browse_btn = QPushButton("Browse…")
        browse_btn.setObjectName("btn_secondary")
        browse_btn.clicked.connect(self._browse_file)
        row.addWidget(browse_btn)
        lay.addLayout(row)

        hint = QLabel("Supported: CSV, Excel (.xlsx), SQLite (.db/.sqlite)")
        hint.setObjectName("muted")
        lay.addWidget(hint)

        # Excel sheet picker (hidden for non-Excel files)
        self._sheet_row_widget = QWidget()
        sheet_row = QFormLayout(self._sheet_row_widget)
        sheet_row.setContentsMargins(0, 0, 0, 0)
        self._sheet_combo = QComboBox()
        self._sheet_combo.setEditable(False)
        self._sheet_combo.currentIndexChanged.connect(self._on_sheet_changed)
        sheet_row.addRow("Sheet (Excel):", self._sheet_combo)
        self._sheet_row_widget.setVisible(False)
        lay.addWidget(self._sheet_row_widget)

        lay.addStretch()
        return w

    def _build_db_tab(self) -> QWidget:
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setContentsMargins(8, 8, 8, 8)
        lay.setSpacing(8)

        # ── Preset picker row ────────────────────────────────────────────
        preset_row = QHBoxLayout()
        preset_row.addWidget(QLabel("Connection:"))
        self._db_preset_combo = QComboBox()
        self._db_preset_combo.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        self._db_preset_combo.currentIndexChanged.connect(self._on_preset_selected)
        preset_row.addWidget(self._db_preset_combo, stretch=1)

        edit_btn = QPushButton("Edit / Add…")
        edit_btn.setObjectName("btn_secondary")
        edit_btn.clicked.connect(self._open_preset_editor)
        preset_row.addWidget(edit_btn)

        del_btn = QPushButton("Delete")
        del_btn.setObjectName("btn_secondary")
        del_btn.clicked.connect(self._delete_preset)
        preset_row.addWidget(del_btn)
        lay.addLayout(preset_row)

        # ── Active connection summary ────────────────────────────────────
        self._db_summary = QLabel("No connection selected.")
        self._db_summary.setObjectName("muted")
        self._db_summary.setWordWrap(True)
        lay.addWidget(self._db_summary)

        # Hidden fields that store the resolved config for the active preset
        self._db_preset_data: dict = {}

        # ── SQL / Table row ──────────────────────────────────────────────
        sql_grp = QGroupBox("Query")
        sql_form = QFormLayout(sql_grp)
        sql_form.setSpacing(6)

        self._db_table = QLineEdit()
        self._db_table.setPlaceholderText("table name")
        sql_form.addRow("Table:", self._db_table)

        self._db_query = QTextEdit()
        self._db_query.setPlaceholderText("Optional: custom SQL query (overrides table)")
        self._db_query.setMaximumHeight(56)
        sql_form.addRow("SQL query:", self._db_query)
        lay.addWidget(sql_grp)

        load_btn = QPushButton("Load from Database")
        load_btn.setObjectName("btn_secondary")
        load_btn.clicked.connect(self._load_from_db)
        lay.addWidget(load_btn)
        lay.addStretch()

        self._refresh_preset_combo()
        return w

    def _build_manual_tab(self) -> QWidget:
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setContentsMargins(8, 8, 8, 8)
        lay.setSpacing(8)

        lbl = QLabel("Paste SMILES (one per line, optional name after a space/comma/tab):")
        lbl.setObjectName("muted")
        lay.addWidget(lbl)

        self._manual_text = QTextEdit()
        self._manual_text.setPlaceholderText(
            "PFOA\tOC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F\n"
            "PFOS\tOS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
        )
        lay.addWidget(self._manual_text, stretch=1)

        parse_btn = QPushButton("Parse Pasted Data")
        parse_btn.setObjectName("btn_secondary")
        parse_btn.clicked.connect(self._parse_manual)
        lay.addWidget(parse_btn)

        return w

    # ── Options panel ─────────────────────────────────────────────────────

    def _build_options_panel(self) -> QWidget:
        panel = QGroupBox("Classification Options")
        lay = QVBoxLayout(panel)
        lay.setSpacing(10)

        form = QFormLayout()
        form.setSpacing(8)

        self._halogen_combo = QComboBox()
        self._halogen_combo.addItems(["F (fluorine only)", "F, Cl", "F, Cl, Br", "F, Cl, Br, I"])
        form.addRow("Halogens:", self._halogen_combo)

        self._saturation_combo = QComboBox()
        self._saturation_combo.addItems(["Both (per + poly)", "Perfluorinated only", "Polyfluorinated only"])
        form.addRow("Saturation filter:", self._saturation_combo)

        lay.addLayout(form)

        self._compute_metrics_cb = QCheckBox("Compute component graph metrics")
        self._compute_metrics_cb.setChecked(True)
        lay.addWidget(self._compute_metrics_cb)

        self._check_defs_cb = QCheckBox("Check PFAS definitions (OECD, EU, OPPT…)")
        self._check_defs_cb.setChecked(True)
        lay.addWidget(self._check_defs_cb)

        lay.addWidget(self._hline())

        defs_lbl = QLabel("PFAS Definitions to check:")
        defs_lbl.setObjectName("muted")
        lay.addWidget(defs_lbl)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setMaximumHeight(180)
        defs_widget = QWidget()
        self._defs_layout = QVBoxLayout(defs_widget)
        self._defs_layout.setSpacing(4)
        scroll.setWidget(defs_widget)
        lay.addWidget(scroll)

        self._populate_definitions()

        lay.addStretch()
        return panel

    def _populate_definitions(self):
        try:
            from PFASGroups import get_PFASDefinitions
            defs = get_PFASDefinitions()
            self._def_checkboxes = {}
            for d in defs:
                cb = QCheckBox(f"{d.id}. {d.name}")
                cb.setChecked(True)
                cb.setToolTip(d.description or "")
                self._defs_layout.addWidget(cb)
                self._def_checkboxes[d.id] = cb
        except Exception:
            warn = QLabel("PFASGroups not found — PFAS definition screening unavailable.")
            warn.setObjectName("muted")
            self._defs_layout.addWidget(warn)
            self._def_checkboxes = {}

    # ── Slots ─────────────────────────────────────────────────────────────

    def _browse_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open Dataset",
            str(Path.home()),
            "Data files (*.csv *.xlsx *.xls *.db *.sqlite *.sqlite3);;All files (*)"
        )
        if not path:
            return
        self._file_path.setText(path)
        # Populate sheet picker for Excel files
        ext = Path(path).suffix.lower()
        if ext in (".xlsx", ".xls", ".xlsm"):
            from gui.utils.io_readers import get_excel_sheets
            sheets = get_excel_sheets(path)
            self._sheet_combo.blockSignals(True)
            self._sheet_combo.clear()
            for s in sheets:
                self._sheet_combo.addItem(s)
            self._sheet_combo.blockSignals(False)
            self._sheet_row_widget.setVisible(bool(sheets))
            # Load the first sheet immediately
            self._load_file(path, sheet=self._sheet_combo.currentText() or 0)
        else:
            self._sheet_row_widget.setVisible(False)
            self._load_file(path)

    def _on_sheet_changed(self, _index: int):
        path = self._file_path.text()
        if path:
            sheet = self._sheet_combo.currentText()
            self._load_file(path, sheet=sheet)

    def _load_file(self, path: str, sheet=0):
        try:
            from gui.utils.io_readers import read_file
            df = read_file(path, sheet=sheet)
            self._set_dataframe(df)
        except Exception as exc:
            QMessageBox.critical(self, "Load Error", str(exc))

    # ── DB preset helpers ─────────────────────────────────────────────────

    def _refresh_preset_combo(self):
        self._db_preset_combo.blockSignals(True)
        prev = self._db_preset_combo.currentText()
        self._db_preset_combo.clear()
        self._db_preset_combo.addItem("— select or create a connection —")
        for p in _load_presets():
            self._db_preset_combo.addItem(p.get("name", "(unnamed)"))
        # Restore selection if still present
        idx = self._db_preset_combo.findText(prev)
        if idx >= 0:
            self._db_preset_combo.setCurrentIndex(idx)
        self._db_preset_combo.blockSignals(False)
        self._on_preset_selected(self._db_preset_combo.currentIndex())

    def _on_preset_selected(self, index: int):
        if index <= 0:
            self._db_preset_data = {}
            self._db_summary.setText("No connection selected.")
            return
        presets = _load_presets()
        p = presets[index - 1]
        self._db_preset_data = p
        dialect = p.get("dialect", "")
        if dialect == "sqlite":
            self._db_summary.setText(f"SQLite · {p.get('sqlite_path', '')}")
        else:
            self._db_summary.setText(
                f"{dialect.upper()} · {p.get('username', '')}@"
                f"{p.get('host', '')}:{p.get('port', '')} / {p.get('database', '')}"
            )

    def _open_preset_editor(self):
        presets = _load_presets()
        idx = self._db_preset_combo.currentIndex()
        existing = presets[idx - 1] if idx > 0 else None
        dlg = _DbPresetDialog(existing, parent=self)
        if dlg.exec() == QDialog.DialogCode.Accepted:
            data = dlg.get_data()
            if existing is None:
                presets.append(data)
            else:
                presets[idx - 1] = data
            _save_presets(presets)
            self._refresh_preset_combo()
            # Re-select the new/edited preset
            new_idx = self._db_preset_combo.findText(data["name"])
            if new_idx >= 0:
                self._db_preset_combo.setCurrentIndex(new_idx)

    def _delete_preset(self):
        idx = self._db_preset_combo.currentIndex()
        if idx <= 0:
            return
        name = self._db_preset_combo.currentText()
        reply = QMessageBox.question(
            self, "Delete preset",
            f"Delete connection '{name}'?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
        )
        if reply == QMessageBox.StandardButton.Yes:
            presets = _load_presets()
            del presets[idx - 1]
            _save_presets(presets)
            self._refresh_preset_combo()

    def _load_from_db(self):
        p = self._db_preset_data
        if not p:
            QMessageBox.warning(self, "No connection",
                                "Please select or create a database connection preset.")
            return
        try:
            from gui.utils.io_readers import read_database
            password = p.get("password", "")
            if password.startswith("env:"):
                password = os.environ.get(password[4:].strip(), "")
            else:
                password = _simple_deobfuscate(password)

            df = read_database(
                dialect=p.get("dialect", "sqlite"),
                sqlite_path=p.get("sqlite_path", ""),
                host=p.get("host", ""),
                port=int(p.get("port", 0)),
                database=p.get("database", ""),
                username=p.get("username", ""),
                password=password,
                table=self._db_table.text().strip(),
                query=self._db_query.toPlainText().strip(),
            )
            self._set_dataframe(df)
        except Exception as exc:
            QMessageBox.critical(self, "Database Error", str(exc))

    def _parse_manual(self):
        """Parse manually pasted SMILES + optional names into a DataFrame."""
        import pandas as pd, re
        text = self._manual_text.toPlainText().strip()
        if not text:
            return
        rows = []
        for line in text.splitlines():
            line = line.strip()
            if not line:
                continue
            parts = re.split(r"[\t,;\s]+", line, maxsplit=1)
            if len(parts) == 2:
                # Decide which part looks like SMILES
                if any(c in parts[0] for c in "()[]=#"):
                    rows.append({"smiles": parts[0], "name": parts[1]})
                else:
                    rows.append({"name": parts[0], "smiles": parts[1]})
            else:
                rows.append({"smiles": parts[0], "name": ""})
        df = pd.DataFrame(rows)
        self._set_dataframe(df)

    def _set_dataframe(self, df):
        import pandas as pd
        from gui.utils.io_readers import get_smiles_column, get_name_column
        self._df = df
        # Update column selectors
        cols = list(df.columns)
        for combo in (self._smiles_col, self._name_col):
            combo.clear()
        self._name_col.addItem("(none)")
        for c in cols:
            self._smiles_col.addItem(c)
            self._name_col.addItem(c)

        smi_guess = get_smiles_column(df)
        if smi_guess:
            self._smiles_col.setCurrentText(smi_guess)
        name_guess = get_name_column(df)
        if name_guess:
            self._name_col.setCurrentText(name_guess)

        self._preview_lbl.setText(
            f"Loaded {len(df):,} rows × {len(cols)} columns."
        )

    def _run(self):
        if self._df is None:
            QMessageBox.warning(self, "No Data", "Please load a dataset first.")
            return

        smiles_col = self._smiles_col.currentText()
        if smiles_col not in self._df.columns:
            QMessageBox.warning(self, "Column Error",
                                f"SMILES column '{smiles_col}' not found in dataset.")
            return

        smiles_list = self._df[smiles_col].fillna("").tolist()
        name_col = self._name_col.currentText()
        names = self._df[name_col].fillna("").tolist() if name_col != "(none)" else []

        halogen_map = {0: "F", 1: ["F", "Cl"], 2: ["F", "Cl", "Br"], 3: ["F", "Cl", "Br", "I"]}
        halogens = halogen_map[self._halogen_combo.currentIndex()]

        self._run_btn.setEnabled(False)
        self._progress.setValue(0)
        self._status_lbl.setText("Starting…")

        self._worker = ClassifyWorker(
            smiles=smiles_list,
            names=names,
            halogens=halogens,
            compute_metrics=self._compute_metrics_cb.isChecked(),
            check_definitions=self._check_defs_cb.isChecked(),
        )
        self._worker.progress.connect(self._progress.setValue)
        self._worker.status.connect(self._status_lbl.setText)
        self._worker.result.connect(self._on_result)
        self._worker.error.connect(self._on_error)
        self._worker.finished.connect(lambda: self._run_btn.setEnabled(True))
        self._worker.start()

    def _on_result(self, embedding_set):
        self._last_embedding = embedding_set
        n = len(embedding_set)
        self._status_lbl.setText(f"Classified {n:,} molecule(s).")
        self.classification_done.emit(embedding_set)

    def _on_error(self, msg: str):
        self._status_lbl.setText("Error.")
        self._progress.setValue(0)
        QMessageBox.critical(self, "Classification Error", msg)

    # ── Helpers ───────────────────────────────────────────────────────────

    @staticmethod
    def _hline() -> QFrame:
        line = QFrame()
        line.setFrameShape(QFrame.Shape.HLine)
        line.setFrameShadow(QFrame.Shadow.Sunken)
        return line


# ── DB Preset Editor Dialog ────────────────────────────────────────────────────

class _DbPresetDialog(QDialog):
    """Pop-up for creating or editing a saved database connection preset."""

    def __init__(self, existing: dict | None = None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Database Connection" if existing is None else "Edit Connection")
        self.setMinimumWidth(460)
        self._existing = existing or {}
        self._build_ui()
        if existing:
            self._populate(existing)

    def _build_ui(self):
        lay = QVBoxLayout(self)
        lay.setSpacing(10)

        form = QFormLayout()
        form.setSpacing(8)

        self._name_edit = QLineEdit()
        self._name_edit.setPlaceholderText("e.g. My local MariaDB")
        form.addRow("Display name:", self._name_edit)

        self._dialect_combo = QComboBox()
        self._dialect_combo.addItems(["SQLite", "MariaDB / MySQL", "PostgreSQL"])
        self._dialect_combo.currentIndexChanged.connect(self._update_fields)
        form.addRow("Database type:", self._dialect_combo)

        # SQLite path row
        self._sqlite_path_edit = QLineEdit()
        self._sqlite_path_edit.setPlaceholderText("Path to .db / .sqlite file")
        sqlite_browse = QPushButton("Browse…")
        sqlite_browse.clicked.connect(self._browse_sqlite)
        sqlite_row = QHBoxLayout()
        sqlite_row.addWidget(self._sqlite_path_edit, stretch=1)
        sqlite_row.addWidget(sqlite_browse)
        sqlite_container = QWidget()
        sqlite_container.setLayout(sqlite_row)
        self._sqlite_path_label = QLabel("File path:")
        form.addRow(self._sqlite_path_label, sqlite_container)

        # Remote connection fields
        self._host_edit = QLineEdit("localhost")
        self._host_label = QLabel("Host:")
        form.addRow(self._host_label, self._host_edit)

        self._port_spin = QSpinBox()
        self._port_spin.setRange(1, 65535)
        self._port_spin.setValue(3306)
        self._port_label = QLabel("Port:")
        form.addRow(self._port_label, self._port_spin)

        self._db_name_edit = QLineEdit()
        self._db_name_edit.setPlaceholderText("database name")
        self._db_name_label = QLabel("Database:")
        form.addRow(self._db_name_label, self._db_name_edit)

        self._user_edit = QLineEdit()
        self._user_label = QLabel("Username:")
        form.addRow(self._user_label, self._user_edit)

        # Password row with env-var toggle
        self._pw_type_grp = QButtonGroup(self)
        self._pw_literal_rb = QRadioButton("Enter password:")
        self._pw_env_rb = QRadioButton("From environment variable:")
        self._pw_type_grp.addButton(self._pw_literal_rb)
        self._pw_type_grp.addButton(self._pw_env_rb)
        self._pw_literal_rb.setChecked(True)
        self._pw_literal_rb.toggled.connect(self._update_pw_fields)

        pw_type_row = QHBoxLayout()
        pw_type_row.addWidget(self._pw_literal_rb)
        pw_type_row.addWidget(self._pw_env_rb)
        pw_type_row.addStretch()
        pw_type_widget = QWidget()
        pw_type_widget.setLayout(pw_type_row)

        self._pw_password_edit = QLineEdit()
        self._pw_password_edit.setEchoMode(QLineEdit.EchoMode.Password)
        self._pw_password_edit.setPlaceholderText("password (stored obfuscated)")

        self._pw_env_edit = QLineEdit()
        self._pw_env_edit.setPlaceholderText("e.g.  DB_PASSWORD")
        self._pw_env_edit.setVisible(False)

        self._pw_label = QLabel("Password:")
        form.addRow(self._pw_label, pw_type_widget)
        form.addRow("", self._pw_password_edit)
        form.addRow("", self._pw_env_edit)

        lay.addLayout(form)

        # Buttons
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self._validate_and_accept)
        buttons.rejected.connect(self.reject)
        lay.addWidget(buttons)

        self._update_fields(self._dialect_combo.currentIndex())

    def _update_fields(self, index: int):
        is_sqlite = (index == 0)
        sqlite_widgets = [self._sqlite_path_label]
        remote_widgets = [self._host_label, self._host_edit,
                          self._port_label, self._port_spin,
                          self._db_name_label, self._db_name_edit,
                          self._user_label, self._user_edit,
                          self._pw_label, self._pw_password_edit, self._pw_env_edit]
        sqlite_containers = []
        # Walk the form to find the SQLite row widget
        form = self.layout().itemAt(0).layout()
        for i in range(form.rowCount()):
            lbl = form.itemAt(i, QFormLayout.ItemRole.LabelRole)
            fld = form.itemAt(i, QFormLayout.ItemRole.FieldRole)
            if lbl and lbl.widget() == self._sqlite_path_label:
                if fld and fld.widget():
                    fld.widget().setVisible(is_sqlite)
            if lbl and lbl.widget() in (self._host_label, self._port_label,
                                         self._db_name_label, self._user_label):
                v = lbl.widget()
                f = form.itemAt(i, QFormLayout.ItemRole.FieldRole)
                v.setVisible(not is_sqlite)
                if f and f.widget():
                    f.widget().setVisible(not is_sqlite)

        self._sqlite_path_label.setVisible(is_sqlite)
        self._host_label.setVisible(not is_sqlite)
        self._host_edit.setVisible(not is_sqlite)
        self._port_label.setVisible(not is_sqlite)
        self._port_spin.setVisible(not is_sqlite)
        self._db_name_label.setVisible(not is_sqlite)
        self._db_name_edit.setVisible(not is_sqlite)
        self._user_label.setVisible(not is_sqlite)
        self._user_edit.setVisible(not is_sqlite)
        self._pw_label.setVisible(not is_sqlite)
        self._pw_password_edit.setVisible(not is_sqlite and self._pw_literal_rb.isChecked())
        self._pw_env_edit.setVisible(not is_sqlite and self._pw_env_rb.isChecked())

        if index == 2:  # PostgreSQL
            self._port_spin.setValue(5432)
        elif index == 1:  # MariaDB
            self._port_spin.setValue(3306)

    def _update_pw_fields(self):
        is_literal = self._pw_literal_rb.isChecked()
        is_sqlite = self._dialect_combo.currentIndex() == 0
        self._pw_password_edit.setVisible(not is_sqlite and is_literal)
        self._pw_env_edit.setVisible(not is_sqlite and not is_literal)

    def _browse_sqlite(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open SQLite Database", str(Path.home()),
            "SQLite (*.db *.sqlite *.sqlite3);;All files (*)"
        )
        if path:
            self._sqlite_path_edit.setText(path)

    def _populate(self, data: dict):
        self._name_edit.setText(data.get("name", ""))
        dialect = data.get("dialect", "sqlite")
        idx = {"sqlite": 0, "mariadb": 1, "mysql": 1, "postgresql": 2}.get(dialect.lower(), 0)
        self._dialect_combo.setCurrentIndex(idx)
        self._sqlite_path_edit.setText(data.get("sqlite_path", ""))
        self._host_edit.setText(data.get("host", "localhost"))
        self._port_spin.setValue(int(data.get("port", 3306)))
        self._db_name_edit.setText(data.get("database", ""))
        self._user_edit.setText(data.get("username", ""))
        raw_pw = data.get("password", "")
        if raw_pw.startswith("env:"):
            self._pw_env_rb.setChecked(True)
            self._pw_env_edit.setText(raw_pw[4:].strip())
        else:
            self._pw_literal_rb.setChecked(True)
            # Don't show the obfuscated form — just leave blank for re-entry
            self._pw_password_edit.setPlaceholderText("(leave blank to keep existing)")

    def _validate_and_accept(self):
        if not self._name_edit.text().strip():
            QMessageBox.warning(self, "Missing name", "Please enter a display name.")
            return
        self.accept()

    def get_data(self) -> dict:
        dialect_map = {0: "sqlite", 1: "mariadb", 2: "postgresql"}
        dialect = dialect_map[self._dialect_combo.currentIndex()]
        data: dict = {
            "name": self._name_edit.text().strip(),
            "dialect": dialect,
        }
        if dialect == "sqlite":
            data["sqlite_path"] = self._sqlite_path_edit.text().strip()
        else:
            data["host"] = self._host_edit.text().strip()
            data["port"] = self._port_spin.value()
            data["database"] = self._db_name_edit.text().strip()
            data["username"] = self._user_edit.text().strip()
            if self._pw_env_rb.isChecked():
                data["password"] = "env:" + self._pw_env_edit.text().strip()
            else:
                pw = self._pw_password_edit.text()
                if pw:
                    data["password"] = _simple_obfuscate(pw)
                else:
                    # Keep existing obfuscated password if not re-entered
                    data["password"] = self._existing.get("password", "")
        return data
