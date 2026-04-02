"""
PFASGroups GUI stylesheet.

All colours are taken from PFASGroups/data/color_scheme.yaml:
  Primary:  orange-red #E15D0B  |  blue #306DBA  |  magenta #9D206C  |  dark-purple #51127C
  Tints (lightest → darkest listed in YAML):
    orange-red: #E15D0B #E67834 #EB935C #F0AE85
    blue:        #306DBA #5385C6 #759ED1 #97B6DD
    magenta:     #9D206C #AD4585 #BE6A9D #CE8FB5
    dark-purple: #51127C #6E3A92 #8B61A8 #A888BD
"""

# ---------------------------------------------------------------------------
# Palette constants (re-used by tabs for matplotlib figures, Plotly, etc.)
# ---------------------------------------------------------------------------

C_ORANGE        = "#E15D0B"
C_BLUE          = "#306DBA"
C_MAGENTA       = "#9D206C"
C_DARK_PURPLE   = "#51127C"

C_ORANGE_LIGHT  = "#F0AE85"
C_BLUE_LIGHT    = "#97B6DD"
C_MAGENTA_LIGHT = "#CE8FB5"
C_PURPLE_LIGHT  = "#A888BD"

C_ORANGE_MID    = "#EB935C"
C_BLUE_MID      = "#759ED1"
C_MAGENTA_MID   = "#BE6A9D"
C_PURPLE_MID    = "#8B61A8"

C_WHITE         = "#FFFFFF"
C_BG            = "#F8F8F9"          # very light grey page background
C_SURFACE       = "#FFFFFF"          # card / panel background
C_BORDER        = "#D4D4D8"          # subtle border
C_TEXT          = "#1A1A2E"          # near-black body text
C_TEXT_MUTED    = "#6B7280"          # secondary text
C_SUCCESS       = "#16A34A"
C_ERROR         = "#DC2626"
C_WARNING       = "#D97706"

PALETTE = [C_ORANGE, C_BLUE, C_MAGENTA, C_DARK_PURPLE,
           C_ORANGE_MID, C_BLUE_MID, C_MAGENTA_MID, C_PURPLE_MID]

# ---------------------------------------------------------------------------
# QSS stylesheet
# ---------------------------------------------------------------------------

QSS = f"""
/* ── Global ─────────────────────────────────────────────────────────────── */
QWidget {{
    background-color: {C_BG};
    color: {C_TEXT};
    font-family: "Segoe UI", Arial, sans-serif;
    font-size: 13px;
}}

/* ── Main window ─────────────────────────────────────────────────────────── */
QMainWindow {{
    background-color: {C_BG};
}}

/* ── Top header bar ──────────────────────────────────────────────────────── */
#header_widget {{
    background-color: {C_DARK_PURPLE};
    border-bottom: 3px solid {C_ORANGE};
}}
#header_title {{
    color: {C_WHITE};
    font-size: 20px;
    font-weight: 700;
    letter-spacing: 0.5px;
}}
#header_subtitle {{
    color: {C_PURPLE_LIGHT};
    font-size: 12px;
}}

/* ── Tab bar ─────────────────────────────────────────────────────────────── */
QTabWidget::pane {{
    border: 1px solid {C_BORDER};
    background-color: {C_SURFACE};
    border-top: none;
}}
QTabWidget::tab-bar {{
    left: 0;
}}
QTabBar::tab {{
    background-color: {C_BG};
    color: {C_TEXT_MUTED};
    border: 1px solid {C_BORDER};
    border-bottom: none;
    padding: 8px 18px;
    margin-right: 2px;
    border-top-left-radius: 5px;
    border-top-right-radius: 5px;
    min-width: 110px;
    font-size: 12px;
    font-weight: 500;
}}
QTabBar::tab:selected {{
    background-color: {C_SURFACE};
    color: {C_DARK_PURPLE};
    border-bottom: 3px solid {C_ORANGE};
    font-weight: 700;
}}
QTabBar::tab:hover:!selected {{
    background-color: {C_PURPLE_LIGHT};
    color: {C_DARK_PURPLE};
}}

/* ── Group boxes ─────────────────────────────────────────────────────────── */
QGroupBox {{
    border: 1px solid {C_BORDER};
    border-radius: 6px;
    margin-top: 10px;
    padding: 8px;
    font-weight: 600;
    font-size: 12px;
    color: {C_DARK_PURPLE};
    background-color: {C_SURFACE};
}}
QGroupBox::title {{
    subcontrol-origin: margin;
    subcontrol-position: top left;
    padding: 0 6px;
    left: 10px;
    color: {C_DARK_PURPLE};
}}

/* ── Buttons ─────────────────────────────────────────────────────────────── */
QPushButton {{
    background-color: {C_BLUE};
    color: {C_WHITE};
    border: none;
    border-radius: 5px;
    padding: 7px 18px;
    font-weight: 600;
    font-size: 12px;
    min-height: 28px;
}}
QPushButton:hover {{
    background-color: {C_BLUE_MID};
}}
QPushButton:pressed {{
    background-color: {C_DARK_PURPLE};
}}
QPushButton:disabled {{
    background-color: {C_BORDER};
    color: {C_TEXT_MUTED};
}}
QPushButton#btn_primary {{
    background-color: {C_ORANGE};
    font-size: 13px;
    padding: 8px 24px;
}}
QPushButton#btn_primary:hover {{
    background-color: {C_ORANGE_MID};
}}
QPushButton#btn_danger {{
    background-color: {C_ERROR};
}}
QPushButton#btn_success {{
    background-color: {C_SUCCESS};
}}
QPushButton#btn_secondary {{
    background-color: {C_SURFACE};
    color: {C_BLUE};
    border: 1.5px solid {C_BLUE};
}}
QPushButton#btn_secondary:hover {{
    background-color: {C_BLUE_LIGHT};
    color: {C_DARK_PURPLE};
}}

/* ── Line edits, text edits, combo boxes ─────────────────────────────────── */
QLineEdit, QTextEdit, QPlainTextEdit {{
    background-color: {C_WHITE};
    border: 1px solid {C_BORDER};
    border-radius: 4px;
    padding: 5px 8px;
    selection-background-color: {C_BLUE_LIGHT};
    color: {C_TEXT};
}}
QLineEdit:focus, QTextEdit:focus, QPlainTextEdit:focus {{
    border: 1.5px solid {C_BLUE};
}}
QComboBox {{
    background-color: {C_WHITE};
    border: 1px solid {C_BORDER};
    border-radius: 4px;
    padding: 5px 8px;
    min-height: 28px;
    color: {C_TEXT};
}}
QComboBox:focus {{
    border: 1.5px solid {C_BLUE};
}}
QComboBox::drop-down {{
    border: none;
    width: 22px;
}}
QComboBox QAbstractItemView {{
    background-color: {C_WHITE};
    border: 1px solid {C_BORDER};
    selection-background-color: {C_BLUE_LIGHT};
}}

/* ── Spin boxes ──────────────────────────────────────────────────────────── */
QDoubleSpinBox, QSpinBox {{
    background-color: {C_WHITE};
    border: 1px solid {C_BORDER};
    border-radius: 4px;
    padding: 4px 8px;
    min-height: 26px;
    color: {C_TEXT};
}}
QDoubleSpinBox:focus, QSpinBox:focus {{
    border: 1.5px solid {C_BLUE};
}}

/* ── Check boxes & radio buttons ─────────────────────────────────────────── */
QCheckBox {{
    spacing: 8px;
    color: {C_TEXT};
}}
QCheckBox::indicator {{
    width: 16px;
    height: 16px;
    border: 1.5px solid {C_BLUE};
    border-radius: 3px;
    background-color: {C_WHITE};
}}
QCheckBox::indicator:checked {{
    background-color: {C_BLUE};
}}
QRadioButton {{
    spacing: 8px;
    color: {C_TEXT};
}}
QRadioButton::indicator {{
    width: 16px;
    height: 16px;
    border-radius: 8px;
    border: 1.5px solid {C_BLUE};
    background-color: {C_WHITE};
}}
QRadioButton::indicator:checked {{
    background-color: {C_BLUE};
}}

/* ── Progress bar ────────────────────────────────────────────────────────── */
QProgressBar {{
    border: 1px solid {C_BORDER};
    border-radius: 4px;
    background-color: {C_BG};
    text-align: center;
    height: 18px;
    font-size: 11px;
    color: {C_TEXT};
}}
QProgressBar::chunk {{
    background-color: {C_ORANGE};
    border-radius: 3px;
}}

/* ── Table views ─────────────────────────────────────────────────────────── */
QTableWidget, QTableView {{
    background-color: {C_WHITE};
    alternate-background-color: {C_BG};
    gridline-color: {C_BORDER};
    border: 1px solid {C_BORDER};
    border-radius: 4px;
    selection-background-color: {C_BLUE_LIGHT};
    selection-color: {C_TEXT};
}}
QHeaderView::section {{
    background-color: {C_DARK_PURPLE};
    color: {C_WHITE};
    padding: 6px 10px;
    font-weight: 700;
    border: 1px solid {C_BORDER};
    font-size: 12px;
}}
QHeaderView::section:hover {{
    background-color: {C_PURPLE_MID};
}}

/* ── Scroll bars ─────────────────────────────────────────────────────────── */
QScrollBar:vertical {{
    background: {C_BG};
    width: 10px;
    margin: 0;
}}
QScrollBar::handle:vertical {{
    background: {C_PURPLE_MID};
    border-radius: 5px;
    min-height: 20px;
}}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{
    height: 0;
}}
QScrollBar:horizontal {{
    background: {C_BG};
    height: 10px;
    margin: 0;
}}
QScrollBar::handle:horizontal {{
    background: {C_PURPLE_MID};
    border-radius: 5px;
    min-width: 20px;
}}
QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {{
    width: 0;
}}

/* ── Labels ──────────────────────────────────────────────────────────────── */
QLabel {{
    background-color: transparent;
}}
QLabel#section_title {{
    font-size: 14px;
    font-weight: 700;
    color: {C_DARK_PURPLE};
    padding-bottom: 4px;
}}
QLabel#muted {{
    color: {C_TEXT_MUTED};
    font-size: 11px;
}}
QLabel#badge_pass {{
    background-color: {C_SUCCESS};
    color: {C_WHITE};
    border-radius: 4px;
    padding: 2px 8px;
    font-weight: 700;
    font-size: 11px;
}}
QLabel#badge_fail {{
    background-color: {C_ERROR};
    color: {C_WHITE};
    border-radius: 4px;
    padding: 2px 8px;
    font-weight: 700;
    font-size: 11px;
}}
QLabel#badge_warn {{
    background-color: {C_WARNING};
    color: {C_WHITE};
    border-radius: 4px;
    padding: 2px 8px;
    font-weight: 700;
    font-size: 11px;
}}

/* ── Splitter ─────────────────────────────────────────────────────────────── */
QSplitter::handle {{
    background-color: {C_BORDER};
}}

/* ── Tool tips ────────────────────────────────────────────────────────────── */
QToolTip {{
    background-color: {C_DARK_PURPLE};
    color: {C_WHITE};
    border: 1px solid {C_PURPLE_MID};
    padding: 4px 8px;
    font-size: 11px;
    border-radius: 4px;
}}

/* ── Status bar ──────────────────────────────────────────────────────────── */
QStatusBar {{
    background-color: {C_DARK_PURPLE};
    color: {C_PURPLE_LIGHT};
    font-size: 11px;
    border-top: 1px solid {C_ORANGE};
}}
QStatusBar::item {{
    border: none;
}}

/* ── List widget ─────────────────────────────────────────────────────────── */
QListWidget {{
    background-color: {C_WHITE};
    border: 1px solid {C_BORDER};
    border-radius: 4px;
    alternate-background-color: {C_BG};
}}
QListWidget::item:selected {{
    background-color: {C_BLUE_LIGHT};
    color: {C_TEXT};
}}
QListWidget::item:hover {{
    background-color: {C_BG};
}}

/* ── Frame dividers ──────────────────────────────────────────────────────── */
QFrame[frameShape="4"],
QFrame[frameShape="5"] {{
    color: {C_BORDER};
}}
"""


def apply(app):
    """Apply the PFASGroups stylesheet to a QApplication."""
    app.setStyleSheet(QSS)
