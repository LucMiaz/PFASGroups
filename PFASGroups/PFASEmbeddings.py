from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple, Union
import os
import re
from pathlib import Path

if TYPE_CHECKING:
    try:
        import sqlalchemy
    except ImportError:
        pass

from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image
from io import BytesIO
import pandas as pd
import numpy as np
from .parser import load_HalogenGroups


def _load_palette() -> List[str]:
    """Load hex colours from color_scheme.yaml (stdlib only, no pyyaml needed)."""
    _defaults = ["#E15D0B", "#306DBA", "#9D206C", "#51127C"]
    try:
        _p = Path(__file__).parent / "data" / "color_scheme.yaml"
        _colors = re.findall(r'"(#[0-9A-Fa-f]{6})"', _p.read_text())
        if len(_colors) >= 4:
            return _colors[:4]
    except Exception:
        pass
    return _defaults


_PALETTE = _load_palette()
# C0=orange, C1=blue (FG table), C2=magenta (metrics table), C3=dark-purple
_C0, _C1, _C2, _C3 = _PALETTE


def _hex_to_rgb_float(h: str) -> Tuple[float, float, float]:
    """Convert a '#RRGGBB' hex string to an RGB float triple in [0, 1]."""
    h = h.lstrip('#')
    return tuple(int(h[i:i+2], 16) / 255.0 for i in (0, 2, 4))  # type: ignore


def _lighter(h: str, factor: float = 0.82) -> str:
    """Return a lighter hex colour by blending *h* with white."""
    r, g, b = _hex_to_rgb_float(h)
    r2 = int((r + (1 - r) * factor) * 255)
    g2 = int((g + (1 - g) * factor) * 255)
    b2 = int((b + (1 - b) * factor) * 255)
    return f'#{r2:02X}{g2:02X}{b2:02X}'

# ---------------------------------------------------------------------------
# Sentinel for "argument not supplied" (distinct from None)
# ---------------------------------------------------------------------------
_UNSET = object()

# ---------------------------------------------------------------------------
# Embedding helpers (used by PFASEmbedding.to_array / PFASEmbeddingSet.to_array)
# ---------------------------------------------------------------------------

def _finite(v) -> bool:
    """True when v is a real, finite number (not NaN/Inf)."""
    try:
        return v == v and abs(v) != float('inf')
    except TypeError:
        return False


def _agg(comps: list, metric: str, aggregation: str) -> float:
    """Aggregate *metric* over component dicts using mean or median."""
    vals = [c[metric] for c in comps if c.get(metric) is not None and _finite(c[metric])]
    if not vals:
        return 0.0
    if aggregation == 'median':
        return float(np.median(vals))
    return float(np.mean(vals))


def _encode_count(match, mode: str) -> float:
    """Scalar encoding for one matched group given count *mode*."""
    if match is None or not match.get('components'):
        return 0.0
    comps = match['components']
    if mode == 'binary':
        return 1.0
    if mode == 'count':
        return float(len(comps))
    if mode == 'max_component':
        return float(max(c.get('size', 0) for c in comps))
    if mode == 'total_component':
        return float(sum(c.get('size', 0) for c in comps))
    return 0.0


def _mol_metric(all_comps: list, metric: str) -> float:
    """Molecule-wide scalar for *metric* over all matched components."""
    if metric == 'n_components':
        return float(len(all_comps))
    if metric == 'total_size':
        return float(sum(c.get('size', 0) or 0 for c in all_comps))
    if metric == 'mean_size':
        sizes = [c.get('size', 0) or 0 for c in all_comps]
        return float(np.mean(sizes)) if sizes else 0.0
    if metric == 'max_size':
        sizes = [c.get('size', 0) or 0 for c in all_comps]
        return float(max(sizes)) if sizes else 0.0
    if metric == 'mean_branching':
        return _agg(all_comps, 'branching', 'mean')
    if metric == 'max_branching':
        vals = [c.get('branching') for c in all_comps if c.get('branching') is not None]
        return float(max(vals)) if vals else 0.0
    if metric == 'mean_eccentricity':
        return _agg(all_comps, 'mean_eccentricity', 'mean')
    if metric == 'max_diameter':
        vals = [c.get('diameter') for c in all_comps if c.get('diameter') is not None]
        return float(max(vals)) if vals else 0.0
    if metric == 'mean_component_fraction':
        return _agg(all_comps, 'component_fraction', 'mean')
    if metric == 'max_component_fraction':
        vals = [c.get('component_fraction') for c in all_comps if c.get('component_fraction') is not None]
        return float(max(vals)) if vals else 0.0
    return 0.0


def _select_groups(group_selection, selected_group_ids, pfas_groups):
    """Translate group_selection / selected_group_ids to a list of group objects."""
    id_to_group = {g.id: g for g in pfas_groups}

    if selected_group_ids is not None:
        return [id_to_group[gid] for gid in selected_group_ids if gid in id_to_group]

    if group_selection is None or group_selection == 'all':
        return list(pfas_groups)
    if group_selection == 'oecd':
        return [id_to_group[gid] for gid in range(1, 29) if gid in id_to_group]
    if group_selection == 'generic':
        return [id_to_group[gid] for gid in range(29, 56) if gid in id_to_group]
    if group_selection == 'telomers':
        return [id_to_group[gid] for gid in range(74, 117) if gid in id_to_group]
    if group_selection == 'generic+telomers':
        ids = list(range(29, 56)) + list(range(74, 117))
        return [id_to_group[gid] for gid in ids if gid in id_to_group]

    from .getter import get_HalogenGroups
    raw = get_HalogenGroups()
    compute_raw = [g for g in raw if g.get('compute', True)]
    matching_ids = {
        g['id'] for g in compute_raw
        if g.get('test', {}).get('category', 'other') == group_selection
    }
    if matching_ids:
        return [id_to_group[gid] for gid in matching_ids if gid in id_to_group]

    raise ValueError(
        f"Unknown group_selection: {group_selection!r}. "
        f"Choose from: 'all', 'oecd', 'generic', 'telomers', 'generic+telomers'"
    )

Color = Tuple[float, float, float]

# ---------------------------------------------------------------------------
# ANSI terminal styling helpers
# ---------------------------------------------------------------------------

_ANSI_RESET  = "\033[0m"
_ANSI_BOLD   = "\033[1m"

_ANSI_HALOGEN: Dict[str, str] = {
    "F":  "\033[96m",   # bright cyan
    "Cl": "\033[92m",   # bright green
    "Br": "\033[93m",   # bright yellow
    "I":  "\033[95m",   # bright magenta
}
_ANSI_FORM: Dict[str, str] = {
    "alkyl":  "\033[34m",  # blue
    "cyclic": "\033[35m",  # magenta
}
_ANSI_SAT: Dict[str, str] = {
    "per":  "\033[31m",  # red
    "poly": "\033[33m",  # yellow
}

def _ansi(text: str, *codes: str) -> str:
    """Wrap *text* with ANSI escape codes and reset afterwards."""
    return "".join(codes) + text + _ANSI_RESET


# ---------------------------------------------------------------------------
# Molecule-highlight colour palettes (RGB float triples, 0–1)
# ---------------------------------------------------------------------------

# Highlight colours by halogen element — mapped to the project colour palette
_HALOGEN_COLORS: Dict[str, Color] = {
    "F":  _hex_to_rgb_float(_C1),   # blue  (most common PFAS element)
    "Cl": _hex_to_rgb_float(_C0),   # orange
    "Br": _hex_to_rgb_float(_C2),   # magenta
    "I":  _hex_to_rgb_float(_C3),   # dark purple
}
_HALOGEN_COLOR_DEFAULT: Color = (0.75, 0.75, 0.75)  # grey for unknown

# Tint modifiers for form (shift hue slightly)
_FORM_TINT: Dict[str, Tuple[float, float, float]] = {
    "alkyl":  (0.0,  0.0,  0.0),   # no change
    "cyclic": (0.12, -0.08, 0.05), # warm shift
}

# Brightness modifiers for saturation
_SAT_BRIGHTNESS: Dict[str, float] = {
    "per":  1.0,   # full saturation → bright
    "poly": 0.65,  # partial saturation → dimmer
}

def _component_color(halogen: Optional[str], form: Optional[str], saturation: Optional[str]) -> Color:
    """Return an RGB highlight colour encoding halogen, form and saturation."""
    base = _HALOGEN_COLORS.get(halogen or "", _HALOGEN_COLOR_DEFAULT)
    tint = _FORM_TINT.get(form or "", (0.0, 0.0, 0.0))
    brightness = _SAT_BRIGHTNESS.get(saturation or "", 0.85)
    r = min(1.0, max(0.0, (base[0] + tint[0]) * brightness))
    g = min(1.0, max(0.0, (base[1] + tint[1]) * brightness))
    b = min(1.0, max(0.0, (base[2] + tint[2]) * brightness))
    return (r, g, b)


# Simple color palette to distinguish PFAS groups in highlight plots (legacy fallback)
_GROUP_COLORS: List[Color] = [
    _hex_to_rgb_float(_C0),  # orange
    _hex_to_rgb_float(_C1),  # blue
    _hex_to_rgb_float(_C2),  # magenta
    _hex_to_rgb_float(_C3),  # dark purple
    (0.00, 0.70, 0.70),      # teal (extra)
    (0.60, 0.60, 0.60),      # grey (extra)
]

# ---------------------------------------------------------------------------
# Component-SMARTS metadata cache
# ---------------------------------------------------------------------------

_COMPONENT_META_CACHE: Optional[Dict[str, Dict[str, Optional[str]]]] = None


def _get_component_meta() -> Dict[str, Dict[str, Optional[str]]]:
    """Return a dict mapping component SMARTS name → {halogen, form, saturation}.

    Built from the same preprocessed component dictionary used by the parser,
    so names are guaranteed to match the SMARTS labels stored in match results.
    """
    global _COMPONENT_META_CACHE
    if _COMPONENT_META_CACHE is not None:
        return _COMPONENT_META_CACHE
    try:
        from .core import get_componentSMARTSs
        raw = get_componentSMARTSs()
        _COMPONENT_META_CACHE = {
            name: {
                "halogen":    info.get("halogen"),
                "form":       info.get("form"),
                "saturation": info.get("saturation"),
            }
            for name, info in raw.items()
        }
    except Exception:
        _COMPONENT_META_CACHE = {}
    return _COMPONENT_META_CACHE


# ---------------------------------------------------------------------------
# Group-info cache and classification helpers
# ---------------------------------------------------------------------------

_GROUP_INFO_CACHE: Optional[Dict[int, Dict[str, str]]] = None


def _get_group_info() -> Dict[int, Dict[str, str]]:
    """Return a dict mapping group_id → {name, category}.

    Category is one of ``'OECD'``, ``'generic'``, ``'telomer'``, or ``'other'``.
    """
    global _GROUP_INFO_CACHE
    if _GROUP_INFO_CACHE is not None:
        return _GROUP_INFO_CACHE
    try:
        from .getter import get_HalogenGroups
        raw = get_HalogenGroups()
        _GROUP_INFO_CACHE = {
            g["id"]: {
                "name":     g.get("name", ""),
                "category": g.get("test", {}).get("category", "other"),
            }
            for g in raw
        }
    except Exception:
        _GROUP_INFO_CACHE = {}
    return _GROUP_INFO_CACHE


# Groups always excluded from the non-OECD category label (perhalogenated /
# polyhalogenated alkyl catch-alls that add no structural specificity).
_CLASSIFY_EXCLUDED_IDS: frozenset = frozenset({51, 52})

# Name subsumption for non-OECD classification: when a more-specific group
# name (key) is present, every name in its list is suppressed.
# Keys and values must match exactly the ``group_name`` values stored in
# match results (i.e. the ``name`` field from the raw group JSON).
_SUBSUMES: Dict[str, List[str]] = {
    "sulfonamide":         ["amine"],
    "amide":               ["amine"],
    "Telomer sulfonamide": ["amine"],
    "phosphonamide":       ["amine"],
}


def _grid_images(imgs: Sequence[Image.Image], buffer: int = 4, ncols: int = 3) -> Tuple[Image.Image, int, int]:
    """Arrange PIL images in a simple grid layout.

    Parameters
    ----------
    imgs : sequence of PIL.Image
        Images to arrange.
    buffer : int, default 4
        Spacing between images in pixels.
    ncols : int, default 3
        Number of columns.
    """
    if not imgs:
        raise ValueError("No images provided to _grid_images")

    # Compute per-row layout
    rows: List[List[Image.Image]] = [list(imgs[i : i + ncols]) for i in range(0, len(imgs), ncols)]
    max_width = 0
    total_height = 0
    row_heights: List[int] = []

    for row in rows:
        row_width = sum(im.width for im in row) + buffer * (len(row) - 1 if len(row) > 0 else 0)
        max_width = max(max_width, row_width)
        h = max((im.height for im in row), default=0)
        row_heights.append(h)
        total_height += h + buffer

    if rows:
        total_height -= buffer  # no buffer after last row

    canvas = Image.new("RGBA", (max_width, total_height), (255, 255, 255, 0))

    y = 0
    for row, h in zip(rows, row_heights):
        x = 0
        for im in row:
            canvas.paste(im, (x, y))
            x += im.width + buffer
        y += h + buffer

    return canvas, max_width, total_height


def _mol_image_with_table(
    mol_img: Image.Image,
    entries: List[Tuple[str, str, str, str, str]],
    comp_metrics: Optional[Tuple[str, ...]] = None,
    mol_label: str = "",
    halogen_label: str = "",
    font_size: int = 9,
) -> Image.Image:
    """Composite a molecule image (PIL) with two formatted matplotlib tables.

    Parameters
    ----------
    mol_img : PIL Image
        The molecule drawing produced by RDKit (no legend text).
    entries : list of tuples
        Each tuple: (group_name, smarts_label, dist_center, dist_periphery)
        FG-specific metrics per matched component row.
    comp_metrics : tuple of str or None
        Component-wide row: (size, branching, eccentricity, diameter, radius,
        eff_graph_resistance, bde_eff_graph_resistance, chain_pct).
        When provided, rendered in a separate second table below the FG table.
    mol_label : str
        Optional header shown above the table (e.g. "mol#1").
    halogen_label : str
        Optional halogen element symbol (e.g. "F") shown as a badge in the
        top-right corner of the molecule image.
    font_size : int
        Font size for table body text.

    Returns
    -------
    PIL Image
        Combined molecule + tables image.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import io as _io

    dpi = 96
    mol_w, mol_h = mol_img.size
    fig_w_in = mol_w / dpi

    # Row height in inches
    row_h_in = (font_size + 5) / 72.0
    header_h_in = row_h_in * 1.4
    n_rows = max(1, len(entries))
    label_h_in = (font_size + 4) / 72.0 * 1.3 if mol_label else 0.0
    # Second table height: header + 1 data row (if comp_metrics provided)
    metrics_tbl_h_in = (header_h_in + row_h_in + 0.04) if comp_metrics else 0.0
    table_h_in = label_h_in + header_h_in + n_rows * row_h_in + 0.06 + metrics_tbl_h_in + 0.06
    mol_h_in = mol_h / dpi
    fig_h_in = mol_h_in + table_h_in

    fig = plt.figure(figsize=(fig_w_in, fig_h_in), dpi=dpi)

    mol_ratio = mol_h_in / fig_h_in
    ax_mol = fig.add_axes([0, 1 - mol_ratio, 1, mol_ratio])
    ax_mol.imshow(mol_img)
    ax_mol.axis('off')
    if halogen_label:
        ax_mol.text(
            0.98, 0.98, halogen_label,
            ha='right', va='top',
            fontsize=font_size + 1, fontweight='bold',
            color='white',
            bbox=dict(boxstyle='round,pad=0.2', facecolor=_C0, edgecolor='none', alpha=0.85),
            transform=ax_mol.transAxes,
        )

    tbl_ratio = 1.0 - mol_ratio

    # --- Layout: split tbl_ratio into label, table1, gap, table2 (bottom-up) ---
    bottom_pad = 0.03 / tbl_ratio  # small padding at very bottom (in axes coords)
    metrics_frac = (metrics_tbl_h_in / table_h_in) if comp_metrics else 0.0
    gap_frac = (0.06 / table_h_in) if comp_metrics else 0.0
    fg_tbl_frac = (header_h_in + n_rows * row_h_in + 0.06) / table_h_in
    label_frac = label_h_in / table_h_in

    ax_tbl = fig.add_axes([0.01, 0, 0.98, tbl_ratio])
    ax_tbl.axis('off')
    ax_tbl.set_xlim(0, 1)
    ax_tbl.set_ylim(0, 1)

    # y positions (axes coords, from bottom=0 to top=1)
    metrics_bottom = bottom_pad
    metrics_top = metrics_bottom + metrics_frac
    gap_top = metrics_top + gap_frac
    fg_tbl_bottom = gap_top
    fg_tbl_top = fg_tbl_bottom + fg_tbl_frac
    label_bottom = fg_tbl_top

    if mol_label:
        label_center_y = label_bottom + label_frac / 2.0
        ax_tbl.text(
            0.5, label_center_y, mol_label,
            ha='center', va='center',
            fontsize=font_size + 1, fontweight='bold',
            transform=ax_tbl.transAxes,
        )

    # --- FG-specific table (table 1) ---
    col_labels_fg = ["Group", "SMARTS", "Dist.Ctr", "Dist.Per"]
    col_widths_fg = [0.35, 0.30, 0.175, 0.175]

    cell_text_fg = [
        [grp, sls or "\u2014", dct, dpr]
        for grp, sls, dct, dpr, *_ in entries
    ] or [["\u2014", "\u2014", "\u2014", "\u2014"]]

    tbl = ax_tbl.table(
        cellText=cell_text_fg,
        colLabels=col_labels_fg,
        colWidths=col_widths_fg,
        loc='upper center',
        cellLoc='left',
        bbox=[0, fg_tbl_bottom, 1, fg_tbl_frac],
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(font_size)

    for j in range(len(col_labels_fg)):
        cell = tbl[0, j]
        cell.set_facecolor(_C1)
        cell.set_text_props(color='white', fontweight='bold')
        cell.set_edgecolor(_C1)

    for i in range(len(cell_text_fg)):
        for j in range(len(col_labels_fg)):
            cell = tbl[i + 1, j]
            cell.set_facecolor(_lighter(_C1) if i % 2 == 0 else 'white')
            cell.set_edgecolor(_lighter(_C1, factor=0.55))

    # --- Component-wide metrics table (table 2) ---
    if comp_metrics:
        col_labels_m = ["Size", "Branch.", "Ecc.", "ø", "Radius", "Eff.Res.", "% C"]
        col_widths_m = [0.11, 0.13, 0.13, 0.10, 0.12, 0.19, 0.22]
        cell_text_m = [list(comp_metrics)]

        tbl2 = ax_tbl.table(
            cellText=cell_text_m,
            colLabels=col_labels_m,
            colWidths=col_widths_m,
            loc='upper center',
            cellLoc='center',
            bbox=[0, metrics_bottom, 1, metrics_frac],
        )
        tbl2.auto_set_font_size(False)
        tbl2.set_fontsize(font_size)

        for j in range(len(col_labels_m)):
            cell = tbl2[0, j]
            cell.set_facecolor(_C2)
            cell.set_text_props(color='white', fontweight='bold')
            cell.set_edgecolor(_C2)

        for j in range(len(col_labels_m)):
            cell = tbl2[1, j]
            cell.set_facecolor(_lighter(_C2))
            cell.set_edgecolor(_lighter(_C2, factor=0.55))

    fig.patch.set_facecolor('white')
    buf = _io.BytesIO()
    fig.savefig(buf, format='png', dpi=dpi, bbox_inches='tight', pad_inches=0.04,
                facecolor='white')
    plt.close(fig)
    buf.seek(0)
    return Image.open(buf).copy()


@dataclass
class ComponentView:
    """Lightweight wrapper around a matched component dict.

    This does not change the underlying structure; it only provides
    nicer attribute-style access where useful.
    """

    data: Dict[str, Any]

    @property
    def atoms(self) -> List[int]:
        return self.data.get("component", [])

    @property
    def smarts_label(self) -> Optional[str]:
        return self.data.get("SMARTS")

    @property
    def size(self) -> int:
        """Number of carbon atoms in the component (falls back to atom-index count)."""
        return self.data.get("size", len(self.atoms))

    @property
    def branching(self) -> Optional[float]:
        """Branching metric: 1.0 = linear, 0.0 = fully branched."""
        return self.data.get("branching")

    @property
    def mean_eccentricity(self) -> Optional[float]:
        """Mean graph eccentricity across nodes in the component."""
        return self.data.get("mean_eccentricity")

    @property
    def min_dist_to_center(self) -> Optional[int]:
        """Minimum graph distance from any SMARTS match atom to the component center."""
        return self.data.get("min_dist_to_center")

    @property
    def max_dist_to_periphery(self) -> Optional[int]:
        """Maximum graph distance from any SMARTS match atom to the component periphery."""
        return self.data.get("max_dist_to_periphery")

    @property
    def component_fraction(self) -> Optional[float]:
        """Fraction of total molecular carbon atoms that belong to this component.

        Computed as (# C atoms in the augmented component) / (total # C atoms in molecule).
        Oxygen, fluorine, and other heteroatoms in the augmented component are excluded
        from both numerator and denominator, so the value is always in [0, 1].
        """
        return self.data.get("component_fraction")

    @property
    def diameter(self) -> Optional[float]:
        """Graph diameter of the component (longest shortest path)."""
        return self.data.get("diameter")

    @property
    def radius(self) -> Optional[float]:
        """Graph radius of the component (minimum eccentricity)."""
        return self.data.get("radius")

    @property
    def effective_graph_resistance(self) -> Optional[float]:
        """Effective graph resistance (Kirchhoff index) of the component."""
        return self.data.get("effective_graph_resistance")

    @property
    def effective_graph_resistance_BDE(self) -> Optional[float]:
        """BDE-weighted effective graph resistance of the component."""
        return self.data.get("effective_graph_resistance_BDE")


class MatchView(dict):
    """Wrapper for a single match dict (PFAS group or definition).

    Behaves as a normal dict but provides helpers for component access.
    """

    @property
    def is_group(self) -> bool:
        return self.get("type") == "HalogenGroup"

    @property
    def is_definition(self) -> bool:
        return self.get("type") == "PFASdefinition"

    @property
    def group_id(self) -> Optional[int]:
        return self.get("id") if self.is_group else None

    @property
    def group_name(self) -> Optional[str]:
        return self.get("group_name") if self.is_group else None

    @property
    def components(self) -> List[ComponentView]:
        return [ComponentView(c) for c in self.get("components", [])]


class EmbeddingArray(np.ndarray):
    """A numpy array subclass that carries molecule identity metadata.

    Returned by :meth:`PFASEmbedding.to_array` (1-D) and
    :meth:`PFASEmbeddingSet.to_array` (2-D).  All standard numpy operations
    work unchanged; the extra attributes allow callers to trace each row back
    to its source molecule.

    Attributes
    ----------
    smiles : str or list of str
        SMILES string(s) for the molecule(s).
    inchi : str or list of str
        InChI string(s).
    inchikey : str or list of str
        InChIKey(s).
    source : PFASEmbedding or PFASEmbeddingSet
        Reference to the originating result object.
    """

    def __new__(cls, array, smiles='', inchi='', inchikey='', source=None):
        obj = np.asarray(array).view(cls)
        obj._emb_smiles   = smiles
        obj._emb_inchi    = inchi
        obj._emb_inchikey = inchikey
        obj._emb_source   = source
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self._emb_smiles   = getattr(obj, '_emb_smiles',   '')
        self._emb_inchi    = getattr(obj, '_emb_inchi',    '')
        self._emb_inchikey = getattr(obj, '_emb_inchikey', '')
        self._emb_source   = getattr(obj, '_emb_source',   None)

    @property
    def smiles(self):
        return self._emb_smiles

    @property
    def inchi(self):
        return self._emb_inchi

    @property
    def inchikey(self):
        return self._emb_inchikey

    @property
    def source(self):
        return self._emb_source


class PFASEmbedding(dict):
    """Single-molecule PFAS result and embedding generator.

    Subclasses :class:`dict` so all existing code accessing ``result['smiles']``,
    ``result['matches']``, etc. continues to work unchanged.  Call
    :meth:`to_array` to produce a numeric embedding vector from the stored
    parsed data.
    """

    @property
    def smiles(self) -> str:
        return self.get("smiles", "")

    @property
    def mol_with_h(self):
        """Get the molecule with explicit hydrogens used for component detection.

        Atom indices in the molblock are guaranteed to match those stored in
        matched component dicts.  The older SMILES path is kept only for
        backward compatibility with cached results (MolToSmiles reorders atoms).
        """
        molblock = self.get("molblock_with_h")
        if molblock:
            return Chem.MolFromMolBlock(molblock, removeHs=False)
        # Backward compatibility: SMILES path (atom ordering may differ)
        smiles_h = self.get("smiles_with_h")
        if smiles_h:
            return Chem.MolFromSmiles(smiles_h)
        return self.get("mol_with_h")

    @property
    def matches(self) -> List[MatchView]:
        return [MatchView(m) for m in self.get("matches", [])]

    def iter_group_matches(self, group_id: Optional[int] = None, group_name: Optional[str] = None) -> Iterator[MatchView]:
        """Iterate over PFAS group matches, optionally filtered by id/name."""

        for m in self.matches:
            if not m.is_group:
                continue
            if group_id is not None and m.group_id != group_id:
                continue
            if group_name is not None and m.group_name != group_name:
                continue
            yield m

    def collect_component_atoms(self, group_id: Optional[int] = None, group_name: Optional[str] = None) -> List[int]:
        """Return all atom indices that belong to matching components."""

        atoms: List[int] = []
        for m in self.iter_group_matches(group_id=group_id, group_name=group_name):
            for comp in m.components:
                atoms.extend(comp.atoms)
        # Deduplicate but keep order stable
        seen = set()
        deduped: List[int] = []
        for idx in atoms:
            if idx not in seen:
                seen.add(idx)
                deduped.append(idx)
        return deduped

    def summarise(self) -> str:
        """Return a coloured text summary of this molecule's results.

        The summary includes:
        - SMILES representation
        - counts of PFAS group and definition matches
        - total number of components across all group matches
        - list of matched PFAS groups (colour-coded by halogen)
        """
        meta = _get_component_meta()

        total_group_matches = 0
        total_definition_matches = 0
        total_components = 0
        group_counts: Dict[str, int] = {}

        for m in self.matches:
            if m.is_group:
                total_group_matches += 1
                total_components += len(m.components)
                name = m.group_name or str(m.get("match_id", ""))
                group_counts[name] = group_counts.get(name, 0) + 1
            elif m.is_definition:
                total_definition_matches += 1

        lines: List[str] = []
        lines.append(_ansi("PFASEmbedding summary", _ANSI_BOLD))
        lines.append(f"- SMILES: {self.smiles}")
        lines.append(f"- PFAS group matches: {total_group_matches}")
        lines.append(f"- PFAS definition matches: {total_definition_matches}")
        lines.append(f"- Total components: {total_components}")

        if group_counts:
            lines.append("- Matched PFAS groups:")
            for name, count in sorted(group_counts.items(), key=lambda kv: kv[1], reverse=True):
                # Determine halogen from first matching component of this group
                halogen: Optional[str] = None
                for m in self.matches:
                    if m.is_group and (m.group_name or "") == name:
                        for comp in m.components:
                            sl = comp.smarts_label
                            sl_str = str(sl) if isinstance(sl, list) else sl
                            info = meta.get(sl_str or "", {})
                            halogen = info.get("halogen")
                            break
                        break
                hal_code = _ANSI_HALOGEN.get(halogen or "", "")
                lines.append(f"  * {hal_code}{_ANSI_BOLD}{name}{_ANSI_RESET}: {count} match(es)")

        return "\n".join(lines)

    def __str__(self) -> str:
        return self.summarise()

    def summary(self) -> None:
        """Print a detailed coloured summary of matched groups and components.

        Each component is shown with graph metrics: size (carbon count),
        branching (1.0 = linear, 0.0 = highly branched) and mean eccentricity.
        """

        print("=" * 80)
        print(f"{_ANSI_BOLD}MOLECULE:{_ANSI_RESET} {self.smiles}")
        print("=" * 80)

        if not self.matches:
            print("No PFAS groups matched.")
            return

        # Collect group information
        groups_info: Dict[Tuple[Optional[int], str], List[ComponentView]] = {}

        for m in self.matches:
            if not m.is_group:
                continue
            key = (m.group_id, m.group_name or "Unknown")
            if key not in groups_info:
                groups_info[key] = []
            groups_info[key].extend(m.components)

        if not groups_info:
            print("No PFAS groups matched.")
            return

        print(f"\nMatched {len(groups_info)} PFAS group(s):")
        print()

        for (group_id, group_name), components in sorted(groups_info.items(), key=lambda x: x[0][0] or 0):
            print(f"{_ANSI_BOLD}Group {group_id}: {group_name}{_ANSI_RESET}")

            # Group components by SMARTS type
            by_smarts: Dict[Optional[str], List[ComponentView]] = {}
            for comp in components:
                smarts = comp.smarts_label
                smarts_key = str(smarts) if isinstance(smarts, list) else smarts
                if smarts_key not in by_smarts:
                    by_smarts[smarts_key] = []
                by_smarts[smarts_key].append(comp)

            # Display components by SMARTS type
            for smarts_label, comps in sorted(by_smarts.items(), key=lambda x: x[0] or ""):
                label = smarts_label or "(no label)"
                print(f"  SMARTS: {label}  ({len(comps)} component(s))")
                for comp in sorted(comps, key=lambda c: c.size, reverse=True):
                    br_str  = f"{comp.branching:.2f}"       if comp.branching        is not None else "\u2014"
                    ecc_str = f"{comp.mean_eccentricity:.2f}" if comp.mean_eccentricity is not None else "\u2014"
                    print(f"    size={_ANSI_BOLD}{comp.size}{_ANSI_RESET}  branching={br_str}  eccentricity={ecc_str}")

            print()

    def table(self) -> str:
        """Return a text table with one row per match.

        Columns:
        - match_index: 1-based match index
        - type: 'group' or 'definition'
        - name: PFAS group name or definition name
        - components: number of components (for groups)
        """

        lines: List[str] = []
        lines.append("match_index\ttype\tname\tcomponents")

        for idx, m in enumerate(self.matches, start=1):
            match_type = "group" if m.is_group else "definition"
            if m.is_group:
                name = m.group_name or str(m.get("match_id", ""))
                components = len(m.components)
            else:
                name = m.get("definition_name") or m.get("name") or str(m.get("match_id", ""))
                components = 0

            lines.append(f"{idx}\t{match_type}\t{name}\t{components}")

        return "\n".join(lines)

    def classify(self) -> Tuple[str, int]:
        """Classify the molecule's PFAS content into a category label.

        Returns
        -------
        (category, total_component_size) : Tuple[str, int]
            *category* — a short label describing the main PFAS groups present:

            * If one or more **OECD** groups are matched, returns their names
              joined by ``", "``.
            * Otherwise, returns a ``"per-"`` or ``"poly-"`` prefixed string
              listing the matched **generic** / **telomeric** group names
              (excluding groups 51 & 52), separated by ``", "``.
              ``"per-"`` is used only when **all** matched component SMARTS
              carry ``saturation='per'``; ``"poly-"`` otherwise.
              Name subsumption is applied: e.g. ``"amine"`` is suppressed
              when ``"sulfonamide"`` or ``"amide"`` is present.

            *total_component_size* — sum of :attr:`ComponentView.size` (C-atom
            count) across **all** matched group components.
        """
        group_info = _get_group_info()
        meta = _get_component_meta()

        oecd_names: List[str] = []
        non_oecd: List[Tuple[int, str]] = []  # (group_id, name)
        total_size: int = 0
        seen_oecd: set = set()
        seen_non_oecd_ids: set = set()

        for m in self.iter_group_matches():
            gid = m.group_id
            if gid is None:
                continue
            info = group_info.get(gid, {})
            cat = info.get("category", "other")
            name = m.group_name or info.get("name", f"group_{gid}")

            # Accumulate total component size across ALL groups
            for comp in m.components:
                total_size += comp.size

            if cat == "OECD":
                if name not in seen_oecd:
                    seen_oecd.add(name)
                    oecd_names.append(name)
            elif cat in ("generic", "telomer"):
                if gid not in _CLASSIFY_EXCLUDED_IDS and gid not in seen_non_oecd_ids:
                    seen_non_oecd_ids.add(gid)
                    non_oecd.append((gid, name))

        # --- OECD priority -----------------------------------------------
        if oecd_names:
            return ", ".join(oecd_names), total_size

        # --- Non-OECD (generic + telomeric) --------------------------------
        if not non_oecd:
            return "unclassified", total_size

        # Apply name subsumption
        matched_name_set = {name for _, name in non_oecd}
        suppressed: set = set()
        for name in matched_name_set:
            for s in _SUBSUMES.get(name, []):
                if s in matched_name_set:
                    suppressed.add(s)

        filtered_names: List[str] = [
            name for _, name in non_oecd if name not in suppressed
        ]

        # Determine per/poly prefix from component saturation metadata
        all_per = True
        any_sat = False
        for m in self.iter_group_matches():
            gid = m.group_id
            if gid is None or gid in _CLASSIFY_EXCLUDED_IDS:
                continue
            if group_info.get(gid, {}).get("category", "other") not in ("generic", "telomer"):
                continue
            for comp in m.components:
                sl = comp.smarts_label
                sl_str = str(sl) if isinstance(sl, list) else (sl or "")
                sat = meta.get(sl_str, {}).get("saturation")
                if sat is not None:
                    any_sat = True
                    if sat != "per":
                        all_per = False

        prefix = "per" if (any_sat and all_per) else "poly"
        label = (
            f"{prefix}-{', '.join(filtered_names)}" if filtered_names
            else "unclassified"
        )
        return label, total_size

    def show(
        self,
        display: bool = True,
        subwidth: int = 350,
        subheight: int = 350,
        ncols: int = 4,
    ) -> Image.Image:
        """Show all component combinations for this molecule in a grid plot.

        Components that share the same highlighted atoms are merged into a
        single panel.  The legend lists every PFAS group (and its
        halogen / form / saturation metadata) that maps to those atoms as a
        bullet-point list, avoiding repeated panels for the same molecular
        fragment.

        Atoms are highlighted with the colour of the first matching entry:
        - **Halogen**: cyan (F), green (Cl), amber (Br), violet (I)
        - **Form**: alkyl (full base colour) vs cyclic (warm-shifted tint)
        - **Saturation**: per- (full brightness) vs poly- (dimmed)

        Parameters
        ----------
        display : bool, default True
            Whether to display the image immediately.
        subwidth : int, default 350
            Width of each sub-image in pixels.
        subheight : int, default 350
            Minimum height of each sub-image in pixels.  Panels with many
            matching groups are automatically made taller.
        ncols : int, default 4
            Number of columns in the grid.

        Returns
        -------
        PIL.Image.Image
            Grid image containing all component visualizations.
        """
        meta = _get_component_meta()

        # Use the molecule with hydrogens that was used during component detection
        mol = self.mol_with_h
        if mol is None:
            # Fallback to reconstructing from SMILES if not available
            mol = Chem.MolFromSmiles(self.smiles)
            if mol is None:
                raise ValueError(f"Cannot parse SMILES: {self.smiles}")
            mol = Chem.AddHs(mol)

        # Collect all (atoms_key -> entries) grouping
        from collections import OrderedDict
        comp_groups: Dict = OrderedDict()

        for match in self.matches:
            if not match.is_group:
                continue
            base_label = match.group_name or match.get("match_id", "")
            for comp in match.components:
                atoms = comp.atoms
                if not atoms:
                    continue
                key = frozenset(atoms)
                sl = comp.smarts_label
                sl_str = str(sl) if isinstance(sl, list) else (sl or "")
                info = meta.get(sl_str, {})
                halogen    = info.get("halogen")
                form       = info.get("form")
                saturation = info.get("saturation")
                colour = _component_color(halogen, form, saturation)
                if key not in comp_groups:
                    # Build component-wide metrics once per unique atom set
                    import math as _math
                    size_str = str(comp.size)
                    br_v   = comp.branching
                    ecc_v  = comp.mean_eccentricity
                    diam_v = comp.diameter
                    rad_v  = comp.radius
                    egr_v  = comp.effective_graph_resistance
                    frc_v  = comp.component_fraction
                    br_str   = f"{br_v:.2f}"   if br_v   is not None else "\u2014"
                    ecc_str  = f"{ecc_v:.2f}"  if ecc_v  is not None else "\u2014"
                    diam_str = (f"{float(diam_v):.0f}"  if diam_v is not None and not (isinstance(diam_v, float) and (_math.isnan(diam_v) or _math.isinf(diam_v))) else "\u2014")
                    rad_str  = (f"{float(rad_v):.0f}"   if rad_v  is not None and not (isinstance(rad_v,  float) and (_math.isnan(rad_v)  or _math.isinf(rad_v)))  else "\u2014")
                    egr_str  = (f"{float(egr_v):.2f}"   if egr_v  is not None and not (isinstance(egr_v,  float) and (_math.isnan(egr_v)  or _math.isinf(egr_v)))  else "\u2014")
                    frc_str  = f"{frc_v*100:.0f}%" if frc_v is not None else "\u2014"
                    comp_groups[key] = {
                        'atoms': sorted(atoms),
                        'colour': colour,
                        'halogen': halogen,
                        'entries': [],
                        'comp_metrics': (size_str, br_str, ecc_str, diam_str, rad_str, egr_str, frc_str),
                    }
                # FG-specific entry: group name, SMARTS type, dist-to-center, dist-to-periphery
                dct_v = comp.min_dist_to_center
                dpr_v = comp.max_dist_to_periphery
                dct_str = str(dct_v) if dct_v is not None else "\u2014"
                dpr_str = str(dpr_v) if dpr_v is not None else "\u2014"
                entry = (base_label, sl_str, dct_str, dpr_str)
                if entry not in comp_groups[key]['entries']:
                    comp_groups[key]['entries'].append(entry)

        imgs: List[Image.Image] = []

        for data in comp_groups.values():
            atoms = data['atoms']
            colour = data['colour']
            entries = data['entries']
            comp_metrics = data.get('comp_metrics')
            halogen_lbl = data.get('halogen') or ""

            atom_colours: Dict[int, Color] = {a: colour for a in atoms}
            d2d = Draw.MolDraw2DCairo(subwidth, subheight)
            dopts = d2d.drawOptions()
            dopts.useBWAtomPalette()
            dopts.fixedBondLength = 20
            dopts.addAtomIndices = True
            dopts.addBondIndices = False
            dopts.maxFontSize = 14
            dopts.minFontSize = 12
            d2d.DrawMolecule(
                mol,
                highlightAtoms=atoms,
                highlightAtomColors=atom_colours,
            )
            d2d.FinishDrawing()
            mol_img = Image.open(BytesIO(d2d.GetDrawingText()))
            imgs.append(_mol_image_with_table(mol_img, entries, comp_metrics=comp_metrics, halogen_label=halogen_lbl))

        if not imgs:
            raise ValueError("No PFAS group components found to display.")

        grid, _, _ = _grid_images(imgs, buffer=4, ncols=ncols)
        if display:
            grid.show()
        return grid

    # Alias so callers can use either mol_result.show() or mol_result.plot()
    plot = show

    def svg(
        self,
        filename: str,
        subwidth: int = 350,
        subheight: int = 350,
        ncols: int = 4,
    ) -> str:
        """Export all component combinations to an SVG file (vector graphics).

        Components that share the same highlighted atoms are merged into a
        single panel with a bullet-point legend listing all matching groups.

        Parameters
        ----------
        filename : str
            Path to the output SVG file.
        subwidth : int, default 350
            Width of each sub-image in pixels.
        subheight : int, default 350
            Minimum height of each sub-image in pixels.
        ncols : int, default 4
            Number of columns in the grid.

        Returns
        -------
        str
            Path to the created SVG file.
        """
        import svgutils.transform as sg

        # Use the molecule with hydrogens that was used during component detection
        mol = self.mol_with_h
        if mol is None:
            # Fallback to reconstructing from SMILES if not available
            mol = Chem.MolFromSmiles(self.smiles)
            if mol is None:
                raise ValueError(f"Cannot parse SMILES: {self.smiles}")
            mol = Chem.AddHs(mol)

        # Group by unique atom set
        from collections import OrderedDict
        comp_groups: Dict = OrderedDict()

        for match in self.matches:
            if not match.is_group:
                continue
            base_label = match.group_name or match.get("match_id", "")
            for comp in match.components:
                atoms = comp.atoms
                if not atoms:
                    continue
                key = frozenset(atoms)
                sl = comp.smarts_label
                sl_str = str(sl) if isinstance(sl, list) else (sl or "")
                if key not in comp_groups:
                    import math as _math
                    br_v   = comp.branching
                    ecc_v  = comp.mean_eccentricity
                    diam_v = comp.diameter
                    rad_v  = comp.radius
                    egr_v  = comp.effective_graph_resistance
                    frc_v  = comp.component_fraction
                    br_str   = f"{br_v:.2f}"  if br_v  is not None else "\u2014"
                    ecc_str  = f"{ecc_v:.2f}" if ecc_v is not None else "\u2014"
                    diam_str = (f"{float(diam_v):.0f}"  if diam_v is not None and not (isinstance(diam_v, float) and (_math.isnan(diam_v) or _math.isinf(diam_v))) else "\u2014")
                    rad_str  = (f"{float(rad_v):.0f}"   if rad_v  is not None and not (isinstance(rad_v,  float) and (_math.isnan(rad_v)  or _math.isinf(rad_v)))  else "\u2014")
                    egr_str  = (f"{float(egr_v):.2f}"   if egr_v  is not None and not (isinstance(egr_v,  float) and (_math.isnan(egr_v)  or _math.isinf(egr_v)))  else "\u2014")
                    frc_str  = f"{frc_v*100:.0f}%" if frc_v is not None else "\u2014"
                    comp_groups[key] = {
                        'atoms': sorted(atoms),
                        'entries': [],
                        'comp_metrics': (str(comp.size), br_str, ecc_str, diam_str, rad_str, egr_str, frc_str),
                    }
                # FG-specific entry
                dct_v = comp.min_dist_to_center
                dpr_v = comp.max_dist_to_periphery
                dct_str = str(dct_v) if dct_v is not None else "\u2014"
                dpr_str = str(dpr_v) if dpr_v is not None else "\u2014"
                entry = (base_label, sl_str, dct_str, dpr_str)
                if entry not in comp_groups[key]['entries']:
                    comp_groups[key]['entries'].append(entry)

        imgs: List[str] = []

        for data in comp_groups.values():
            atoms = data['atoms']
            entries = data['entries']
            comp_metrics = data.get('comp_metrics')
            n = len(entries)

            lines: List[str] = []
            for grp, sls, dct, dpr in entries:
                line = f"\u2022 {grp}"
                if sls:
                    line += f" | {sls}"
                line += f"  dct={dct}  dpr={dpr}"
                lines.append(line)
            if comp_metrics:
                sz, br, ecc, diam, rad, egr, frc = comp_metrics
                lines.append(f"  size={sz}  br={br}  ecc={ecc}  diam={diam}  rad={rad}  egr={egr}  chain={frc}")
            legend = "\n".join(lines)

            effective_height = max(subheight, 220 + n * 38)

            d2d = Draw.MolDraw2DSVG(subwidth, effective_height)
            dopts = d2d.drawOptions()
            dopts.useBWAtomPalette()
            dopts.fixedBondLength = 20
            dopts.addAtomIndices = True
            dopts.addBondIndices = False
            dopts.maxFontSize = 14
            dopts.minFontSize = 12
            d2d.DrawMolecule(mol, legend=legend, highlightAtoms=atoms)
            d2d.FinishDrawing()
            imgs.append(d2d.GetDrawingText())

        if not imgs:
            raise ValueError("No PFAS group components found to display.")

        # Convert SVG strings to svgutils figures
        svg_figs = [sg.fromstring(img) for img in imgs]

        # Merge into grid
        from .draw_mols import merge_svg
        grid, _, _ = merge_svg(svg_figs, buffer=4, ncols=ncols)

        grid.save(filename)
        return filename

    # ------------------------------------------------------------------
    # Factory constructors
    # ------------------------------------------------------------------

    @classmethod
    def from_smiles(cls, smiles: str, **kwargs) -> "PFASEmbedding":
        """Parse a SMILES string and return a single :class:`PFASEmbedding`.

        Parameters
        ----------
        smiles : str
            SMILES string for one molecule.
        **kwargs
            Forwarded to :func:`~PFASGroups.parser.parse_smiles`
            (e.g. ``halogens``, ``saturation``, ``progress``).
        """
        from .parser import parse_smiles
        return parse_smiles(smiles, **kwargs)[0]

    @classmethod
    def from_mol(cls, mol, **kwargs) -> "PFASEmbedding":
        """Parse an RDKit molecule and return a single :class:`PFASEmbedding`.

        Parameters
        ----------
        mol : rdkit.Chem.Mol
            RDKit molecule object.
        **kwargs
            Forwarded to :func:`~PFASGroups.parser.parse_mols`.
        """
        from .parser import parse_mols
        return parse_mols([mol], **kwargs)[0]

    @classmethod
    def from_inchi(cls, inchi: str, **kwargs) -> "PFASEmbedding":
        """Parse an InChI string and return a single :class:`PFASEmbedding`.

        Parameters
        ----------
        inchi : str
            InChI string for one molecule.
        **kwargs
            Forwarded to :func:`~PFASGroups.parser.parse_mols`.
        """
        from rdkit.Chem.inchi import MolFromInchi
        mol = MolFromInchi(inchi)
        if mol is None:
            raise ValueError(f"Cannot parse InChI: {inchi!r}")
        return cls.from_mol(mol, **kwargs)

    def to_fingerprint(
        self,
        group_selection: str = 'all',
        component_metrics: Optional[List[str]] = None,
        selected_group_ids: Optional[List[int]] = None,
        halogens: Union[str, List[str]] = 'F',
        saturation: Optional[str] = 'per',
        molecule_metrics: Optional[List[str]] = None,
        pfas_groups: Optional[List[Dict]] = None,
        preset: Optional[str] = None,
        count_mode: Optional[str] = None,
        graph_metrics: Optional[List[str]] = None,
        progress: bool = False,
        **kwargs,
    ) -> np.ndarray:
        """Deprecated. Use :meth:`to_array` instead."""
        import warnings
        warnings.warn(
            "to_fingerprint() is deprecated; use to_array() instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        if count_mode is not None or graph_metrics is not None:
            if component_metrics is None:
                component_metrics = [count_mode or 'binary'] + list(graph_metrics or [])
        return self.to_array(
            component_metrics=component_metrics,
            molecule_metrics=molecule_metrics,
            group_selection=group_selection,
            selected_group_ids=selected_group_ids,
            preset=preset,
            pfas_groups=pfas_groups,
        )

    def to_array(
        self,
        component_metrics=_UNSET,
        molecule_metrics=_UNSET,
        group_selection=_UNSET,
        selected_group_ids=_UNSET,
        aggregation=_UNSET,
        preset=_UNSET,
        pfas_groups=_UNSET,
    ) -> np.ndarray:
        """Generate a 1-D embedding vector for this molecule.

        When called with no arguments, returns the last cached embedding (or
        binary by default on the first call).  Pass explicit arguments to
        override and update the cache.

        Parameters
        ----------
        component_metrics : list of str, default ['binary']
            Per-component metrics.  Count modes: 'binary', 'count',
            'max_component', 'total_component'.  Graph metrics:
            'effective_graph_resistance', 'effective_graph_resistance_BDE',
            'branching', 'mean_eccentricity', etc.
        molecule_metrics : list of str, optional
            Molecule-wide scalars appended after all component columns:
            'n_components', 'total_size', 'mean_branching', etc.
        group_selection : str, default 'all'
            'all', 'oecd', 'generic', 'telomers', or 'generic+telomers'.
        selected_group_ids : list of int, optional
            Explicit group IDs (overrides group_selection).
        aggregation : str, default 'mean'
            How to aggregate multiple matched components per group:
            'mean' or 'median'.
        preset : str, optional
            Named configuration from ``EMBEDDING_PRESETS``.
        pfas_groups : list, optional
            Custom group list (loaded from defaults when None).

        Returns
        -------
        np.ndarray
            1-D float array of length ``n_groups × len(component_metrics)
            + len(molecule_metrics)``.
        """
        # Return last cached result when called with no arguments
        _no_args = (
            component_metrics is _UNSET and molecule_metrics is _UNSET and
            group_selection is _UNSET and selected_group_ids is _UNSET and
            aggregation is _UNSET and preset is _UNSET and pfas_groups is _UNSET
        )
        if _no_args and getattr(self, '_last_array', None) is not None:
            return self._last_array

        # Resolve sentinels to actual defaults
        if component_metrics is _UNSET: component_metrics = None
        if molecule_metrics is _UNSET:  molecule_metrics = None
        if group_selection is _UNSET:   group_selection = 'all'
        if selected_group_ids is _UNSET: selected_group_ids = None
        if aggregation is _UNSET:       aggregation = 'mean'
        if preset is _UNSET:            preset = None
        if pfas_groups is _UNSET:       pfas_groups = None

        from .embeddings import FINGERPRINT_PRESETS, _COUNT_MODES
        from .getter import get_compiled_HalogenGroups

        resolved_cm: List[str] = list(component_metrics) if component_metrics else ['binary']
        resolved_mm: List[str] = list(molecule_metrics) if molecule_metrics else []

        if preset is not None:
            if preset not in FINGERPRINT_PRESETS:
                raise ValueError(f"Unknown preset: {preset!r}. Available: {sorted(FINGERPRINT_PRESETS)}")
            _p = FINGERPRINT_PRESETS[preset]
            if _p.get('component_metrics') is not None:
                resolved_cm = list(_p['component_metrics'])
            if _p.get('molecule_metrics') is not None:
                resolved_mm = list(_p['molecule_metrics'])

        if pfas_groups is None:
            pfas_groups = get_compiled_HalogenGroups()

        sel_groups = _select_groups(group_selection, selected_group_ids, pfas_groups)
        match_by_id = {m['id']: m for m in self.get('matches', [])}

        row: List[float] = []
        for m in resolved_cm:
            for g in sel_groups:
                match = match_by_id.get(g.id)
                if match is None:
                    row.append(0.0)
                elif m in _COUNT_MODES:
                    row.append(_encode_count(match, m))
                else:
                    comps = match.get('components', [])
                    row.append(_agg(comps, m, aggregation))

        all_comps = [c for m in match_by_id.values() for c in m.get('components', [])]
        for m in resolved_mm:
            row.append(_mol_metric(all_comps, m))

        result = EmbeddingArray(
            np.array(row, dtype=float),
            smiles=self.get('smiles', ''),
            inchi=self.get('inchi', ''),
            inchikey=self.get('inchikey', ''),
            source=self,
        )
        if _no_args:
            self._last_array = result
        return result

    def column_names(
        self,
        component_metrics: Optional[List[str]] = None,
        molecule_metrics: Optional[List[str]] = None,
        group_selection: str = 'all',
        selected_group_ids: Optional[List[int]] = None,
        preset: Optional[str] = None,
        pfas_groups=None,
    ) -> List[str]:
        """Return the list of column labels for :meth:`to_array` without computing values.

        Parameters match those of :meth:`to_array` (``aggregation`` is not
        relevant for column names).
        """
        from .embeddings import FINGERPRINT_PRESETS
        from .getter import get_compiled_HalogenGroups

        resolved_cm: List[str] = list(component_metrics) if component_metrics else ['binary']
        resolved_mm: List[str] = list(molecule_metrics) if molecule_metrics else []

        if preset is not None:
            _p = FINGERPRINT_PRESETS.get(preset, {})
            if _p.get('component_metrics') is not None:
                resolved_cm = list(_p['component_metrics'])
            if _p.get('molecule_metrics') is not None:
                resolved_mm = list(_p['molecule_metrics'])

        if pfas_groups is None:
            pfas_groups = get_compiled_HalogenGroups()

        sel_groups = _select_groups(group_selection, selected_group_ids, pfas_groups)

        names = [f"{g.name} [{m}]" for m in resolved_cm for g in sel_groups]
        names += [f"mol:{m}" for m in resolved_mm]
        return names

    def to_sql(
        self,
        filename: Optional[str] = None,
        dbname: Optional[str] = None,
        user: Optional[str] = None,
        password: Optional[str] = None,
        host: Optional[str] = None,
        port: Optional[int] = None,
        components_table: str = "components",
        groups_table: str = "pfas_groups_in_compound",
        if_exists: str = "append",
    ) -> None:
        """Export this molecule result to a SQL database.

        Can write to either SQLite (via filename) or PostgreSQL/MySQL (via connection parameters).

        Parameters
        ----------
        filename : str, optional
            Path to SQLite database file. If provided, uses SQLite.
        dbname : str, optional
            Database name (for PostgreSQL/MySQL).
        user : str, optional
            Database username. Defaults to os.environ['DB_USER'] if not provided.
        password : str, optional
            Database password. Defaults to os.environ['DB_PASSWORD'] if not provided.
        host : str, optional
            Database host. Defaults to os.environ.get('DB_HOST', 'localhost').
        port : int, optional
            Database port. Defaults to os.environ.get('DB_PORT', 5432 for PostgreSQL).
        components_table : str, default "components"
            Name of the table to store component-level data.
        groups_table : str, default "pfas_groups_in_compound"
            Name of the table to store PFAS group matches.
        if_exists : str, default "append"
            How to behave if tables exist: 'fail', 'replace', or 'append'.
        """
        try:
            import pandas as pd
            import sqlalchemy
        except ImportError as exc:
            raise ImportError("pandas and sqlalchemy are required for to_sql. Install with: pip install pandas sqlalchemy") from exc
        # Determine connection
        if filename:
            engine = sqlalchemy.create_engine(f"sqlite:///{filename}")
        elif dbname:
            # Get credentials from environment if not provided
            if user is None:
                user = os.environ.get('DB_USER')
            if password is None:
                password = os.environ.get('DB_PASSWORD')
            if host is None:
                host = os.environ.get('DB_HOST', 'localhost')
            if port is None:
                port = int(os.environ.get('DB_PORT', 5432))

            if not user or not password:
                raise ValueError("Database credentials required. Provide user/password or set DB_USER/DB_PASSWORD environment variables.")

            # Assuming PostgreSQL; adjust for MySQL if needed
            connection_string = f"postgresql://{user}:{password}@{host}:{port}/{dbname}"
            engine = sqlalchemy.create_engine(connection_string)
        else:
            raise ValueError("Either filename (for SQLite) or dbname (for PostgreSQL) must be provided.")

        # Prepare components data
        components_data = []
        for match in self.matches:
            if not match.is_group:
                continue
            for comp in match.components:
                smarts = comp.smarts_label
                components_data.append({
                    'smiles': self.smiles,
                    'group_id': match.group_id,
                    'group_name': match.group_name,
                    'smarts_label': str(smarts) if isinstance(smarts, list) else smarts,
                    'component_atoms': ','.join(map(str, comp.atoms)),
                })

        # Prepare groups data
        groups_data = []
        group_counts: Dict[Tuple[Optional[int], str], int] = {}
        for match in self.matches:
            if not match.is_group:
                continue
            key = (match.group_id, match.group_name or '')
            group_counts[key] = group_counts.get(key, 0) + 1

        for (group_id, group_name), count in group_counts.items():
            groups_data.append({
                'smiles': self.smiles,
                'group_id': group_id,
                'group_name': group_name,
                'match_count': count,
            })

        # Write to database
        if components_data:
            df_components = pd.DataFrame(components_data)
            df_components.to_sql(components_table, engine, if_exists=if_exists, index=False)

        if groups_data:
            df_groups = pd.DataFrame(groups_data)
            df_groups.to_sql(groups_table, engine, if_exists=if_exists, index=False)


class PFASEmbeddingSet(list):
    """List-like container for multiple :class:`PFASEmbedding` results.

    Subclasses :class:`list` so existing code that iterates over results
    continues to work.  Call :meth:`to_array` to produce a
    ``(n_molecules, n_columns)`` matrix from all stored results.
    """

    def __init__(self, iterable: Iterable[Dict[str, Any]] = ()):  # type: ignore[override]
        super().__init__(PFASEmbedding(m) if not isinstance(m, PFASEmbedding) else m for m in iterable)

    @property
    def matches(self) -> List[MatchView]:
        """Flattened list of all MatchView objects across all molecules.

        Some older code expects a ``matches`` attribute on a ResultsModel
        instance. Provide a read-only aggregated view by concatenating the
        per-molecule match lists.
        """
        out: List[MatchView] = []
        for mol_res in self:  # type: ignore[assignment]
            out.extend(mol_res.matches)
        return out

    @classmethod
    def from_raw(cls, results: Iterable[Dict[str, Any]]) -> "PFASEmbeddingSet":
        """Wrap an existing list of result dicts without changing them."""

        return cls(results)

    @classmethod
    def from_smiles(cls, smiles: Union[str, List[str]], **kwargs) -> "PFASEmbeddingSet":
        """Parse SMILES string(s) and return a :class:`PFASEmbeddingSet`.

        Parameters
        ----------
        smiles : str or list of str
            One or more SMILES strings.
        **kwargs
            Forwarded to :func:`~PFASGroups.parser.parse_smiles`
            (e.g. ``halogens``, ``saturation``, ``progress``).
        """
        from .parser import parse_smiles
        return parse_smiles(smiles, **kwargs)

    @classmethod
    def from_mols(cls, mols, **kwargs) -> "PFASEmbeddingSet":
        """Parse RDKit molecules and return a :class:`PFASEmbeddingSet`.

        Parameters
        ----------
        mols : list of rdkit.Chem.Mol
            List of RDKit molecule objects.
        **kwargs
            Forwarded to :func:`~PFASGroups.parser.parse_mols`.
        """
        from .parser import parse_mols
        return parse_mols(list(mols), **kwargs)

    @classmethod
    def from_inchis(cls, inchis: List[str], **kwargs) -> "PFASEmbeddingSet":
        """Parse InChI strings and return a :class:`PFASEmbeddingSet`.

        Parameters
        ----------
        inchis : list of str
            List of InChI strings.
        **kwargs
            Forwarded to :func:`~PFASGroups.parser.parse_mols`.
        """
        from rdkit.Chem.inchi import MolFromInchi
        mols = []
        for inchi in inchis:
            mol = MolFromInchi(inchi)
            if mol is None:
                raise ValueError(f"Cannot parse InChI: {inchi!r}")
            mols.append(mol)
        return cls.from_mols(mols, **kwargs)



    def iter_group_matches(
        self,
        group_id: Optional[int] = None,
        group_name: Optional[str] = None,
    ) -> Iterator[Tuple["PFASEmbedding", MatchView]]:
        """Iterate over all PFAS group matches across all molecules."""

        for mol_res in self:  # type: ignore[assignment]
            for match in mol_res.iter_group_matches(group_id=group_id, group_name=group_name):
                yield mol_res, match

    # --- Plotting helpers ---------------------------------------------------

    def _draw_single_molecule(
        self,
        mol: Chem.Mol,
        highlight_atoms: List[int],
        legend: str = "",
        subwidth: int = 300,
        subheight: int = 300,
        atom_colours: Optional[Dict[int, Color]] = None,
        maxFontSize: int = 14,
        minFontSize: int = 11,
    ) -> Image.Image:
        """Draw a single molecule with highlighted atoms as a PIL image.

        Parameters
        ----------
        atom_colours : dict, optional
            Mapping of atom index → (R, G, B) float colour. When provided,
            individual atoms are highlighted with their assigned colour.
            When ``None`` the default highlight colour is used for all atoms.
        maxFontSize : int, default 14
            Maximum font size for legend text.
        minFontSize : int, default 11
            Minimum font size for legend text.
        """

        d2d = Draw.MolDraw2DCairo(subwidth, subheight)
        dopts = d2d.drawOptions()
        dopts.useBWAtomPalette()
        dopts.fixedBondLength = 20
        dopts.addAtomIndices = True
        dopts.addBondIndices = False
        dopts.maxFontSize = maxFontSize
        dopts.minFontSize = minFontSize
        if atom_colours:
            d2d.DrawMolecule(
                mol,
                legend=legend,
                highlightAtoms=highlight_atoms,
                highlightAtomColors=atom_colours,
            )
        else:
            d2d.DrawMolecule(mol, legend=legend, highlightAtoms=highlight_atoms)
        d2d.FinishDrawing()
        png = d2d.GetDrawingText()
        return Image.open(BytesIO(png))

    def plot_components_for_group(
        self,
        group_id: Optional[int] = None,
        group_name: Optional[str] = None,
        max_molecules: Optional[int] = None,
        subwidth: int = 300,
        subheight: int = 300,
        ncols: int = 3,
    ) -> Tuple[Image.Image, int, int]:
        """Plot all components for a specific PFAS group across molecules.

        Either ``group_id`` or ``group_name`` (or both) can be provided to
        select the target group. Each panel corresponds to one molecule,
        with all its components for that group highlighted together.
        """

        imgs: List[Image.Image] = []
        count = 0

        for mol_res in self:  # type: ignore[assignment]
            if max_molecules is not None and count >= max_molecules:
                break

            atoms = mol_res.collect_component_atoms(group_id=group_id, group_name=group_name)
            if not atoms:
                continue

            # Use mol_with_h to preserve atom ordering from parse_groups_in_mol
            mol = mol_res.mol_with_h
            if mol is None:
                smiles = mol_res.smiles
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    continue
                mol = Chem.AddHs(mol)

            # Try to infer a label from the first matching group
            label = ""
            for m in mol_res.iter_group_matches(group_id=group_id, group_name=group_name):
                if m.group_name is not None:
                    label = m.group_name
                break

            img = self._draw_single_molecule(mol, atoms, legend=label, subwidth=subwidth, subheight=subheight)
            imgs.append(img)
            count += 1

        if not imgs:
            raise ValueError("No matching components found for the requested group.")

        return _grid_images(imgs, buffer=4, ncols=ncols)

    def show(
        self,
        display: bool = True,
        subwidth: int = 350,
        subheight: int = 350,
        ncols: int = 4,
    ) -> Image.Image:
        """Show all component combinations in a grid plot.

        Components that share the same highlighted atoms within a molecule are
        merged into a single panel.  The table below each structure lists the
        matched PFAS group, the component SMARTS type, and three graph metrics:
        size (C-atom count), branching (1.0 = linear) and mean eccentricity.

        Atoms are highlighted with the colour derived from the component SMARTS
        metadata (halogen / form / saturation) of the first entry in each panel.
        """
        meta = _get_component_meta()

        imgs: List[Image.Image] = []

        for mol_index, mol_res in enumerate(self):  # type: ignore[assignment]
            # Use the molecule with hydrogens from parsing
            mol = mol_res.mol_with_h
            if mol is None:
                smiles = mol_res.smiles
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    continue
                mol = Chem.AddHs(mol)

            # Group by unique atom set within this molecule
            from collections import OrderedDict
            comp_groups: Dict = OrderedDict()

            for match in mol_res.matches:
                if not match.is_group:
                    continue
                base_label = match.group_name or match.get("match_id", "")
                for comp in match.components:
                    atoms = comp.atoms
                    if not atoms:
                        continue
                    key = frozenset(atoms)
                    sl = comp.smarts_label
                    sl_str = str(sl) if isinstance(sl, list) else (sl or "")
                    info = meta.get(sl_str, {})
                    halogen    = info.get("halogen")
                    form       = info.get("form")
                    saturation = info.get("saturation")
                    colour = _component_color(halogen, form, saturation)
                    if key not in comp_groups:
                        import math as _math
                        br_v   = comp.branching
                        ecc_v  = comp.mean_eccentricity
                        diam_v = comp.diameter
                        rad_v  = comp.radius
                        egr_v  = comp.effective_graph_resistance
                        frc_v  = comp.component_fraction
                        br_str   = f"{br_v:.2f}"  if br_v  is not None else "\u2014"
                        ecc_str  = f"{ecc_v:.2f}" if ecc_v is not None else "\u2014"
                        diam_str = (f"{float(diam_v):.0f}"  if diam_v is not None and not (isinstance(diam_v, float) and (_math.isnan(diam_v) or _math.isinf(diam_v))) else "\u2014")
                        rad_str  = (f"{float(rad_v):.0f}"   if rad_v  is not None and not (isinstance(rad_v,  float) and (_math.isnan(rad_v)  or _math.isinf(rad_v)))  else "\u2014")
                        egr_str  = (f"{float(egr_v):.2f}"   if egr_v  is not None and not (isinstance(egr_v,  float) and (_math.isnan(egr_v)  or _math.isinf(egr_v)))  else "\u2014")
                        frc_str  = f"{frc_v*100:.0f}%" if frc_v is not None else "\u2014"
                        comp_groups[key] = {
                            'atoms': sorted(atoms),
                            'colour': colour,
                            'halogen': halogen,
                            'entries': [],
                            'comp_metrics': (str(comp.size), br_str, ecc_str, diam_str, rad_str, egr_str, frc_str),
                        }
                    dct_v = comp.min_dist_to_center
                    dpr_v = comp.max_dist_to_periphery
                    dct_str = str(dct_v) if dct_v is not None else "\u2014"
                    dpr_str = str(dpr_v) if dpr_v is not None else "\u2014"
                    entry = (base_label, sl_str, dct_str, dpr_str)
                    if entry not in comp_groups[key]['entries']:
                        comp_groups[key]['entries'].append(entry)

            for data in comp_groups.values():
                atoms = data['atoms']
                colour = data['colour']
                entries = data['entries']
                comp_metrics = data.get('comp_metrics')
                halogen_lbl = data.get('halogen') or ""

                atom_colours: Dict[int, Color] = {a: colour for a in atoms}
                d2d = Draw.MolDraw2DCairo(subwidth, subheight)
                dopts = d2d.drawOptions()
                dopts.useBWAtomPalette()
                dopts.fixedBondLength = 20
                dopts.addAtomIndices = True
                dopts.addBondIndices = False
                dopts.maxFontSize = 14
                dopts.minFontSize = 12
                d2d.DrawMolecule(
                    mol,
                    highlightAtoms=atoms,
                    highlightAtomColors=atom_colours,
                )
                d2d.FinishDrawing()
                mol_img = Image.open(BytesIO(d2d.GetDrawingText()))
                imgs.append(_mol_image_with_table(
                    mol_img, entries, comp_metrics=comp_metrics,
                    mol_label=f"mol#{mol_index + 1}", halogen_label=halogen_lbl
                ))

        if not imgs:
            raise ValueError("No PFAS group components found to display.")

        grid, _, _ = _grid_images(imgs, buffer=4, ncols=ncols)
        if display:
            grid.show()
        return grid

    # Alias so callers can use either results.show() or results.plot()
    plot = show

    def to_sql(
        self,
        filename: Optional[str] = None,
        dbname: Optional[str] = None,
        user: Optional[str] = None,
        password: Optional[str] = None,
        host: Optional[str] = None,
        port: Optional[int] = None,
        components_table: str = "components",
        groups_table: str = "pfas_groups_in_compound",
        if_exists: str = "append",
    ) -> None:
        """Export this molecule result to a SQL database.

        Can write to either SQLite (via filename) or PostgreSQL/MySQL (via connection parameters).

        Parameters
        ----------
        filename : str, optional
            Path to SQLite database file. If provided, uses SQLite.
        dbname : str, optional
            Database name (for PostgreSQL/MySQL).
        user : str, optional
            Database username. Defaults to os.environ['DB_USER'] if not provided.
        password : str, optional
            Database password. Defaults to os.environ['DB_PASSWORD'] if not provided.
        host : str, optional
            Database host. Defaults to os.environ.get('DB_HOST', 'localhost').
        port : int, optional
            Database port. Defaults to os.environ.get('DB_PORT', 5432 for PostgreSQL).
        components_table : str, default "components"
            Name of the table to store component-level data.
        groups_table : str, default "pfas_groups_in_compound"
            Name of the table to store PFAS group matches.
        if_exists : str, default "append"
            How to behave if tables exist: 'fail', 'replace', or 'append'.
        """
        try:
            import pandas as pd
            import sqlalchemy
        except ImportError as exc:
            raise ImportError("pandas and sqlalchemy are required for to_sql. Install with: pip install pandas sqlalchemy") from exc
        # Determine connection
        if filename:
            engine = sqlalchemy.create_engine(f"sqlite:///{filename}")
        elif dbname:
            # Get credentials from environment if not provided
            if user is None:
                user = os.environ.get('DB_USER')
            if password is None:
                password = os.environ.get('DB_PASSWORD')
            if host is None:
                host = os.environ.get('DB_HOST', 'localhost')
            if port is None:
                port = int(os.environ.get('DB_PORT', 5432))

            if not user or not password:
                raise ValueError("Database credentials required. Provide user/password or set DB_USER/DB_PASSWORD environment variables.")

            # Assuming PostgreSQL; adjust for MySQL if needed
            connection_string = f"postgresql://{user}:{password}@{host}:{port}/{dbname}"
            engine = sqlalchemy.create_engine(connection_string)
        else:
            raise ValueError("Either filename (for SQLite) or dbname (for PostgreSQL) must be provided.")

        # Prepare components data across all molecules in this ResultsModel
        components_data = []
        groups_data = []

        for mol_res in self:  # type: ignore[assignment]
            # local counts per molecule
            local_group_counts: Dict[Tuple[Optional[int], str], int] = {}
            for match in mol_res.matches:
                if not match.is_group:
                    continue
                for comp in match.components:
                    smarts = comp.smarts_label
                    components_data.append({
                        'smiles': mol_res.smiles,
                        'group_id': match.group_id,
                        'group_name': match.group_name,
                        'smarts_label': str(smarts) if isinstance(smarts, list) else smarts,
                        'component_atoms': ','.join(map(str, comp.atoms)),
                    })

                key = (match.group_id, match.group_name or '')
                local_group_counts[key] = local_group_counts.get(key, 0) + 1

            for (group_id, group_name), count in local_group_counts.items():
                groups_data.append({
                    'smiles': mol_res.smiles,
                    'group_id': group_id,
                    'group_name': group_name,
                    'match_count': count,
                })

        # Write to database
        if components_data:
            df_components = pd.DataFrame(components_data)
            df_components.to_sql(components_table, engine, if_exists=if_exists, index=False)

        if groups_data:
            df_groups = pd.DataFrame(groups_data)
            df_groups.to_sql(groups_table, engine, if_exists=if_exists, index=False)

    def svg(
        self,
        filename: str,
        subwidth: int = 350,
        subheight: int = 350,
        ncols: int = 4,
    ) -> str:
        """Export all component combinations to an SVG file (vector graphics).

        Components that share the same highlighted atoms within a molecule are
        merged into a single panel with a bullet-point legend.

        Parameters
        ----------
        filename : str
            Path to the output SVG file.
        subwidth : int, default 350
            Width of each sub-image in pixels.
        subheight : int, default 350
            Minimum height of each sub-image in pixels.
        ncols : int, default 4
            Number of columns in the grid.

        Returns
        -------
        str
            Path to the created SVG file.
        """
        import svgutils.transform as sg

        imgs: List[str] = []

        for mol_index, mol_res in enumerate(self):  # type: ignore[assignment]
            # Use the molecule with hydrogens from parsing
            mol = mol_res.mol_with_h
            if mol is None:
                smiles = mol_res.smiles
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    continue
                mol = Chem.AddHs(mol)

            # Group by unique atom set within this molecule
            from collections import OrderedDict
            comp_groups: Dict = OrderedDict()

            for match in mol_res.matches:
                if not match.is_group:
                    continue
                base_label = match.group_name or match.get("match_id", "")
                for comp in match.components:
                    atoms = comp.atoms
                    if not atoms:
                        continue
                    key = frozenset(atoms)
                    sl = comp.smarts_label
                    sl_str = str(sl) if isinstance(sl, list) else (sl or "")
                    if key not in comp_groups:
                        import math as _math
                        br_v   = comp.branching
                        ecc_v  = comp.mean_eccentricity
                        diam_v = comp.diameter
                        rad_v  = comp.radius
                        egr_v  = comp.effective_graph_resistance
                        frc_v  = comp.component_fraction
                        br_str   = f"{br_v:.2f}"  if br_v  is not None else "\u2014"
                        ecc_str  = f"{ecc_v:.2f}" if ecc_v is not None else "\u2014"
                        diam_str = (f"{float(diam_v):.0f}"  if diam_v is not None and not (isinstance(diam_v, float) and (_math.isnan(diam_v) or _math.isinf(diam_v))) else "\u2014")
                        rad_str  = (f"{float(rad_v):.0f}"   if rad_v  is not None and not (isinstance(rad_v,  float) and (_math.isnan(rad_v)  or _math.isinf(rad_v)))  else "\u2014")
                        egr_str  = (f"{float(egr_v):.2f}"   if egr_v  is not None and not (isinstance(egr_v,  float) and (_math.isnan(egr_v)  or _math.isinf(egr_v)))  else "\u2014")
                        frc_str  = f"{frc_v*100:.0f}%" if frc_v is not None else "\u2014"
                        comp_groups[key] = {
                            'atoms': sorted(atoms),
                            'entries': [],
                            'comp_metrics': (str(comp.size), br_str, ecc_str, diam_str, rad_str, egr_str, frc_str),
                        }
                    dct_v = comp.min_dist_to_center
                    dpr_v = comp.max_dist_to_periphery
                    dct_str = str(dct_v) if dct_v is not None else "\u2014"
                    dpr_str = str(dpr_v) if dpr_v is not None else "\u2014"
                    entry = (base_label, sl_str, dct_str, dpr_str)
                    if entry not in comp_groups[key]['entries']:
                        comp_groups[key]['entries'].append(entry)

            for data in comp_groups.values():
                atoms = data['atoms']
                entries = data['entries']
                comp_metrics = data.get('comp_metrics')
                n = len(entries)

                lines: List[str] = [f"mol#{mol_index + 1}"]
                for grp, sls, dct, dpr in entries:
                    line = f"\u2022 {grp}"
                    if sls:
                        line += f" | {sls}"
                    line += f"  dct={dct}  dpr={dpr}"
                    lines.append(line)
                if comp_metrics:
                    sz, br, ecc, diam, rad, egr, frc = comp_metrics
                    lines.append(f"  size={sz}  br={br}  ecc={ecc}  diam={diam}  rad={rad}  egr={egr}  chain={frc}")
                legend = "\n".join(lines)

                effective_height = max(subheight, 220 + (n + 1) * 38)

                d2d = Draw.MolDraw2DSVG(subwidth, effective_height)
                dopts = d2d.drawOptions()
                dopts.useBWAtomPalette()
                dopts.fixedBondLength = 20
                dopts.addAtomIndices = True
                dopts.addBondIndices = False
                dopts.maxFontSize = 14
                dopts.minFontSize = 11
                d2d.DrawMolecule(mol, legend=legend, highlightAtoms=atoms)
                d2d.FinishDrawing()
                imgs.append(d2d.GetDrawingText())

        if not imgs:
            raise ValueError("No PFAS group components found to display.")

        # Convert SVG strings to svgutils figures
        svg_figs = [sg.fromstring(img) for img in imgs]

        # Merge into grid
        from .draw_mols import merge_svg
        grid, _, _ = merge_svg(svg_figs, buffer=4, ncols=ncols)

        grid.save(filename)
        return filename

    def summarise(self) -> str:
        """Return a coloured text summary of the results.

        The summary includes:
        - number of molecules
        - counts of PFAS group and definition matches
        - total number of components across all group matches
        - the most frequent PFAS groups (colour-coded by halogen)
        """
        meta = _get_component_meta()

        total_molecules = len(self)
        total_group_matches = 0
        total_definition_matches = 0
        total_components = 0
        group_counts: Dict[str, int] = {}
        # map group name → first halogen seen
        group_halogen: Dict[str, Optional[str]] = {}

        for mol_res in self:  # type: ignore[assignment]
            for m in mol_res.matches:
                if m.is_group:
                    total_group_matches += 1
                    total_components += len(m.components)
                    name = m.group_name or str(m.get("match_id", ""))
                    group_counts[name] = group_counts.get(name, 0) + 1
                    if name not in group_halogen:
                        for comp in m.components:
                            sl = comp.smarts_label
                            sl_str = str(sl) if isinstance(sl, list) else sl
                            info = meta.get(sl_str or "", {})
                            group_halogen[name] = info.get("halogen")
                            break
                elif m.is_definition:
                    total_definition_matches += 1

        unique_groups = len(group_counts)

        lines: List[str] = []
        lines.append(_ansi("PFASEmbeddingSet summary", _ANSI_BOLD))
        lines.append(f"- Molecules: {total_molecules}")
        lines.append(
            f"- PFAS group matches: {total_group_matches} (unique groups: {unique_groups})"
        )
        lines.append(f"- PFAS definition matches: {total_definition_matches}")
        lines.append(
            f"- Total components across all PFAS group matches: {total_components}"
        )

        if group_counts:
            lines.append("- Top PFAS groups (by number of matches):")
            for name, count in sorted(
                group_counts.items(), key=lambda kv: kv[1], reverse=True
            )[:10]:
                halogen = group_halogen.get(name)
                hal_code = _ANSI_HALOGEN.get(halogen or "", "")
                lines.append(f"  * {hal_code}{_ANSI_BOLD}{name}{_ANSI_RESET}: {count} match(es)")

        return "\n".join(lines)

    def __str__(self) -> str:
        return self.summarise()

    def table(self) -> str:
        """Return a more detailed text table with one row per molecule.

        Columns:
        - index: 1-based molecule index in this ResultsModel
        - smiles: molecule SMILES
        - group_matches: number of PFAS group matches
        - definition_matches: number of PFAS definition matches
        - groups: per-molecule PFAS groups with counts, e.g.
          "Perfluoroalkyl (2); Polyfluoroalkyl (1)".
        """

        lines: List[str] = []
        lines.append(
            "index\tsmiles\tgroup_matches\tdefinition_matches\tgroups"
        )

        for idx, mol_res in enumerate(self, start=1):  # type: ignore[assignment]
            smiles = mol_res.smiles
            group_matches = [m for m in mol_res.matches if m.is_group]
            definition_matches = [m for m in mol_res.matches if m.is_definition]

            per_mol_group_counts: Dict[str, int] = {}
            for m in group_matches:
                name = m.group_name or str(m.get("match_id", ""))
                per_mol_group_counts[name] = per_mol_group_counts.get(name, 0) + 1

            if per_mol_group_counts:
                groups_summary = "; ".join(
                    f"{name} ({count})"
                    for name, count in sorted(
                        per_mol_group_counts.items(), key=lambda kv: kv[1], reverse=True
                    )
                )
            else:
                groups_summary = "None"

            lines.append(
                f"{idx}\t{smiles}\t{len(group_matches)}\t{len(definition_matches)}\t{groups_summary}"
            )

        return "\n".join(lines)

    def classify(self) -> pd.DataFrame:
        """Return a classification DataFrame with one row per molecule.

        Each molecule is classified by :meth:`MoleculeResult.classify`.

        Returns
        -------
        pandas.DataFrame
            Columns:

            * ``smiles`` — molecule SMILES.
            * ``category`` — classification label: OECD group name(s) if
              matched, otherwise ``"per-"``/``"poly-"`` + generic/telomeric
              group names (comma-separated).
            * ``total_component_size`` — sum of C-atom counts across all
              matched group components.
        """
        rows = []
        for mol_res in self:  # type: ignore[assignment]
            category, total_size = mol_res.classify()
            rows.append({
                "smiles":               mol_res.smiles,
                "category":             category,
                "total_component_size": total_size,
            })
        return pd.DataFrame(rows, columns=["smiles", "category", "total_component_size"])

    def summary(self) -> None:
        """Print a detailed coloured summary of matched groups across all molecules.

        For each group, shows the component SMARTS type and, per component,
        the graph metrics: size (C-atom count), branching and mean eccentricity.
        Component size statistics (min, max, mean) are also shown.
        """

        print("=" * 80)
        print(f"{_ANSI_BOLD}RESULTS SUMMARY:{_ANSI_RESET} {len(self)} molecule(s)")
        print("=" * 80)

        if not self:
            print("No molecules in results.")
            return

        # Collect all group matches across all molecules
        all_groups_info: Dict[Tuple[Optional[int], str], List[ComponentView]] = {}

        for mol_res in self:  # type: ignore[assignment]
            for m in mol_res.matches:
                if not m.is_group:
                    continue
                key = (m.group_id, m.group_name or "Unknown")
                if key not in all_groups_info:
                    all_groups_info[key] = []
                all_groups_info[key].extend(m.components)

        if not all_groups_info:
            print("\nNo PFAS groups matched across all molecules.")
            return

        print(f"\nMatched {len(all_groups_info)} unique PFAS group(s) across all molecules:")
        print()

        for (group_id, group_name), components in sorted(all_groups_info.items(), key=lambda x: x[0][0] or 0):
            print(f"{_ANSI_BOLD}Group {group_id}: {group_name}{_ANSI_RESET}")

            # Group components by SMARTS type
            by_smarts: Dict[Optional[str], List[ComponentView]] = {}
            for comp in components:
                smarts = comp.smarts_label
                smarts_key = str(smarts) if isinstance(smarts, list) else smarts
                if smarts_key not in by_smarts:
                    by_smarts[smarts_key] = []
                by_smarts[smarts_key].append(comp)

            # Display components by SMARTS type
            for smarts_label, comps in sorted(by_smarts.items(), key=lambda x: x[0] or ""):
                label = smarts_label or "(no label)"
                sizes = [comp.size for comp in comps]
                min_sz = min(sizes) if sizes else 0
                max_sz = max(sizes) if sizes else 0
                mean_sz = sum(sizes) / len(sizes) if sizes else 0.0
                print(f"  SMARTS: {label}  ({len(comps)} component(s))  "
                      f"size {min_sz}\u2013{max_sz} (mean {mean_sz:.1f})")
                for comp in sorted(comps, key=lambda c: c.size, reverse=True):
                    br_str  = f"{comp.branching:.2f}"       if comp.branching        is not None else "\u2014"
                    ecc_str = f"{comp.mean_eccentricity:.2f}" if comp.mean_eccentricity is not None else "\u2014"
                    print(f"    size={_ANSI_BOLD}{comp.size}{_ANSI_RESET}  branching={br_str}  eccentricity={ecc_str}")

            print()

        print("=" * 80)

    def plot_all_components_with_group_colours(
        self,
        max_molecules: Optional[int] = None,
        subwidth: int = 300,
        subheight: int = 300,
        ncols: int = 3,
    ) -> Tuple[Image.Image, int, int]:
        """Plot all matched components, coloured by PFAS group.

        Each panel corresponds to one molecule; atoms are highlighted with
        colours assigned per PFAS group. The legend lists the groups found
        in that molecule.
        """

        imgs: List[Image.Image] = []
        count = 0

        for mol_res in self:  # type: ignore[assignment]
            if max_molecules is not None and count >= max_molecules:
                break

            # Use the molecule with hydrogens from parsing
            mol = mol_res.mol_with_h
            if mol is None:
                # Fallback to SMILES
                smiles = mol_res.smiles
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    continue
                mol = Chem.AddHs(mol)

            # Build colour map per group
            atom_colours: Dict[int, Color] = {}
            group_labels: List[str] = []
            group_index: Dict[str, int] = {}

            for match in mol_res.matches:
                if not match.is_group:
                    continue
                gid = match.get("match_id") or f"G{match.group_id}"  # type: ignore[operator]
                if gid not in group_index:
                    group_index[gid] = len(group_index)
                colour = _GROUP_COLORS[group_index[gid] % len(_GROUP_COLORS)]

                label_name = match.group_name or gid
                if label_name not in group_labels:
                    group_labels.append(label_name)

                for comp in match.components:
                    for atom_idx in comp.atoms:
                        atom_colours.setdefault(atom_idx, colour)

            if not atom_colours:
                continue

            # Legend lists all group names present in this molecule
            legend = ", ".join(group_labels)

            d2d = Draw.MolDraw2DCairo(subwidth, subheight)
            dopts = d2d.drawOptions()
            dopts.useBWAtomPalette()
            dopts.fixedBondLength = 20
            dopts.addAtomIndices = True
            dopts.addBondIndices = False
            dopts.maxFontSize = 16
            dopts.minFontSize = 13

            highlight_atoms = list(atom_colours.keys())
            d2d.DrawMolecule(mol, legend=legend, highlightAtoms=highlight_atoms, highlightAtomColors=atom_colours)
            d2d.FinishDrawing()
            png = d2d.GetDrawingText()
            imgs.append(Image.open(BytesIO(png)))
            count += 1

        if not imgs:
            raise ValueError("No PFAS group components found in results.")

        return _grid_images(imgs, buffer=4, ncols=ncols)

    def to_sql_all(
        self,
        conn: Optional[Union[str, 'sqlalchemy.engine.Engine']] = None,
        filename: Optional[str] = None,
        components_table: str = "components",
        groups_table: str = "pfas_groups_in_compound",
        if_exists: str = "append",
    ) -> None:
        """Export all molecule results to a SQL database.

        This method efficiently batches all molecules into the database in a single operation.

        Parameters
        ----------
        conn : str or sqlalchemy.engine.Engine, optional
            Database connection. Can be:
            - SQLAlchemy Engine object
            - Connection string (e.g., 'postgresql://user:pass@host:port/db')
            - SQLite path with 'sqlite:///' prefix
        filename : str, optional
            Path to SQLite database file (legacy parameter, use conn instead).
        components_table : str, default "components"
            Name of the table to store component-level data.
        groups_table : str, default "pfas_groups_in_compound"
            Name of the table to store PFAS group matches.
        if_exists : str, default "append"
            How to behave if tables exist: 'fail', 'replace', or 'append'.

        Examples
        --------
        >>> # Using connection string
        >>> results.to_sql(conn='postgresql://user:pass@localhost/pfas_db')
        >>>
        >>> # Using SQLAlchemy engine
        >>> from sqlalchemy import create_engine
        >>> engine = create_engine('sqlite:///pfas.db')
        >>> results.to_sql_all(conn=engine)
        >>>
        >>> # Using filename (legacy)
        >>> results.to_sql(filename='pfas.db')
        """
        try:
            import pandas as pd
            import sqlalchemy
        except ImportError as exc:
            raise ImportError("pandas and sqlalchemy are required for to_sql. Install with: pip install pandas sqlalchemy") from exc
        # Determine connection
        if conn is None and filename is None:
            raise ValueError("Either 'conn' or 'filename' must be provided.")

        if conn is not None:
            # Handle conn parameter
            if isinstance(conn, str):
                # If it's a string, create engine from connection string
                engine = sqlalchemy.create_engine(conn)
            else:
                # Assume it's already a SQLAlchemy engine
                engine = conn
        else:
            # Legacy filename parameter
            engine = sqlalchemy.create_engine(f"sqlite:///{filename}")

        # Prepare components data for all molecules
        components_data = []
        for mol_res in self:  # type: ignore[assignment]
            for match in mol_res.matches:
                if not match.is_group:
                    continue
                for comp in match.components:
                    smarts = comp.smarts_label
                    components_data.append({
                        'smiles': mol_res.smiles,
                        'group_id': match.group_id,
                        'group_name': match.group_name,
                        'smarts_label': str(smarts) if isinstance(smarts, list) else smarts,
                        'component_atoms': ','.join(map(str, comp.atoms)),
                    })

        # Prepare groups data for all molecules
        groups_data = []
        for mol_res in self:  # type: ignore[assignment]
            group_counts: Dict[Tuple[Optional[int], str], int] = {}
            for match in mol_res.matches:
                if not match.is_group:
                    continue
                key = (match.group_id, match.group_name or '')
                group_counts[key] = group_counts.get(key, 0) + 1

            for (group_id, group_name), count in group_counts.items():
                groups_data.append({
                    'smiles': mol_res.smiles,
                    'group_id': group_id,
                    'group_name': group_name,
                    'match_count': count,
                })

        # Write to database
        if components_data:
            df_components = pd.DataFrame(components_data)
            df_components.to_sql(components_table, engine, if_exists=if_exists, index=False)

        if groups_data:
            df_groups = pd.DataFrame(groups_data)
            df_groups.to_sql(groups_table, engine, if_exists=if_exists, index=False)
    def to_fingerprint(
        self,
        group_selection: str = 'all',
        component_metrics: Optional[List[str]] = None,
        selected_group_ids: Optional[List[int]] = None,
        halogens: Union[str, List[str]] = 'F',
        saturation: Optional[str] = 'per',
        molecule_metrics: Optional[List[str]] = None,
        pfas_groups: Optional[List[Dict]] = None,
        preset: Optional[str] = None,
        count_mode: Optional[str] = None,
        graph_metrics: Optional[List[str]] = None,
        progress: bool = False,
        **kwargs,
    ) -> np.ndarray:
        """Deprecated. Use :meth:`to_array` instead."""
        import warnings
        warnings.warn(
            "to_fingerprint() is deprecated; use to_array() instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        if count_mode is not None or graph_metrics is not None:
            if component_metrics is None:
                component_metrics = [count_mode or 'binary'] + list(graph_metrics or [])
        return self.to_array(
            component_metrics=component_metrics,
            molecule_metrics=molecule_metrics,
            group_selection=group_selection,
            selected_group_ids=selected_group_ids,
            preset=preset,
            pfas_groups=pfas_groups,
        )

    def to_array(
        self,
        component_metrics=_UNSET,
        molecule_metrics=_UNSET,
        group_selection=_UNSET,
        selected_group_ids=_UNSET,
        aggregation=_UNSET,
        preset=_UNSET,
        pfas_groups=_UNSET,
        progress: bool = True,
    ) -> "EmbeddingArray":
        """Stack per-molecule embedding rows into a ``(n_mols, n_cols)`` matrix.

        When called with no arguments, returns the last cached embedding (or
        binary by default on the first call).  Pass explicit arguments to
        override and update the cache.

        Parameters match those of :meth:`PFASEmbedding.to_array`, plus:

        progress : bool, default True
            Show a tqdm progress bar while computing embeddings.
        """
        _no_args = (
            component_metrics is _UNSET and molecule_metrics is _UNSET and
            group_selection is _UNSET and selected_group_ids is _UNSET and
            aggregation is _UNSET and preset is _UNSET and pfas_groups is _UNSET
        )
        if _no_args and getattr(self, '_last_array', None) is not None:
            return self._last_array

        # Resolve sentinels to defaults
        if component_metrics is _UNSET: component_metrics = None
        if molecule_metrics is _UNSET:  molecule_metrics = None
        if group_selection is _UNSET:   group_selection = 'all'
        if selected_group_ids is _UNSET: selected_group_ids = None
        if aggregation is _UNSET:       aggregation = 'mean'
        if preset is _UNSET:            preset = None
        if pfas_groups is _UNSET:       pfas_groups = None

        from .getter import get_compiled_HalogenGroups
        if pfas_groups is None:
            pfas_groups = get_compiled_HalogenGroups()

        _iter = self
        if progress and len(self) > 1:
            try:
                from tqdm.auto import tqdm as _tqdm
            except ImportError:
                from tqdm import tqdm as _tqdm
            _iter = _tqdm(self, desc='Computing embeddings', total=len(self))

        rows = [
            mol.to_array(
                component_metrics=component_metrics,
                molecule_metrics=molecule_metrics,
                group_selection=group_selection,
                selected_group_ids=selected_group_ids,
                aggregation=aggregation,
                preset=preset,
                pfas_groups=pfas_groups,
            )
            for mol in _iter
        ]
        if not rows:
            mat = np.zeros((0, 0), dtype=float)
        else:
            mat = np.vstack(rows)

        all_smiles   = [m.get('smiles', '')   for m in self]
        all_inchi    = [m.get('inchi', '')    for m in self]
        all_inchikey = [m.get('inchikey', '') for m in self]
        result = EmbeddingArray(mat, smiles=all_smiles, inchi=all_inchi,
                                inchikey=all_inchikey, source=self)
        if _no_args:
            self._last_array = result
        return result

    # ------------------------------------------------------------------
    # Analysis methods
    # ------------------------------------------------------------------

    def compare_kld(
        self,
        other: "PFASEmbeddingSet",
        method: str = 'minmax',
    ) -> float:
        """Compare two sets using KL divergence on group-occurrence frequencies.

        Parameters
        ----------
        other : PFASEmbeddingSet
            Second set to compare against.
        method : str, default ``'minmax'``
            ``'forward'``, ``'reverse'``, ``'symmetric'``, or ``'minmax'``
            (normalised symmetric KLD).

        Returns
        -------
        float
            KL divergence value (lower = more similar).
        """
        from scipy.stats import entropy as _entropy

        p_mat = self.to_array()
        q_mat = other.to_array()

        if p_mat.shape[1] != q_mat.shape[1]:
            raise ValueError(
                f"Column count mismatch: {p_mat.shape[1]} vs {q_mat.shape[1]}. "
                "Call to_array() with explicit arguments on both sets to ensure "
                "the same embedding configuration."
            )

        eps = 1e-10
        p = np.sum(np.atleast_2d(np.asarray(p_mat)) > 0, axis=0).astype(float) + eps
        q = np.sum(np.atleast_2d(np.asarray(q_mat)) > 0, axis=0).astype(float) + eps
        p /= p.sum()
        q /= q.sum()

        if method == 'forward':
            return float(_entropy(p, q))
        if method == 'reverse':
            return float(_entropy(q, p))
        if method == 'symmetric':
            return float((_entropy(p, q) + _entropy(q, p)) / 2)
        if method == 'minmax':
            kl_fwd = _entropy(p, q)
            kl_rev = _entropy(q, p)
            kl_sym = (kl_fwd + kl_rev) / 2
            max_kl = np.log(len(p))
            return float(kl_sym / max_kl) if max_kl > 0 else 0.0
        raise ValueError(f"Unknown method: {method!r}. Choose from 'forward', 'reverse', 'symmetric', 'minmax'.")

    def perform_pca(
        self,
        n_components: int = 2,
        plot: bool = True,
        output_file: Optional[str] = None,
    ) -> Dict:
        """Perform PCA on the embedding matrix.

        Parameters
        ----------
        n_components : int, default 2
        plot : bool, default True
        output_file : str, optional

        Returns
        -------
        dict
            Keys: ``'transformed'``, ``'explained_variance'``, ``'components'``,
            ``'pca_model'``, ``'scaler'``.
        """
        try:
            from sklearn.decomposition import PCA
            from sklearn.preprocessing import StandardScaler
            import matplotlib
            import matplotlib.pyplot as plt
        except ImportError as exc:
            raise ImportError(
                "scikit-learn and matplotlib required: pip install scikit-learn matplotlib"
            ) from exc

        mat = np.atleast_2d(np.asarray(self.to_array()))
        scaler = StandardScaler()
        X = scaler.fit_transform(mat)
        pca = PCA(n_components=n_components)
        X_pca = pca.fit_transform(X)

        if plot and n_components >= 2:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
            ax1.scatter(X_pca[:, 0], X_pca[:, 1], alpha=0.6, s=50)
            ax1.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
            ax1.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
            ax1.set_title('PCA of PFAS Group Embeddings')
            ax1.grid(True, alpha=0.3)
            ax2.bar(range(1, n_components + 1), pca.explained_variance_ratio_)
            ax2.set_xlabel('Principal Component')
            ax2.set_ylabel('Explained Variance Ratio')
            ax2.set_title('Scree Plot')
            ax2.grid(True, alpha=0.3)
            plt.tight_layout()
            if output_file:
                plt.savefig(output_file, dpi=300, bbox_inches='tight')
            elif matplotlib.get_backend() != 'agg':
                plt.show()
            plt.close()

        return {
            'transformed': X_pca,
            'explained_variance': pca.explained_variance_ratio_,
            'components': pca.components_,
            'pca_model': pca,
            'scaler': scaler,
        }

    def perform_kernel_pca(
        self,
        n_components: int = 2,
        kernel: str = 'rbf',
        gamma: Optional[float] = None,
        plot: bool = True,
        output_file: Optional[str] = None,
    ) -> Dict:
        """Perform kernel PCA on the embedding matrix.

        Parameters
        ----------
        n_components : int, default 2
        kernel : str, default ``'rbf'``
        gamma : float, optional
        plot : bool, default True
        output_file : str, optional

        Returns
        -------
        dict
            Keys: ``'transformed'``, ``'kpca_model'``, ``'scaler'``, ``'kernel'``, ``'gamma'``.
        """
        try:
            from sklearn.decomposition import KernelPCA
            from sklearn.preprocessing import StandardScaler
            import matplotlib
            import matplotlib.pyplot as plt
        except ImportError as exc:
            raise ImportError(
                "scikit-learn and matplotlib required: pip install scikit-learn matplotlib"
            ) from exc

        mat = np.atleast_2d(np.asarray(self.to_array()))
        scaler = StandardScaler()
        X = scaler.fit_transform(mat)
        gamma = gamma if gamma is not None else 1.0 / X.shape[1]
        kpca = KernelPCA(n_components=n_components, kernel=kernel, gamma=gamma)
        X_kpca = kpca.fit_transform(X)

        if plot and n_components >= 2:
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.scatter(X_kpca[:, 0], X_kpca[:, 1], alpha=0.6, s=50)
            ax.set_xlabel('Kernel PC1')
            ax.set_ylabel('Kernel PC2')
            ax.set_title(f'Kernel PCA ({kernel}) of PFAS Group Embeddings')
            ax.grid(True, alpha=0.3)
            plt.tight_layout()
            if output_file:
                plt.savefig(output_file, dpi=300, bbox_inches='tight')
            elif matplotlib.get_backend() != 'agg':
                plt.show()
            plt.close()

        return {
            'transformed': X_kpca,
            'kpca_model': kpca,
            'scaler': scaler,
            'kernel': kernel,
            'gamma': gamma,
        }

    def perform_tsne(
        self,
        n_components: int = 2,
        perplexity: float = 30.0,
        learning_rate: float = 200.0,
        max_iter: int = 1000,
        plot: bool = True,
        output_file: Optional[str] = None,
    ) -> Dict:
        """Perform t-SNE dimensionality reduction on the embedding matrix.

        Parameters
        ----------
        n_components : int, default 2
        perplexity : float, default 30.0
        learning_rate : float, default 200.0
        max_iter : int, default 1000
        plot : bool, default True
        output_file : str, optional

        Returns
        -------
        dict
            Keys: ``'transformed'``, ``'tsne_model'``, ``'scaler'``, ``'perplexity'``.
        """
        try:
            from sklearn.manifold import TSNE
            from sklearn.preprocessing import StandardScaler
            import matplotlib
            import matplotlib.pyplot as plt
        except ImportError as exc:
            raise ImportError(
                "scikit-learn and matplotlib required: pip install scikit-learn matplotlib"
            ) from exc

        mat = np.atleast_2d(np.asarray(self.to_array()))
        scaler = StandardScaler()
        X = scaler.fit_transform(mat)
        tsne = TSNE(
            n_components=n_components,
            perplexity=perplexity,
            learning_rate=learning_rate,
            max_iter=max_iter,
            random_state=42,
        )
        X_tsne = tsne.fit_transform(X)

        if plot and n_components >= 2:
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.scatter(X_tsne[:, 0], X_tsne[:, 1], alpha=0.6, s=50)
            ax.set_xlabel('t-SNE 1')
            ax.set_ylabel('t-SNE 2')
            ax.set_title(f't-SNE (perplexity={perplexity}) of PFAS Group Embeddings')
            ax.grid(True, alpha=0.3)
            plt.tight_layout()
            if output_file:
                plt.savefig(output_file, dpi=300, bbox_inches='tight')
            elif matplotlib.get_backend() != 'agg':
                plt.show()
            plt.close()

        return {
            'transformed': X_tsne,
            'tsne_model': tsne,
            'scaler': scaler,
            'perplexity': perplexity,
        }

    def perform_umap(
        self,
        n_components: int = 2,
        n_neighbors: int = 15,
        min_dist: float = 0.1,
        metric: str = 'euclidean',
        plot: bool = True,
        output_file: Optional[str] = None,
    ) -> Dict:
        """Perform UMAP dimensionality reduction on the embedding matrix.

        Parameters
        ----------
        n_components : int, default 2
        n_neighbors : int, default 15
        min_dist : float, default 0.1
        metric : str, default ``'euclidean'``
        plot : bool, default True
        output_file : str, optional

        Returns
        -------
        dict
            Keys: ``'transformed'``, ``'umap_model'``, ``'scaler'``, ``'n_neighbors'``, ``'min_dist'``.
        """
        try:
            import umap
            from sklearn.preprocessing import StandardScaler
            import matplotlib
            import matplotlib.pyplot as plt
        except ImportError as exc:
            raise ImportError(
                "umap-learn and matplotlib required: pip install umap-learn matplotlib"
            ) from exc

        import warnings as _warnings
        import os as _os
        _os.environ.setdefault('KMP_DUPLICATE_LIB_OK', 'TRUE')

        mat = np.atleast_2d(np.asarray(self.to_array()))
        scaler = StandardScaler()
        X = scaler.fit_transform(mat)

        with _warnings.catch_warnings():
            _warnings.filterwarnings('ignore', message=r'.*n_jobs.*overridden.*', category=UserWarning)
            _warnings.filterwarnings('ignore', message=r'.*Intel OpenMP.*LLVM OpenMP.*', category=RuntimeWarning)
            reducer = umap.UMAP(
                n_components=n_components,
                n_neighbors=n_neighbors,
                min_dist=min_dist,
                metric=metric,
                random_state=42,
            )
            X_umap = reducer.fit_transform(X)

        if plot and n_components >= 2:
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.scatter(X_umap[:, 0], X_umap[:, 1], alpha=0.6, s=50)
            ax.set_xlabel('UMAP 1')
            ax.set_ylabel('UMAP 2')
            ax.set_title(f'UMAP (n_neighbors={n_neighbors}) of PFAS Group Embeddings')
            ax.grid(True, alpha=0.3)
            plt.tight_layout()
            if output_file:
                plt.savefig(output_file, dpi=300, bbox_inches='tight')
            elif matplotlib.get_backend() != 'agg':
                plt.show()
            plt.close()

        return {
            'transformed': X_umap,
            'umap_model': reducer,
            'scaler': scaler,
            'n_neighbors': n_neighbors,
            'min_dist': min_dist,
        }

    def column_names(
        self,
        component_metrics: Optional[List[str]] = None,
        molecule_metrics: Optional[List[str]] = None,
        group_selection: str = 'all',
        selected_group_ids: Optional[List[int]] = None,
        preset: Optional[str] = None,
        pfas_groups=None,
    ) -> List[str]:
        """Return column labels (delegates to first element)."""
        if self:
            return self[0].column_names(
                component_metrics=component_metrics,
                molecule_metrics=molecule_metrics,
                group_selection=group_selection,
                selected_group_ids=selected_group_ids,
                preset=preset,
                pfas_groups=pfas_groups,
            )
        return []

    @classmethod
    def from_sql(
        cls,
        conn: Optional[Union[str, 'sqlalchemy.engine.Engine']] = None,
        filename: Optional[str] = None,
        components_table: str = "components",
        groups_table: str = "pfas_groups_in_compound",
        limit: Optional[int] = None,
    ) -> "PFASEmbeddingSet":
        """Load results from SQL database.

        Parameters
        ----------
        conn : str or SQLAlchemy Engine, optional
            Database connection string or engine
        filename : str, optional
            SQLite database filename (alternative to conn)
        components_table : str, default "components"
            Name of the components table
        groups_table : str, default "pfas_groups_in_compound"
            Name of the groups table
        limit : int, optional
            Limit number of molecules to load

        Returns
        -------
        ResultsModel
            Loaded results
        """
        try:
            import sqlalchemy
        except ImportError as exc:
            raise ImportError("sqlalchemy is required for SQL operations. Install with: pip install sqlalchemy") from exc
        # Create engine
        if conn is not None:
            if isinstance(conn, str):
                engine = sqlalchemy.create_engine(conn)
            else:
                engine = conn
        elif filename is not None:
            engine = sqlalchemy.create_engine(f"sqlite:///{filename}")
        else:
            raise ValueError("Either conn or filename must be provided")

        # Load groups data
        query = f"SELECT * FROM {groups_table}"
        if limit is not None:
            query += f" LIMIT {limit}"

        df_groups = pd.read_sql(query, engine)

        # Reconstruct results
        results = []
        for smiles in df_groups['smiles'].unique():
            mol_groups = df_groups[df_groups['smiles'] == smiles]

            matches = []
            for _, row in mol_groups.iterrows():
                matches.append({
                    'type': 'HalogenGroup',
                    'match_id': row['group_id'],
                    'group_id': row['group_id'],
                    'group_name': row['group_name'],
                    'match_count': row['match_count'],
                    'components': [],  # Components not stored in basic SQL format
                })

            results.append({
                'smiles': smiles,
                'matches': matches,
            })

        return cls(results)


# Backward-compatible aliases
MoleculeResult = PFASEmbedding
ResultsModel = PFASEmbeddingSet
