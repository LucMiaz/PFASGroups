from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple, Union
import os

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
from .parser import load_HalogenGroups

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

# Highlight colours by halogen element
_HALOGEN_COLORS: Dict[str, Color] = {
    "F":  (0.20, 0.70, 0.95),  # cyan-blue
    "Cl": (0.20, 0.80, 0.30),  # green
    "Br": (0.95, 0.75, 0.10),  # amber
    "I":  (0.80, 0.20, 0.85),  # violet
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
    (0.90, 0.10, 0.10),  # red
    (0.10, 0.40, 0.90),  # blue
    (0.10, 0.70, 0.10),  # green
    (0.90, 0.60, 0.10),  # orange
    (0.60, 0.10, 0.70),  # purple
    (0.00, 0.70, 0.70),  # teal
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
    comp_metrics: Optional[Tuple[str, str, str, str]] = None,
    mol_label: str = "",
    font_size: int = 9,
) -> Image.Image:
    """Composite a molecule image (PIL) with two formatted matplotlib tables.

    Parameters
    ----------
    mol_img : PIL Image
        The molecule drawing produced by RDKit (no legend text).
    entries : list of tuples
        Each tuple: (group_name, smarts_label, dist_center, dist_periphery, chain_pct)
        FG-specific metrics per matched component row.
    comp_metrics : tuple of str or None
        Component-wide row: (size, branching, eccentricity, diameter).
        When provided, rendered in a separate second table below the FG table.
    mol_label : str
        Optional header shown above the table (e.g. "mol#1").
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
    col_labels_fg = ["Group", "SMARTS", "Dist.Ctr", "Dist.Per", "% chain"]
    col_widths_fg = [0.35, 0.30, 0.11, 0.11, 0.13]

    cell_text_fg = [
        [grp, sls or "\u2014", dct, dpr, frc]
        for grp, sls, dct, dpr, frc in entries
    ] or [["\u2014", "\u2014", "\u2014", "\u2014", "\u2014"]]

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
        cell.set_facecolor('#3B5BA5')
        cell.set_text_props(color='white', fontweight='bold')
        cell.set_edgecolor('#3B5BA5')

    for i in range(len(cell_text_fg)):
        for j in range(len(col_labels_fg)):
            cell = tbl[i + 1, j]
            cell.set_facecolor('#EEF2FF' if i % 2 == 0 else 'white')
            cell.set_edgecolor('#C0C8E8')

    # --- Component-wide metrics table (table 2) ---
    if comp_metrics:
        col_labels_m = ["Size", "Branch.", "Ecc.", "Diam."]
        col_widths_m = [0.22, 0.26, 0.26, 0.26]
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
            cell.set_facecolor('#2E7D32')
            cell.set_text_props(color='white', fontweight='bold')
            cell.set_edgecolor('#2E7D32')

        for j in range(len(col_labels_m)):
            cell = tbl2[1, j]
            cell.set_facecolor('#E8F5E9')
            cell.set_edgecolor('#A5D6A7')

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


class MoleculeResult(dict):
    """Wrapper for a single molecule result entry.

    This subclasses ``dict`` so existing code that expects a mapping
    continues to work unchanged.
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
        lines.append(_ansi("MoleculeResult summary", _ANSI_BOLD))
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
                    br_v  = comp.branching
                    ecc_v = comp.mean_eccentricity
                    diam_v = comp.diameter
                    br_str   = f"{br_v:.2f}"   if br_v   is not None else "\u2014"
                    ecc_str  = f"{ecc_v:.2f}"  if ecc_v  is not None else "\u2014"
                    diam_str = (f"{float(diam_v):.0f}" if diam_v is not None and not _math.isnan(float(diam_v)) else "\u2014")
                    comp_groups[key] = {
                        'atoms': sorted(atoms),
                        'colour': colour,
                        'entries': [],
                        'comp_metrics': (size_str, br_str, ecc_str, diam_str),
                    }
                # FG-specific entry: group name, SMARTS type, dist-to-center, dist-to-periphery, % chain
                dct_v = comp.min_dist_to_center
                dpr_v = comp.max_dist_to_periphery
                frc_v = comp.component_fraction
                dct_str = str(dct_v) if dct_v is not None else "\u2014"
                dpr_str = str(dpr_v) if dpr_v is not None else "\u2014"
                frc_str = f"{frc_v*100:.0f}%" if frc_v is not None else "\u2014"
                entry = (base_label, sl_str, dct_str, dpr_str, frc_str)
                if entry not in comp_groups[key]['entries']:
                    comp_groups[key]['entries'].append(entry)

        imgs: List[Image.Image] = []

        for data in comp_groups.values():
            atoms = data['atoms']
            colour = data['colour']
            entries = data['entries']
            comp_metrics = data.get('comp_metrics')

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
            imgs.append(_mol_image_with_table(mol_img, entries, comp_metrics=comp_metrics))

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
                    br_str   = f"{br_v:.2f}"  if br_v  is not None else "\u2014"
                    ecc_str  = f"{ecc_v:.2f}" if ecc_v is not None else "\u2014"
                    diam_str = (f"{float(diam_v):.0f}" if diam_v is not None and not _math.isnan(float(diam_v)) else "\u2014")
                    comp_groups[key] = {
                        'atoms': sorted(atoms),
                        'entries': [],
                        'comp_metrics': (str(comp.size), br_str, ecc_str, diam_str),
                    }
                # FG-specific entry
                dct_v = comp.min_dist_to_center
                dpr_v = comp.max_dist_to_periphery
                frc_v = comp.component_fraction
                dct_str = str(dct_v) if dct_v is not None else "\u2014"
                dpr_str = str(dpr_v) if dpr_v is not None else "\u2014"
                frc_str = f"{frc_v*100:.0f}%" if frc_v is not None else "\u2014"
                entry = (base_label, sl_str, dct_str, dpr_str, frc_str)
                if entry not in comp_groups[key]['entries']:
                    comp_groups[key]['entries'].append(entry)

        imgs: List[str] = []

        for data in comp_groups.values():
            atoms = data['atoms']
            entries = data['entries']
            comp_metrics = data.get('comp_metrics')
            n = len(entries)

            lines: List[str] = []
            for grp, sls, dct, dpr, frc in entries:
                line = f"\u2022 {grp}"
                if sls:
                    line += f" | {sls}"
                line += f"  dct={dct}  dpr={dpr}  chain={frc}"
                lines.append(line)
            if comp_metrics:
                sz, br, ecc, diam = comp_metrics
                lines.append(f"  size={sz}  br={br}  ecc={ecc}  diam={diam}")
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

    @load_HalogenGroups()
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
        # Backward-compat aliases (deprecated)
        count_mode: Optional[str] = None,
        graph_metrics: Optional[List[str]] = None,
        progress: bool = False,
        **kwargs,
    ) -> 'ResultsFingerprint':
        """Convert this molecule result to a ResultsFingerprint.

        Parameters
        ----------
        group_selection : str, default 'all'
            Which groups to include: 'all', 'oecd', 'generic', 'telomers',
            or 'generic+telomers'.
        component_metrics : list of str, default ['binary']
            Ordered list of per-component metrics.  Each entry produces one
            block of n_groups columns.  Valid values are count modes
            ('binary', 'count', 'max_component', 'total_component') and
            per-group graph metric names ('effective_graph_resistance',
            'min_dist_to_barycenter', 'branching', ...).
        selected_group_ids : list of int, optional
            Explicit group IDs (overrides group_selection).
        halogens : str or list of str, default 'F'
            Which halogen(s) to match SMARTS against.
        saturation : str or None, default 'per'
            Saturation filter: 'per', 'poly', or None.
        molecule_metrics : list of str, optional
            Molecule-wide scalar metrics appended as final columns.
        pfas_groups : list, optional
            Custom group objects (injected by decorator).
        preset : str, optional
            Named configuration overriding component_metrics and molecule_metrics.

        Returns
        -------
        ResultsFingerprint
        """
        from .fingerprints import PFASFingerprint
        return PFASFingerprint(
            self,
            preset=preset,
            component_metrics=component_metrics,
            group_selection=group_selection,
            selected_group_ids=selected_group_ids,
            halogens=halogens,
            saturation=saturation,
            molecule_metrics=molecule_metrics,
            count_mode=count_mode,
            graph_metrics=graph_metrics,
            progress=progress,
        )

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


class ResultsModel(list):
    """List-like container for HalogenGroups results.

    This subclasses ``list`` of molecule result dicts, so existing code
    that iterates over a list of results continues to work. It adds
    convenience methods for navigating and visualising components.
    """

    def __init__(self, iterable: Iterable[Dict[str, Any]] = ()):  # type: ignore[override]
        super().__init__(MoleculeResult(m) if not isinstance(m, MoleculeResult) else m for m in iterable)

    @classmethod
    def from_raw(cls, results: Iterable[Dict[str, Any]]) -> "ResultsModel":
        """Wrap an existing list of result dicts without changing them."""

        return cls(results)

    # --- Navigation helpers -------------------------------------------------

    def iter_group_matches(
        self,
        group_id: Optional[int] = None,
        group_name: Optional[str] = None,
    ) -> Iterator[Tuple[MoleculeResult, MatchView]]:
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
                        br_str   = f"{br_v:.2f}"  if br_v  is not None else "\u2014"
                        ecc_str  = f"{ecc_v:.2f}" if ecc_v is not None else "\u2014"
                        diam_str = (f"{float(diam_v):.0f}" if diam_v is not None and not _math.isnan(float(diam_v)) else "\u2014")
                        comp_groups[key] = {
                            'atoms': sorted(atoms),
                            'colour': colour,
                            'entries': [],
                            'comp_metrics': (str(comp.size), br_str, ecc_str, diam_str),
                        }
                    dct_v = comp.min_dist_to_center
                    dpr_v = comp.max_dist_to_periphery
                    frc_v = comp.component_fraction
                    dct_str = str(dct_v) if dct_v is not None else "\u2014"
                    dpr_str = str(dpr_v) if dpr_v is not None else "\u2014"
                    frc_str = f"{frc_v*100:.0f}%" if frc_v is not None else "\u2014"
                    entry = (base_label, sl_str, dct_str, dpr_str, frc_str)
                    if entry not in comp_groups[key]['entries']:
                        comp_groups[key]['entries'].append(entry)

            for data in comp_groups.values():
                atoms = data['atoms']
                colour = data['colour']
                entries = data['entries']
                comp_metrics = data.get('comp_metrics')

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
                    mol_img, entries, comp_metrics=comp_metrics, mol_label=f"mol#{mol_index + 1}"
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
        group_counts: Dict[Tuple[int, str], int] = {}
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
                        br_str   = f"{br_v:.2f}"  if br_v  is not None else "\u2014"
                        ecc_str  = f"{ecc_v:.2f}" if ecc_v is not None else "\u2014"
                        diam_str = (f"{float(diam_v):.0f}" if diam_v is not None and not _math.isnan(float(diam_v)) else "\u2014")
                        comp_groups[key] = {
                            'atoms': sorted(atoms),
                            'entries': [],
                            'comp_metrics': (str(comp.size), br_str, ecc_str, diam_str),
                        }
                    dct_v = comp.min_dist_to_center
                    dpr_v = comp.max_dist_to_periphery
                    frc_v = comp.component_fraction
                    dct_str = str(dct_v) if dct_v is not None else "\u2014"
                    dpr_str = str(dpr_v) if dpr_v is not None else "\u2014"
                    frc_str = f"{frc_v*100:.0f}%" if frc_v is not None else "\u2014"
                    entry = (base_label, sl_str, dct_str, dpr_str, frc_str)
                    if entry not in comp_groups[key]['entries']:
                        comp_groups[key]['entries'].append(entry)

            for data in comp_groups.values():
                atoms = data['atoms']
                entries = data['entries']
                comp_metrics = data.get('comp_metrics')
                n = len(entries)

                lines: List[str] = [f"mol#{mol_index + 1}"]
                for grp, sls, dct, dpr, frc in entries:
                    line = f"\u2022 {grp}"
                    if sls:
                        line += f" | {sls}"
                    line += f"  dct={dct}  dpr={dpr}  chain={frc}"
                    lines.append(line)
                if comp_metrics:
                    sz, br, ecc, diam = comp_metrics
                    lines.append(f"  size={sz}  br={br}  ecc={ecc}  diam={diam}")
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
        lines.append(_ansi("ResultsModel summary", _ANSI_BOLD))
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
    @load_HalogenGroups()
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
        # Backward-compat aliases (deprecated)
        count_mode: Optional[str] = None,
        graph_metrics: Optional[List[str]] = None,
        progress: bool = False,
        **kwargs,
    ) -> 'ResultsFingerprint':
        """Convert ResultsModel to ResultsFingerprint.

        Parameters
        ----------
        group_selection : str, default 'all'
            Which groups to include: 'all', 'oecd', 'generic', 'telomers',
            or 'generic+telomers'.
        component_metrics : list of str, default ['binary']
            Ordered list of per-component metrics.  Each entry produces one
            block of n_groups columns.  Valid values are count modes
            ('binary', 'count', 'max_component', 'total_component') and
            per-group graph metric names ('effective_graph_resistance',
            'min_dist_to_barycenter', 'branching', ...).
            Multiple halogens produce stacked blocks with ``[F]``, ``[Cl]`` suffixes.
        selected_group_ids : list of int, optional
            Explicit group IDs (overrides group_selection).
        halogens : str or list of str, default 'F'
            Which halogen(s) to match SMARTS against.
        saturation : str or None, default 'per'
            'per', 'poly', or None.
        molecule_metrics : list of str, optional
            Molecule-wide scalar metrics appended after all component-metric columns.
            E.g. ['n_components', 'mean_branching', 'max_diameter'].
        pfas_groups : list, optional
            Custom group objects (injected by decorator).
        preset : str, optional
            Named configuration overriding component_metrics and molecule_metrics.

        Returns
        -------
        ResultsFingerprint
        """
        from .fingerprints import PFASFingerprint
        return PFASFingerprint(
            self,
            preset=preset,
            component_metrics=component_metrics,
            group_selection=group_selection,
            selected_group_ids=selected_group_ids,
            halogens=halogens,
            saturation=saturation,
            molecule_metrics=molecule_metrics,
            count_mode=count_mode,
            graph_metrics=graph_metrics,
            progress=progress,
        )

    @classmethod
    def from_sql(
        cls,
        conn: Optional[Union[str, 'sqlalchemy.engine.Engine']] = None,
        filename: Optional[str] = None,
        components_table: str = "components",
        groups_table: str = "pfas_groups_in_compound",
        limit: Optional[int] = None,
    ) -> "ResultsModel":
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


# ResultsFingerprint is now PFASFingerprint — kept as alias.
from .fingerprints import PFASFingerprint as ResultsFingerprint
