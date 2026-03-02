"""HomologueSeries – result container for :func:`~PFASGroups.generate_homologues`.

Backward-compatible with the raw ``{InChIKey: {formula: mol}}`` dict that
``generate_homologues`` previously returned.  All existing code that iterates
or indexes the dict continues to work unchanged.

The class adds:

- :meth:`summary` / :meth:`summarise` – coloured text summary
- :meth:`show` / :meth:`plot` – grid image of all homologues
- :meth:`svg` – export to an SVG file
- :meth:`to_sql` – persist to SQLite or PostgreSQL

Example
-------
>>> from rdkit import Chem
>>> from HalogenGroups import generate_homologues
>>> pfoa = Chem.MolFromSmiles('OC(=O)' + 'C(F)(F)' * 7 + 'F')
>>> series = generate_homologues(pfoa)
>>> series.summary()
>>> img = series.show()
"""

from __future__ import annotations

from dataclasses import dataclass
from io import BytesIO
from typing import Any, List, Optional, Tuple, Union

from PIL import Image
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

# Re-use ANSI helpers and grid utility from results_model
from .results_model import _ansi, _ANSI_BOLD, _ANSI_HALOGEN, _grid_images

Color = Tuple[float, float, float]

# Colours for parent vs homologue panels
_PARENT_COLOR: Color = (0.95, 0.50, 0.10)   # orange
_HOMOLOGUE_COLOR: Color = (0.20, 0.60, 0.95) # blue


# ---------------------------------------------------------------------------
# HomologueEntry – lightweight descriptor for a single homologue
# ---------------------------------------------------------------------------

@dataclass
class HomologueEntry:
    """One homologue molecule from a :class:`HomologueSeries`.

    Parameters
    ----------
    mol :
        RDKit molecule.
    formula :
        Molecular formula string (e.g. ``'C6HF11O2'``).
    inchikey :
        InChIKey of this homologue.
    n_removed :
        Number of CX2 units removed relative to the parent.
    """

    mol: Chem.Mol
    formula: str
    inchikey: str
    n_removed: int

    @property
    def smiles(self) -> str:
        return Chem.MolToSmiles(self.mol)

    @property
    def n_carbons(self) -> int:
        return sum(1 for a in self.mol.GetAtoms() if a.GetSymbol() == 'C')

    def __repr__(self) -> str:  # pragma: no cover
        return (
            f"HomologueEntry(formula={self.formula!r}, "
            f"n_removed={self.n_removed}, smiles={self.smiles!r})"
        )


# ---------------------------------------------------------------------------
# HomologueSeries – dict subclass (backward-compatible) with rich API
# ---------------------------------------------------------------------------

class HomologueSeries(dict):
    """Result container for :func:`~PFASGroups.generate_homologues`.

    Inherits from :class:`dict` so that ``series[inchikey][formula]`` access
    and ``for inchikey, inner in series.items()`` still work exactly as before.

    Additional attributes set by :func:`~PFASGroups.generate_homologues`:

    Attributes
    ----------
    parent_mol : rdkit.Chem.Mol or None
        The original molecule that was shortened.
    parent_smiles : str or None
        Canonical SMILES of the parent.
    parent_formula : str or None
        Molecular formula of the parent.
    halogen : str
        Halogen element symbol used (``'F'``, ``'Cl'``, …).
    component_name : str
        Name of the component SMARTS used (e.g. ``'Perfluoroalkyl'``).

    These are set post-construction via :meth:`_set_metadata` (called
    internally by ``generate_homologues``).
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Metadata – populated by generate_homologues after construction
        self.parent_mol: Optional[Chem.Mol] = None
        self.parent_smiles: Optional[str] = None
        self.parent_formula: Optional[str] = None
        self.halogen: str = 'F'
        self.component_name: str = 'Perfluoroalkyl'

    def _set_metadata(
        self,
        parent_mol: Chem.Mol,
        halogen: str,
        component_name: str,
    ) -> None:
        """Called by ``generate_homologues`` to attach provenance info."""
        self.parent_mol = parent_mol
        self.parent_smiles = Chem.MolToSmiles(parent_mol)
        self.parent_formula = CalcMolFormula(parent_mol)
        self.halogen = halogen
        self.component_name = component_name

    # ------------------------------------------------------------------
    # Convenience iteration helpers
    # ------------------------------------------------------------------

    def entries(self) -> List[HomologueEntry]:
        """Return all homologues as a sorted list of :class:`HomologueEntry`.

        Sorted by the number of CX2 units removed (ascending), i.e. from
        the longest (closest to parent) to the shortest homologue.
        """
        parent_c = (
            sum(1 for a in self.parent_mol.GetAtoms() if a.GetSymbol() == 'C')
            if self.parent_mol is not None else None
        )
        result: List[HomologueEntry] = []
        for inchikey, inner in self.items():
            for formula, mol in inner.items():
                child_c = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == 'C')
                n_removed = (parent_c - child_c) if parent_c is not None else 0
                result.append(
                    HomologueEntry(
                        mol=mol,
                        formula=formula,
                        inchikey=inchikey,
                        n_removed=n_removed,
                    )
                )
        result.sort(key=lambda e: e.n_removed)
        return result

    def mols(self) -> List[Chem.Mol]:
        """Return all homologue molecules in order (shortest first)."""
        return [e.mol for e in self.entries()]

    # ------------------------------------------------------------------
    # Text output
    # ------------------------------------------------------------------

    def summarise(self) -> str:
        """Return a coloured text summary of the homologous series."""
        hal_code = _ANSI_HALOGEN.get(self.halogen, '')
        lines: List[str] = []
        lines.append(_ansi("HomologueSeries summary", _ANSI_BOLD))
        if self.parent_smiles:
            lines.append(
                f"- Parent  : {self.parent_smiles}  "
                f"({_ansi(self.parent_formula or '?', hal_code)})"
            )
        lines.append(
            f"- Halogen : {_ansi(self.halogen, hal_code)}"
            f"  |  Component : {self.component_name}"
        )
        if not self:
            lines.append("- No homologues found (no repeating CX2 units).")
            return "\n".join(lines)

        lines.append(f"- Homologues ({len(self)}):")
        for entry in self.entries():
            line = (
                f"    [{entry.n_removed:+d} CX2]  "
                f"{_ansi(entry.formula, hal_code)}  "
                f"{entry.smiles}"
            )
            lines.append(line)
        return "\n".join(lines)

    def summary(self) -> None:
        """Print the coloured text summary."""
        print(self.summarise())

    # ------------------------------------------------------------------
    # Visualisation helpers
    # ------------------------------------------------------------------

    def _draw_mol_image(
        self,
        mol: Chem.Mol,
        legend: str = "",
        subwidth: int = 300,
        subheight: int = 300,
        highlight_atoms: Optional[List[int]] = None,
        highlight_color: Optional[Color] = None,
    ) -> Image.Image:
        """Render *mol* as a PIL image using Cairo."""
        d2d = Draw.MolDraw2DCairo(subwidth, subheight)
        dopts = d2d.drawOptions()
        dopts.useBWAtomPalette()
        dopts.fixedBondLength = 20
        dopts.addAtomIndices = False
        dopts.maxFontSize = 14
        dopts.minFontSize = 11

        hl_atoms = highlight_atoms or []
        if hl_atoms and highlight_color:
            atom_colors = {idx: highlight_color for idx in hl_atoms}
            d2d.DrawMolecule(
                mol,
                legend=legend,
                highlightAtoms=hl_atoms,
                highlightAtomColors=atom_colors,
            )
        else:
            d2d.DrawMolecule(mol, legend=legend, highlightAtoms=hl_atoms)

        d2d.FinishDrawing()
        return Image.open(BytesIO(d2d.GetDrawingText()))

    def _draw_mol_svg(
        self,
        mol: Chem.Mol,
        legend: str = "",
        subwidth: int = 300,
        subheight: int = 300,
    ) -> str:
        """Render *mol* as an SVG string."""
        d2d = Draw.MolDraw2DSVG(subwidth, subheight)
        dopts = d2d.drawOptions()
        dopts.useBWAtomPalette()
        dopts.fixedBondLength = 20
        dopts.addAtomIndices = False
        dopts.maxFontSize = 14
        dopts.minFontSize = 11
        d2d.DrawMolecule(mol, legend=legend)
        d2d.FinishDrawing()
        return d2d.GetDrawingText()

    def show(
        self,
        display: bool = True,
        subwidth: int = 300,
        subheight: int = 300,
        ncols: int = 4,
        show_parent: bool = True,
    ) -> Image.Image:
        """Render all homologues as a grid image.

        Parameters
        ----------
        display :
            Whether to call ``image.show()`` immediately.
        subwidth, subheight :
            Pixel dimensions of each panel.
        ncols :
            Number of columns in the grid.
        show_parent :
            If ``True``, prepend a panel showing the parent molecule
            (highlighted in orange) before the homologues.

        Returns
        -------
        PIL.Image.Image
            Grid image.
        """
        if not self and not (show_parent and self.parent_mol is not None):
            raise ValueError("No molecules to display (empty HomologueSeries).")

        imgs: List[Image.Image] = []

        if show_parent and self.parent_mol is not None:
            parent_label = f"Parent  {self.parent_formula or ''}"
            img = self._draw_mol_image(
                self.parent_mol,
                legend=parent_label,
                subwidth=subwidth,
                subheight=subheight,
            )
            imgs.append(img)

        for entry in self.entries():
            label = f"-{entry.n_removed} CX2  {entry.formula}"
            img = self._draw_mol_image(
                entry.mol,
                legend=label,
                subwidth=subwidth,
                subheight=subheight,
            )
            imgs.append(img)

        if not imgs:
            raise ValueError("No molecules to display (empty HomologueSeries).")

        grid, _, _ = _grid_images(imgs, buffer=4, ncols=ncols)
        if display:
            grid.show()
        return grid

    # Alias
    plot = show

    def svg(
        self,
        filename: str,
        subwidth: int = 300,
        subheight: int = 300,
        ncols: int = 4,
        show_parent: bool = True,
    ) -> str:
        """Export all homologues to an SVG file.

        Parameters
        ----------
        filename :
            Destination path (e.g. ``'series.svg'``).
        subwidth, subheight :
            Pixel dimensions of each panel.
        ncols :
            Number of columns.
        show_parent :
            Prepend a parent-molecule panel.

        Returns
        -------
        str
            The path written to (*filename*).
        """
        import svgutils.transform as sg
        from .draw_mols import merge_svg

        svg_strs: List[str] = []

        if show_parent and self.parent_mol is not None:
            parent_label = f"Parent  {self.parent_formula or ''}"
            svg_strs.append(
                self._draw_mol_svg(self.parent_mol, legend=parent_label,
                                   subwidth=subwidth, subheight=subheight)
            )

        for entry in self.entries():
            label = f"-{entry.n_removed} CX2  {entry.formula}"
            svg_strs.append(
                self._draw_mol_svg(entry.mol, legend=label,
                                   subwidth=subwidth, subheight=subheight)
            )

        if not svg_strs:
            raise ValueError("No molecules to export (empty HomologueSeries).")

        figs = [sg.fromstring(s) for s in svg_strs]
        grid, _, _ = merge_svg(figs, buffer=4, ncols=ncols)
        grid.save(filename)
        return filename

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------

    def to_sql(
        self,
        conn: Optional[Union[str, Any]] = None,
        filename: Optional[str] = None,
        table_name: str = "homologue_series",
        if_exists: str = "append",
    ) -> None:
        """Save the homologous series to a SQL database.

        Each row represents one homologue.  The parent molecule is stored in
        the ``parent_smiles`` column for grouping.

        Parameters
        ----------
        conn :
            SQLAlchemy engine, connection string, or ``None``.
        filename :
            SQLite database path (used when *conn* is ``None``).
        table_name :
            Name of the table to write to.
        if_exists :
            ``'append'``, ``'replace'``, or ``'fail'``.

        Raises
        ------
        ValueError
            If neither *conn* nor *filename* is supplied.
        ImportError
            If ``pandas`` or ``sqlalchemy`` are not installed.

        Examples
        --------
        >>> series.to_sql(filename='pfas_homologues.db')
        >>> series.to_sql(conn='postgresql://user:pass@localhost/pfasdb')
        """
        try:
            import pandas as pd
            import sqlalchemy
        except ImportError as exc:
            raise ImportError(
                "pandas and sqlalchemy are required for to_sql. "
                "Install with: pip install pandas sqlalchemy"
            ) from exc

        if conn is None and filename is None:
            raise ValueError("Either 'conn' or 'filename' must be provided.")

        if conn is not None:
            engine = (
                sqlalchemy.create_engine(conn)
                if isinstance(conn, str)
                else conn
            )
        else:
            engine = sqlalchemy.create_engine(f"sqlite:///{filename}")

        rows = []
        for entry in self.entries():
            rows.append(
                {
                    "parent_smiles":  self.parent_smiles or "",
                    "parent_formula": self.parent_formula or "",
                    "halogen":        self.halogen,
                    "component_name": self.component_name,
                    "inchikey":       entry.inchikey,
                    "formula":        entry.formula,
                    "smiles":         entry.smiles,
                    "n_removed":      entry.n_removed,
                    "n_carbons":      entry.n_carbons,
                }
            )

        if not rows:
            return  # nothing to write

        df = pd.DataFrame(rows)
        df.to_sql(table_name, engine, if_exists=if_exists, index=False)

    # ------------------------------------------------------------------
    # Dunder niceties
    # ------------------------------------------------------------------

    def __repr__(self) -> str:  # pragma: no cover
        parent = self.parent_smiles or "?"
        return (
            f"HomologueSeries(parent={parent!r}, "
            f"halogen={self.halogen!r}, "
            f"n_homologues={len(self)})"
        )
