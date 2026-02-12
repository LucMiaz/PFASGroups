from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple, Union
import os

from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image
from io import BytesIO


Color = Tuple[float, float, float]


# Simple color palette to distinguish PFAS groups in highlight plots
_GROUP_COLORS: List[Color] = [
    (0.90, 0.10, 0.10),  # red
    (0.10, 0.40, 0.90),  # blue
    (0.10, 0.70, 0.10),  # green
    (0.90, 0.60, 0.10),  # orange
    (0.60, 0.10, 0.70),  # purple
    (0.00, 0.70, 0.70),  # teal
]


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


class MatchView(dict):
    """Wrapper for a single match dict (PFAS group or definition).

    Behaves as a normal dict but provides helpers for component access.
    """

    @property
    def is_group(self) -> bool:
        return self.get("type") == "PFASgroup"

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
        
        This reconstructs the molecule from the stored SMILES with explicit hydrogens,
        ensuring atom indices match those used during component detection.
        """
        smiles_h = self.get("smiles_with_h")
        if smiles_h:
            return Chem.MolFromSmiles(smiles_h)
        # Fallback for backwards compatibility
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
        """Return a text summary of this molecule's results.

        The summary includes:
        - SMILES representation
        - counts of PFAS group and definition matches
        - total number of components across all group matches
        - list of matched PFAS groups
        """

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
        lines.append("MoleculeResult summary")
        lines.append(f"- SMILES: {self.smiles}")
        lines.append(f"- PFAS group matches: {total_group_matches}")
        lines.append(f"- PFAS definition matches: {total_definition_matches}")
        lines.append(f"- Total components: {total_components}")

        if group_counts:
            lines.append("- Matched PFAS groups:")
            for name, count in sorted(group_counts.items(), key=lambda kv: kv[1], reverse=True):
                lines.append(f"  * {name}: {count} match(es)")

        return "\n".join(lines)

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

    def show(
        self,
        display: bool = True,
        subwidth: int = 300,
        subheight: int = 300,
        ncols: int = 4,
    ) -> Image.Image:
        """Show all component combinations for this molecule in a grid plot.

        Each panel corresponds to one combination of:
        - matched PFAS group
        - SMARTS/component type
        - individual component instance

        Parameters
        ----------
        display : bool, default True
            Whether to display the image immediately.
        subwidth : int, default 300
            Width of each sub-image in pixels.
        subheight : int, default 300
            Height of each sub-image in pixels.
        ncols : int, default 4
            Number of columns in the grid.

        Returns
        -------
        PIL.Image.Image
            Grid image containing all component visualizations.
        """

        # Use the molecule with hydrogens that was used during component detection
        mol = self.mol_with_h
        if mol is None:
            # Fallback to reconstructing from SMILES if not available
            mol = Chem.MolFromSmiles(self.smiles)
            if mol is None:
                raise ValueError(f"Cannot parse SMILES: {self.smiles}")
            mol = Chem.AddHs(mol)

        imgs: List[Image.Image] = []

        for match in self.matches:
            if not match.is_group:
                continue
            base_label = match.group_name or match.get("match_id", "")
            for comp_idx, comp in enumerate(match.components, start=1):
                atoms = comp.atoms
                if not atoms:
                    continue
                smarts_label = comp.smarts_label or ""
                legend = f"{base_label} | {smarts_label} | comp#{comp_idx}"

                d2d = Draw.MolDraw2DCairo(subwidth, subheight)
                dopts = d2d.drawOptions()
                dopts.useBWAtomPalette()
                dopts.fixedBondLength = 20
                dopts.addAtomIndices = True
                dopts.addBondIndices = False
                dopts.maxFontSize = 16
                dopts.minFontSize = 13
                d2d.DrawMolecule(mol, legend=legend, highlightAtoms=atoms)
                d2d.FinishDrawing()
                png = d2d.GetDrawingText()
                imgs.append(Image.open(BytesIO(png)))

        if not imgs:
            raise ValueError("No PFAS group components found to display.")

        grid, _, _ = _grid_images(imgs, buffer=4, ncols=ncols)
        if display:
            grid.show()
        return grid

    def svg(
        self,
        filename: str,
        subwidth: int = 300,
        subheight: int = 300,
        ncols: int = 4,
    ) -> str:
        """Export all component combinations to an SVG file (vector graphics).

        Parameters
        ----------
        filename : str
            Path to the output SVG file.
        subwidth : int, default 300
            Width of each sub-image in pixels.
        subheight : int, default 300
            Height of each sub-image in pixels.
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

        imgs: List[str] = []

        for match in self.matches:
            if not match.is_group:
                continue
            base_label = match.group_name or match.get("match_id", "")
            for comp_idx, comp in enumerate(match.components, start=1):
                atoms = comp.atoms
                if not atoms:
                    continue
                smarts_label = comp.smarts_label or ""
                legend = f"{base_label} | {smarts_label} | comp#{comp_idx}"

                d2d = Draw.MolDraw2DSVG(subwidth, subheight)
                dopts = d2d.drawOptions()
                dopts.useBWAtomPalette()
                dopts.fixedBondLength = 20
                dopts.addAtomIndices = True
                dopts.addBondIndices = False
                dopts.maxFontSize = 16
                dopts.minFontSize = 13
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
        except ImportError:
            raise ImportError("pandas and sqlalchemy are required for to_sql. Install with: pip install pandas sqlalchemy")

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
                components_data.append({
                    'smiles': self.smiles,
                    'group_id': match.group_id,
                    'group_name': match.group_name,
                    'smarts_label': comp.smarts_label,
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
    """List-like container for PFASGroups results.

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
    ) -> Image.Image:
        """Draw a single molecule with highlighted atoms as a PIL image."""

        d2d = Draw.MolDraw2DCairo(subwidth, subheight)
        dopts = d2d.drawOptions()
        dopts.useBWAtomPalette()
        dopts.fixedBondLength = 20
        dopts.addAtomIndices = True
        dopts.addBondIndices = False
        dopts.maxFontSize = 16
        dopts.minFontSize = 13
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

            smiles = mol_res.smiles
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            # IMPORTANT: Add hydrogens to match atom numbering used in parse_groups_in_mol
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
        subwidth: int = 300,
        subheight: int = 300,
        ncols: int = 4,
    ) -> Image.Image:
        """Show all component combinations in a grid plot.

        Each panel corresponds to one combination of:
        - molecule
        - matched PFAS group
        - SMARTS/component type
        - individual component instance

        This is intended for detailed inspection of how PFASgroups and
        component SMARTSs decompose the molecules.
        """

        imgs: List[Image.Image] = []

        for mol_index, mol_res in enumerate(self):  # type: ignore[assignment]
            # Use the molecule with hydrogens from parsing
            mol = mol_res.mol_with_h
            if mol is None:
                # Fallback to SMILES
                smiles = mol_res.smiles
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    continue
                mol = Chem.AddHs(mol)

            for match in mol_res.matches:
                if not match.is_group:
                    continue
                base_label = match.group_name or match.get("match_id", "")
                for comp_idx, comp in enumerate(match.components, start=1):
                    atoms = comp.atoms
                    if not atoms:
                        continue
                    smarts_label = comp.smarts_label or ""
                    legend = f"{base_label} | {smarts_label} | mol#{mol_index+1} comp#{comp_idx}"
                    img = self._draw_single_molecule(
                        mol,
                        atoms,
                        legend=legend,
                        subwidth=subwidth,
                        subheight=subheight,
                    )
                    imgs.append(img)

        if not imgs:
            raise ValueError("No PFAS group components found to display.")

        grid, _, _ = _grid_images(imgs, buffer=4, ncols=ncols)
        if display:
            grid.show()
        return grid

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
        except ImportError:
            raise ImportError("pandas and sqlalchemy are required for to_sql. Install with: pip install pandas sqlalchemy")

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
                components_data.append({
                    'smiles': self.smiles,
                    'group_id': match.group_id,
                    'group_name': match.group_name,
                    'smarts_label': comp.smarts_label,
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
        subwidth: int = 300,
        subheight: int = 300,
        ncols: int = 4,
    ) -> str:
        """Export all component combinations to an SVG file (vector graphics).

        Parameters
        ----------
        filename : str
            Path to the output SVG file.
        subwidth : int, default 300
            Width of each sub-image in pixels.
        subheight : int, default 300
            Height of each sub-image in pixels.
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
                # Fallback to SMILES
                smiles = mol_res.smiles
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    continue
                mol = Chem.AddHs(mol)

            for match in mol_res.matches:
                if not match.is_group:
                    continue
                base_label = match.group_name or match.get("match_id", "")
                for comp_idx, comp in enumerate(match.components, start=1):
                    atoms = comp.atoms
                    if not atoms:
                        continue
                    smarts_label = comp.smarts_label or ""
                    legend = f"{base_label} | {smarts_label} | mol#{mol_index+1} comp#{comp_idx}"
                    
                    d2d = Draw.MolDraw2DSVG(subwidth, subheight)
                    dopts = d2d.drawOptions()
                    dopts.useBWAtomPalette()
                    dopts.fixedBondLength = 20
                    dopts.addAtomIndices = True
                    dopts.addBondIndices = False
                    dopts.maxFontSize = 16
                    dopts.minFontSize = 13
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
        """Return a more detailed text summary of the results.

        The summary includes:
        - number of molecules
        - counts of PFAS group and definition matches
        - total number of components across all group matches
        - the most frequent PFAS groups by number of matches
        """

        total_molecules = len(self)
        total_group_matches = 0
        total_definition_matches = 0
        total_components = 0
        group_counts: Dict[str, int] = {}

        for mol_res in self:  # type: ignore[assignment]
            for m in mol_res.matches:
                if m.is_group:
                    total_group_matches += 1
                    total_components += len(m.components)
                    name = m.group_name or str(m.get("match_id", ""))
                    group_counts[name] = group_counts.get(name, 0) + 1
                elif m.is_definition:
                    total_definition_matches += 1

        unique_groups = len(group_counts)

        lines: List[str] = []
        lines.append("ResultsModel summary")
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
                lines.append(f"  * {name}: {count} matches")

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
        """Export all molecule results to a SQL database.

        Can write to either SQLite (via filename) or PostgreSQL/MySQL (via connection parameters).
        This method efficiently batches all molecules into the database in a single operation.

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
        except ImportError:
            raise ImportError("pandas and sqlalchemy are required for to_sql. Install with: pip install pandas sqlalchemy")

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

        # Prepare components data for all molecules
        components_data = []
        for mol_res in self:  # type: ignore[assignment]
            for match in mol_res.matches:
                if not match.is_group:
                    continue
                for comp in match.components:
                    components_data.append({
                        'smiles': mol_res.smiles,
                        'group_id': match.group_id,
                        'group_name': match.group_name,
                        'smarts_label': comp.smarts_label,
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
