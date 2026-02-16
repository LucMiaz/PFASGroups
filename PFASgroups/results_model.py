from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple, Union
import os
import warnings

from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image
from io import BytesIO
import numpy as np
import pandas as pd
from scipy.stats import entropy


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

    def summary(self) -> None:
        """Print a detailed summary of matched groups and components.
        
        The summary includes:
        - List of matched group IDs and names
        - Number of components by SMARTS type for each group
        - Component sizes (number of atoms)
        """
        print("=" * 80)
        print(f"MOLECULE: {self.smiles}")
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
            print(f"Group {group_id}: {group_name}")
            
            # Group components by SMARTS type
            by_smarts: Dict[Optional[str], List[ComponentView]] = {}
            for comp in components:
                smarts = comp.smarts_label
                if smarts not in by_smarts:
                    by_smarts[smarts] = []
                by_smarts[smarts].append(comp)
            
            # Display components by SMARTS type
            for smarts_label, comps in sorted(by_smarts.items(), key=lambda x: x[0] or ""):
                if smarts_label:
                    print(f"  SMARTS: {smarts_label}")
                else:
                    print(f"  SMARTS: (no label)")
                
                print(f"    Components: {len(comps)}")
                
                # Show sizes of each component
                sizes = [len(comp.atoms) for comp in comps]
                if sizes:
                    sizes_str = ", ".join(map(str, sorted(sizes, reverse=True)))
                    print(f"    Sizes (atoms): {sizes_str}")
            
            print()
        
        print("=" * 80)

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

    def summary(self) -> None:
        """Print a detailed summary of matched groups and components across all molecules.
        
        The summary includes:
        - Total number of molecules
        - List of all matched group IDs and names
        - For each group: number of components by SMARTS type
        - Component size statistics (min, max, mean)
        """
        print("=" * 80)
        print(f"RESULTS SUMMARY: {len(self)} molecule(s)")
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
            print(f"Group {group_id}: {group_name}")
            
            # Group components by SMARTS type
            by_smarts: Dict[Optional[str], List[ComponentView]] = {}
            for comp in components:
                smarts = comp.smarts_label
                if smarts not in by_smarts:
                    by_smarts[smarts] = []
                by_smarts[smarts].append(comp)
            
            # Display components by SMARTS type
            for smarts_label, comps in sorted(by_smarts.items(), key=lambda x: x[0] or ""):
                if smarts_label:
                    print(f"  SMARTS: {smarts_label}")
                else:
                    print(f"  SMARTS: (no label)")
                
                print(f"    Total components: {len(comps)}")
                
                # Calculate size statistics
                sizes = [len(comp.atoms) for comp in comps]
                if sizes:
                    min_size = min(sizes)
                    max_size = max(sizes)
                    mean_size = sum(sizes) / len(sizes)
                    print(f"    Size range (atoms): {min_size} - {max_size} (mean: {mean_size:.1f})")
            
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

    def to_sql(
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
        >>> results.to_sql(conn=engine)
        >>> 
        >>> # Using filename (legacy)
        >>> results.to_sql(filename='pfas.db')
        """
        try:
            import pandas as pd
            import sqlalchemy
        except ImportError:
            raise ImportError("pandas and sqlalchemy are required for to_sql. Install with: pip install pandas sqlalchemy")

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

    def to_fingerprint(
        self,
        group_selection: str = 'all',
        count_mode: str = 'binary',
        selected_group_ids: Optional[List[int]] = None,
    ) -> 'ResultsFingerprint':
        """Convert ResultsModel to ResultsFingerprint for dimensionality reduction.

        Parameters
        ----------
        group_selection : str, default 'all'
            Which groups to include in fingerprint:
            - 'all': All 55 groups
            - 'generic': Generic groups only (groups 29-55)
            - 'oecd': OECD groups only (groups 1-28)
            - 'telomers': Telomer-related groups only
            - 'generic+telomers': Combination of generic and telomer groups
        count_mode : str, default 'binary'
            How to encode group matches:
            - 'binary': 1 if present, 0 if absent
            - 'count': Number of matches
            - 'max_component': Maximum component size
        selected_group_ids : list of int, optional
            Explicit list of group IDs to include (overrides group_selection)

        Returns
        -------
        ResultsFingerprint
            Fingerprint representation of the results
        """
        from .fingerprints import generate_fingerprint
        
        # Extract SMILES from results
        smiles_list = [mol_res.smiles for mol_res in self]  # type: ignore[assignment]
        
        # Define group selections
        if selected_group_ids is not None:
            selected_groups = selected_group_ids
        elif group_selection == 'all':
            selected_groups = None  # Use all groups
        elif group_selection == 'oecd':
            selected_groups = list(range(1, 29))  # OECD groups 1-28
        elif group_selection == 'generic':
            selected_groups = list(range(29, 56))  # Generic groups 29-55
        elif group_selection == 'telomers':
            # Telomer-related group IDs (example - adjust based on actual group definitions)
            selected_groups = [24, 25, 26]  # Adjust as needed
        elif group_selection == 'generic+telomers':
            selected_groups = list(range(29, 56)) + [24, 25, 26]  # Adjust as needed
        else:
            raise ValueError(
                f"Unknown group_selection: {group_selection}. "
                f"Choose from: 'all', 'generic', 'oecd', 'telomers', 'generic+telomers'"
            )
        
        # Generate fingerprints
        fps, info = generate_fingerprint(
            smiles_list,
            selected_groups=selected_groups,
            representation='vector',
            count_mode=count_mode,
        )
        
        return ResultsFingerprint(
            fingerprints=fps,
            smiles=smiles_list,
            group_names=info['group_names'],
            group_selection=group_selection,
            count_mode=count_mode,
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
        except ImportError:
            raise ImportError("sqlalchemy is required for SQL operations. Install with: pip install sqlalchemy")
        
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
                    'type': 'PFASgroup',
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


class ResultsFingerprint:
    """Container for PFAS group fingerprints with dimensionality reduction analysis.
    
    This class encapsulates fingerprint representations of PFASGroups results and
    provides methods for performing dimensionality reduction (PCA, kernel-PCA, t-SNE, UMAP)
    and comparison between fingerprint sets using KL divergence.
    
    Parameters
    ----------
    fingerprints : np.ndarray
        Fingerprint matrix (n_molecules, n_groups)
    smiles : list of str
        SMILES strings for each molecule
    group_names : list of str
        Names of the PFAS groups in fingerprint
    group_selection : str
        Type of group selection used
    count_mode : str
        Encoding mode used (binary, count, max_component)
    """
    
    def __init__(
        self,
        fingerprints: np.ndarray,
        smiles: List[str],
        group_names: List[str],
        group_selection: str = 'all',
        count_mode: str = 'binary',
    ):
        self.fingerprints = np.array(fingerprints)
        self.smiles = smiles
        self.group_names = group_names
        self.group_selection = group_selection
        self.count_mode = count_mode
        
        if len(self.fingerprints) != len(self.smiles):
            raise ValueError(
                f"Fingerprints and SMILES length mismatch: "
                f"{len(self.fingerprints)} != {len(self.smiles)}"
            )
    
    def __len__(self) -> int:
        return len(self.fingerprints)
    
    def __repr__(self) -> str:
        return (
            f"ResultsFingerprint(n_molecules={len(self)}, "
            f"n_groups={len(self.group_names)}, "
            f"group_selection='{self.group_selection}', "
            f"count_mode='{self.count_mode}')"
        )
    
    def summary(self) -> str:
        """Return text summary of fingerprint statistics."""
        lines = []
        lines.append("ResultsFingerprint Summary")
        lines.append("=" * 50)
        lines.append(f"Molecules: {len(self)}")
        lines.append(f"Groups: {len(self.group_names)}")
        lines.append(f"Group selection: {self.group_selection}")
        lines.append(f"Count mode: {self.count_mode}")
        lines.append(f"Fingerprint shape: {self.fingerprints.shape}")
        lines.append(f"Non-zero entries: {np.count_nonzero(self.fingerprints)}")
        lines.append(f"Sparsity: {1 - np.count_nonzero(self.fingerprints) / self.fingerprints.size:.2%}")
        
        # Group frequency statistics
        group_frequencies = np.sum(self.fingerprints > 0, axis=0)
        lines.append("\nMost common groups:")
        top_indices = np.argsort(group_frequencies)[::-1][:10]
        for idx in top_indices:
            if group_frequencies[idx] > 0:
                lines.append(f"  {self.group_names[idx]}: {group_frequencies[idx]} molecules "
                           f"({100 * group_frequencies[idx] / len(self):.1f}%)")
        
        return "\n".join(lines)
    
    def perform_pca(
        self,
        n_components: int = 2,
        plot: bool = True,
        output_file: Optional[str] = None,
    ) -> Dict[str, Any]:
        """Perform PCA analysis on fingerprints.
        
        Parameters
        ----------
        n_components : int, default 2
            Number of principal components
        plot : bool, default True
            Whether to create visualization
        output_file : str, optional
            Path to save plot (if None, displays interactively)
        
        Returns
        -------
        dict
            Dictionary containing:
            - 'transformed': PCA-transformed data
            - 'explained_variance': Explained variance ratio
            - 'components': Principal component vectors
            - 'pca_model': Fitted PCA model
        """
        try:
            from sklearn.decomposition import PCA
            from sklearn.preprocessing import StandardScaler
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError(
                "scikit-learn and matplotlib required. Install with: "
                "pip install scikit-learn matplotlib"
            )
        
        # Standardize features
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(self.fingerprints)
        
        # Perform PCA
        pca = PCA(n_components=n_components)
        X_pca = pca.fit_transform(X_scaled)
        
        # Create plot if requested
        if plot and n_components >= 2:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
            
            # Scatter plot
            scatter = ax1.scatter(X_pca[:, 0], X_pca[:, 1], alpha=0.6, s=50)
            ax1.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
            ax1.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
            ax1.set_title('PCA of PFAS Group Fingerprints')
            ax1.grid(True, alpha=0.3)
            
            # Scree plot
            ax2.bar(range(1, n_components + 1), pca.explained_variance_ratio_)
            ax2.set_xlabel('Principal Component')
            ax2.set_ylabel('Explained Variance Ratio')
            ax2.set_title('Scree Plot')
            ax2.grid(True, alpha=0.3)
            
            plt.tight_layout()
            
            if output_file:
                plt.savefig(output_file, dpi=300, bbox_inches='tight')
                print(f"PCA plot saved to {output_file}")
            else:
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
    ) -> Dict[str, Any]:
        """Perform kernel PCA analysis on fingerprints.
        
        Parameters
        ----------
        n_components : int, default 2
            Number of components
        kernel : str, default 'rbf'
            Kernel type: 'linear', 'poly', 'rbf', 'sigmoid', 'cosine'
        gamma : float, optional
            Kernel coefficient (if None, uses 1/n_features)
        plot : bool, default True
            Whether to create visualization
        output_file : str, optional
            Path to save plot
        
        Returns
        -------
        dict
            Dictionary containing transformed data and model
        """
        try:
            from sklearn.decomposition import KernelPCA
            from sklearn.preprocessing import StandardScaler
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError(
                "scikit-learn and matplotlib required. Install with: "
                "pip install scikit-learn matplotlib"
            )
        
        # Standardize features
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(self.fingerprints)
        
        # Perform kernel PCA
        if gamma is None:
            gamma = 1.0 / X_scaled.shape[1]
        
        kpca = KernelPCA(n_components=n_components, kernel=kernel, gamma=gamma)
        X_kpca = kpca.fit_transform(X_scaled)
        
        # Create plot if requested
        if plot and n_components >= 2:
            fig, ax = plt.subplots(figsize=(8, 6))
            
            scatter = ax.scatter(X_kpca[:, 0], X_kpca[:, 1], alpha=0.6, s=50)
            ax.set_xlabel('Kernel PC1')
            ax.set_ylabel('Kernel PC2')
            ax.set_title(f'Kernel PCA ({kernel} kernel) of PFAS Group Fingerprints')
            ax.grid(True, alpha=0.3)
            
            plt.tight_layout()
            
            if output_file:
                plt.savefig(output_file, dpi=300, bbox_inches='tight')
                print(f"Kernel PCA plot saved to {output_file}")
            else:
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
    ) -> Dict[str, Any]:
        """Perform t-SNE analysis on fingerprints.
        
        Parameters
        ----------
        n_components : int, default 2
            Number of dimensions for embedding
        perplexity : float, default 30.0
            t-SNE perplexity parameter (5-50 typical)
        learning_rate : float, default 200.0
            Learning rate for optimization
        max_iter : int, default 1000
            Maximum number of iterations
        plot : bool, default True
            Whether to create visualization
        output_file : str, optional
            Path to save plot
        
        Returns
        -------
        dict
            Dictionary containing transformed data and model
        """
        try:
            from sklearn.manifold import TSNE
            from sklearn.preprocessing import StandardScaler
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError(
                "scikit-learn and matplotlib required. Install with: "
                "pip install scikit-learn matplotlib"
            )
        
        # Standardize features
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(self.fingerprints)
        
        # Perform t-SNE
        tsne = TSNE(
            n_components=n_components,
            perplexity=perplexity,
            learning_rate=learning_rate,
            max_iter=max_iter,
            random_state=42,
        )
        X_tsne = tsne.fit_transform(X_scaled)
        
        # Create plot if requested
        if plot and n_components >= 2:
            fig, ax = plt.subplots(figsize=(8, 6))
            
            scatter = ax.scatter(X_tsne[:, 0], X_tsne[:, 1], alpha=0.6, s=50)
            ax.set_xlabel('t-SNE 1')
            ax.set_ylabel('t-SNE 2')
            ax.set_title(f't-SNE of PFAS Group Fingerprints (perplexity={perplexity})')
            ax.grid(True, alpha=0.3)
            
            plt.tight_layout()
            
            if output_file:
                plt.savefig(output_file, dpi=300, bbox_inches='tight')
                print(f"t-SNE plot saved to {output_file}")
            else:
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
    ) -> Dict[str, Any]:
        """Perform UMAP analysis on fingerprints.
        
        Parameters
        ----------
        n_components : int, default 2
            Number of dimensions for embedding
        n_neighbors : int, default 15
            Number of neighbors to consider
        min_dist : float, default 0.1
            Minimum distance between points
        metric : str, default 'euclidean'
            Distance metric
        plot : bool, default True
            Whether to create visualization
        output_file : str, optional
            Path to save plot
        
        Returns
        -------
        dict
            Dictionary containing transformed data and model
        """
        try:
            import umap
            from sklearn.preprocessing import StandardScaler
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError(
                "umap-learn and matplotlib required. Install with: "
                "pip install umap-learn matplotlib"
            )
        
        # Standardize features
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(self.fingerprints)
        
        # Perform UMAP
        reducer = umap.UMAP(
            n_components=n_components,
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            metric=metric,
            random_state=42,
        )
        X_umap = reducer.fit_transform(X_scaled)
        
        # Create plot if requested
        if plot and n_components >= 2:
            fig, ax = plt.subplots(figsize=(8, 6))
            
            scatter = ax.scatter(X_umap[:, 0], X_umap[:, 1], alpha=0.6, s=50)
            ax.set_xlabel('UMAP 1')
            ax.set_ylabel('UMAP 2')
            ax.set_title(f'UMAP of PFAS Group Fingerprints (n_neighbors={n_neighbors})')
            ax.grid(True, alpha=0.3)
            
            plt.tight_layout()
            
            if output_file:
                plt.savefig(output_file, dpi=300, bbox_inches='tight')
                print(f"UMAP plot saved to {output_file}")
            else:
                plt.show()
            
            plt.close()
        
        return {
            'transformed': X_umap,
            'umap_model': reducer,
            'scaler': scaler,
            'n_neighbors': n_neighbors,
            'min_dist': min_dist,
        }
    
    def compare_kld(self, other: 'ResultsFingerprint', method: str = 'minmax') -> float:
        """Compare two fingerprint sets using KL divergence.
        
        This implements the minmaxKLd method for comparing distributions of
        PFAS group fingerprints between two datasets.
        
        Parameters
        ----------
        other : ResultsFingerprint
            Other fingerprint set to compare against
        method : str, default 'minmax'
            Comparison method:
            - 'minmax': Min-max normalized KL divergence (symmetric)
            - 'forward': KL(self || other)
            - 'reverse': KL(other || self)
            - 'symmetric': (KL(self || other) + KL(other || self)) / 2
        
        Returns
        -------
        float
            KL divergence value (lower = more similar)
        """
        if len(self.group_names) != len(other.group_names):
            raise ValueError(
                f"Fingerprints must have same number of groups: "
                f"{len(self.group_names)} != {len(other.group_names)}"
            )
        
        # Compute frequency distributions for each group
        p_dist = self._compute_group_distribution()
        q_dist = other._compute_group_distribution()
        
        # Add small constant to avoid log(0)
        epsilon = 1e-10
        p_dist = p_dist + epsilon
        q_dist = q_dist + epsilon
        
        # Normalize
        p_dist = p_dist / p_dist.sum()
        q_dist = q_dist / q_dist.sum()
        
        if method == 'forward':
            return entropy(p_dist, q_dist)
        elif method == 'reverse':
            return entropy(q_dist, p_dist)
        elif method == 'symmetric':
            return (entropy(p_dist, q_dist) + entropy(q_dist, p_dist)) / 2
        elif method == 'minmax':
            # Min-max normalized symmetric KL divergence
            kl_forward = entropy(p_dist, q_dist)
            kl_reverse = entropy(q_dist, p_dist)
            kl_sym = (kl_forward + kl_reverse) / 2
            
            # Normalize by maximum possible KL divergence
            # Max KL occurs when distributions are maximally different
            max_kl = np.log(len(p_dist))  # Theoretical maximum
            
            return kl_sym / max_kl if max_kl > 0 else 0.0
        else:
            raise ValueError(f"Unknown method: {method}")
    
    def _compute_group_distribution(self) -> np.ndarray:
        """Compute frequency distribution of groups across molecules."""
        # Sum across molecules to get total count per group
        group_counts = np.sum(self.fingerprints > 0, axis=0)
        return group_counts.astype(float)
    
    def to_sql(
        self,
        conn: Optional[Union[str, 'sqlalchemy.engine.Engine']] = None,
        filename: Optional[str] = None,
        table_name: str = "fingerprints",
        metadata_table: str = "fingerprint_metadata",
        if_exists: str = "append",
    ) -> None:
        """Save fingerprints to SQL database.
        
        Parameters
        ----------
        conn : str or SQLAlchemy Engine, optional
            Database connection string or engine
        filename : str, optional
            SQLite database filename (alternative to conn)
        table_name : str, default "fingerprints"
            Name of the table for fingerprint data
        metadata_table : str, default "fingerprint_metadata"
            Name of the table for metadata
        if_exists : str, default "append"
            How to behave if table exists: 'fail', 'replace', 'append'
        """
        try:
            import sqlalchemy
        except ImportError:
            raise ImportError("sqlalchemy required. Install with: pip install sqlalchemy")
        
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
        
        # Prepare fingerprint data
        data = []
        for i, smiles in enumerate(self.smiles):
            for j, group_name in enumerate(self.group_names):
                value = self.fingerprints[i, j]
                if value > 0:  # Only store non-zero values for efficiency
                    data.append({
                        'smiles': smiles,
                        'group_name': group_name,
                        'group_index': j,
                        'value': float(value),
                    })
        
        df_fingerprints = pd.DataFrame(data)
        df_fingerprints.to_sql(table_name, engine, if_exists=if_exists, index=False)
        
        # Save metadata
        metadata = pd.DataFrame([{
            'group_selection': self.group_selection,
            'count_mode': self.count_mode,
            'n_molecules': len(self),
            'n_groups': len(self.group_names),
        }])
        metadata.to_sql(metadata_table, engine, if_exists=if_exists, index=False)
        
        print(f"Saved {len(df_fingerprints)} fingerprint entries to {table_name}")
        print(f"Saved metadata to {metadata_table}")
    
    @classmethod
    def from_sql(
        cls,
        conn: Optional[Union[str, 'sqlalchemy.engine.Engine']] = None,
        filename: Optional[str] = None,
        table_name: str = "fingerprints",
        metadata_table: str = "fingerprint_metadata",
        limit: Optional[int] = None,
    ) -> 'ResultsFingerprint':
        """Load fingerprints from SQL database.
        
        Parameters
        ----------
        conn : str or SQLAlchemy Engine, optional
            Database connection string or engine
        filename : str, optional
            SQLite database filename
        table_name : str, default "fingerprints"
            Name of the fingerprints table
        metadata_table : str, default "fingerprint_metadata"
            Name of the metadata table
        limit : int, optional
            Limit number of molecules to load
        
        Returns
        -------
        ResultsFingerprint
            Loaded fingerprints
        """
        try:
            import sqlalchemy
        except ImportError:
            raise ImportError("sqlalchemy required. Install with: pip install sqlalchemy")
        
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
        
        # Load metadata
        df_metadata = pd.read_sql(f"SELECT * FROM {metadata_table} LIMIT 1", engine)
        group_selection = df_metadata['group_selection'].iloc[0]
        count_mode = df_metadata['count_mode'].iloc[0]
        
        # Load fingerprint data
        query = f"SELECT * FROM {table_name}"
        if limit is not None:
            query += f" LIMIT {limit}"
        
        df_fingerprints = pd.read_sql(query, engine)
        
        # Reconstruct fingerprint matrix
        smiles_list = sorted(df_fingerprints['smiles'].unique())
        group_names = sorted(df_fingerprints['group_name'].unique())
        
        fingerprints = np.zeros((len(smiles_list), len(group_names)))
        
        for _, row in df_fingerprints.iterrows():
            i = smiles_list.index(row['smiles'])
            j = group_names.index(row['group_name'])
            fingerprints[i, j] = row['value']
        
        return cls(
            fingerprints=fingerprints,
            smiles=smiles_list,
            group_names=group_names,
            group_selection=group_selection,
            count_mode=count_mode,
        )

