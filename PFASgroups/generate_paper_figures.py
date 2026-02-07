"""Generate publication figures for PFASgroups.

This small driver script creates example figures for the PFASgroups
manuscript and supplementary information using the existing
plotting and parsing utilities.

It produces two main outputs:

- A main-text example figure showing a representative PFAS molecule
  with its fluorinated component highlighted.
- A supplementary figure illustrating different fluorinated
  components (e.g., linear vs branched) side by side, and
  optionally a fluorotelomer example if benchmark data is present.

The script is designed to be run from a local clone of the
PFASgroups repository, alongside the overleaf_PFASgroups_article
project. By default, it writes PNG files directly into the
"imgs" directory of the Overleaf project so that the LaTeX
manuscript can include them without manual renaming.
"""

from __future__ import annotations

from pathlib import Path
from typing import List

from rdkit import Chem

from .draw_mols import plot_pfasgroups, draw_images


# ---------------------------------------------------------------------------
# Output locations
# ---------------------------------------------------------------------------

def _get_overleaf_imgs_dir() -> Path:
    """Return the expected imgs/ directory in the Overleaf project.

    The expected directory layout (as on the author's machine) is::

        .../git/
            PFASGroups/
            overleaf_PFASgroups_article/

    This function walks up from this file until it reaches the
    shared parent directory and then constructs the path to
    ``overleaf_PFASgroups_article/imgs``. If that directory does
    not exist, it is created.
    """

    # This file lives in .../PFASGroups/PFASgroups/generate_paper_figures.py
    this_file = Path(__file__).resolve()
    pfasgroups_repo = this_file.parents[1]  # .../PFASGroups
    workspace_root = pfasgroups_repo.parent  # .../git

    overleaf_root = workspace_root / "overleaf_PFASgroups_article"
    imgs_dir = overleaf_root / "imgs"
    imgs_dir.mkdir(parents=True, exist_ok=True)
    return imgs_dir


# ---------------------------------------------------------------------------
# Figure generation helpers
# ---------------------------------------------------------------------------


def generate_main_text_figure(output_dir: Path) -> Path:
    """Generate the main-text PFAS example figure.

    Uses a linear perfluoroalkyl carboxylic acid (PFOA-like) example
    that is already employed in the test suite and highlights the
    fluorinated component detected by PFASgroups.

    Parameters
    ----------
    output_dir : Path
        Directory where the PNG file should be written.

    Returns
    -------
    Path
        Path to the generated PNG file.
    """

    # PFOA-like linear perfluoroalkyl carboxylic acid (from tests)
    smiles = "C(=O)(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"

    fig, _, _ = plot_pfasgroups(
        smiles,
        svg=False,
        subwidth=400,
        subheight=300,
        ncols=1,
        addAtomIndices=False,
        addBondIndices=False,
        bycomponent=True,
        panel_labels=["A: PFCA example"],
    )

    output_path = output_dir / "pfasgroups_main_example.png"
    # fig is a PIL.Image when svg=False
    fig.save(output_path)
    return output_path


def _get_telomer_example_from_sdf() -> str | None:
    """Return a SMILES string for a fluorotelomer example if available.

    Tries to read the first valid molecule from the benchmark
    ``PubChem_fluorotelomers.sdf`` file. If the file is missing
    or no molecule can be read, returns ``None``.
    """

    # The benchmark data lives in the PFASGroups repository
    this_file = Path(__file__).resolve()
    pfasgroups_repo = this_file.parents[1]
    sdf_path = pfasgroups_repo / "benchmark" / "data" / "PubChem_fluorotelomers.sdf"

    if not sdf_path.exists():
        return None

    from rdkit import Chem as _Chem  # local import to avoid hard dependency at import time

    supplier = _Chem.SDMolSupplier(str(sdf_path))
    for mol in supplier:
        if mol is not None:
            return _Chem.MolToSmiles(mol)
    return None


def generate_si_figure(output_dir: Path) -> Path:
    """Generate the supplementary figure with multiple PFAS components.

    The default layout includes:

    - A linear perfluoroalkyl carboxylic acid (same as main figure).
    - A branched perfluoroalkyl compound.
    - Optionally, a fluorotelomer example if benchmark data are
      available (first valid entry from ``PubChem_fluorotelomers.sdf``).

    Parameters
    ----------
    output_dir : Path
        Directory where the PNG file should be written.

    Returns
    -------
    Path
        Path to the generated PNG file.
    """

    smiles_list: List[str] = []

    # Linear PFAS (PFOA-like)
    smiles_linear = "C(=O)(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
    smiles_list.append(smiles_linear)

    # Branched PFAS (from tests)
    smiles_branched = "C(C(F)(F)F)(C(F)(F)F)C(F)(F)C(F)(F)C(F)(F)C(=O)O"
    smiles_list.append(smiles_branched)

    # Optional fluorotelomer example from benchmark SDF
    telomer_smiles = _get_telomer_example_from_sdf()
    if telomer_smiles is not None:
        smiles_list.append(telomer_smiles)

    fig, _, _ = plot_pfasgroups(
        smiles_list,
        svg=False,
        subwidth=400,
        subheight=300,
        ncols=len(smiles_list),
        addAtomIndices=False,
        addBondIndices=False,
        bycomponent=True,
        panel_labels=["A: Linear PFAS", "B: Branched PFAS"]
        if len(smiles_list) == 2
        else ["A: Linear PFAS", "B: Branched PFAS", "C: Fluorotelomer"],
    )

    output_path = output_dir / "pfasgroups_SI_components.png"
    fig.save(output_path)
    return output_path


def generate_advanced_si_figure(output_dir: Path) -> Path:
    """Generate an advanced SI figure illustrating components and path types.

    Panels:
    - A: Molecule with two disjoint fluorinated components.
    - B: Same-type example highlighting only perfluoroalkyl paths.
    - C: Same-type example highlighting only polyfluoroalkyl paths.
    """

    # Molecule with two fluorinated components attached to an aromatic ring
    smiles_two_components = "FC(F)(F)C(F)(F)C1=CC(=CC=C1)C(F)(F)C(F)(F)F"

    fig_A, _, _ = plot_pfasgroups(
        [smiles_two_components],
        svg=False,
        subwidth=350,
        subheight=280,
        ncols=1,
        addAtomIndices=False,
        addBondIndices=False,
        bycomponent=True,
        panel_labels=["A: Two fluorinated components"],
    )

    # Mixed perfluoro/polyfluoro example (one CF2 and one CHF in the chain)
    smiles_mixed = "C(=O)(O)C(F)(F)C(F)(F)C(F)(H)F"

    fig_B, _, _ = plot_pfasgroups(
        [smiles_mixed],
        svg=False,
        subwidth=350,
        subheight=280,
        ncols=1,
        addAtomIndices=False,
        addBondIndices=False,
        bycomponent=True,
        paths=["Perfluoroalkyl"],
        panel_labels=["B: Perfluoroalkyl paths"],
    )

    fig_C, _, _ = plot_pfasgroups(
        [smiles_mixed],
        svg=False,
        subwidth=350,
        subheight=280,
        ncols=1,
        addAtomIndices=False,
        addBondIndices=False,
        bycomponent=True,
        paths=["Polyfluoroalkyl"],
        panel_labels=["C: Polyfluoroalkyl paths"],
    )

    # Combine the three panels into a single row
    combined_fig, _, _ = draw_images([fig_A, fig_B, fig_C], buffer=10, ncols=3, svg=False)

    output_path = output_dir / "pfasgroups_SI_components_advanced.png"
    combined_fig.save(output_path)
    return output_path


def main() -> None:
    """Entry point for script-style execution.

    When run as ``python -m PFASgroups.generate_paper_figures`` or
    ``python PFASgroups/generate_paper_figures.py``, this will generate
    both the main-text and SI figures in the Overleaf ``imgs``
    directory and print their locations.
    """

    imgs_dir = _get_overleaf_imgs_dir()

    main_fig = generate_main_text_figure(imgs_dir)
    si_fig = generate_si_figure(imgs_dir)
    si_adv_fig = generate_advanced_si_figure(imgs_dir)

    print("Generated main-text figure:", main_fig)
    print("Generated SI figure:", si_fig)
    print("Generated advanced SI figure:", si_adv_fig)


if __name__ == "__main__":  # pragma: no cover
    main()
