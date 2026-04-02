"""
Molecule renderer — wraps PFASGroups drawing to produce SVG bytes
that can be displayed inside a QSvgWidget or converted to a QPixmap.
"""
from __future__ import annotations

from io import BytesIO
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass


def embedding_to_svg(embedding, width: int = 300, height: int = 300) -> bytes | None:
    """Render a single PFASEmbedding to SVG bytes with highlighted atoms.

    Returns None if rendering fails.
    """
    try:
        from PFASGroups.draw_mols import plot_HalogenGroups
        from rdkit import Chem

        smi = embedding.get("smiles") or embedding.get("smi") or ""
        if not smi:
            return None

        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None

        fig, w, h = plot_HalogenGroups(
            mol,
            embedding,
            svg=True,
            subwidth=width,
            subheight=height,
        )
        buf = BytesIO()
        fig.save(buf)
        return buf.getvalue()
    except Exception:
        return _plain_mol_svg(smi, width, height)


def embedding_set_to_svgs(embedding_set, width: int = 280, height: int = 280) -> list[bytes | None]:
    """Render every molecule in a PFASEmbeddingSet to a list of SVG bytes."""
    return [embedding_to_svg(emb, width, height) for emb in embedding_set]


def _plain_mol_svg(smiles: str, width: int = 300, height: int = 300) -> bytes | None:
    """Fallback: draw a plain molecule (no highlighting) as SVG."""
    try:
        from rdkit import Chem
        from rdkit.Chem.Draw import rdMolDraw2D

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText().encode()
    except Exception:
        return None


def smiles_to_svg(smiles: str, width: int = 300, height: int = 300) -> bytes | None:
    """Render a SMILES string to SVG bytes (no highlighting)."""
    return _plain_mol_svg(smiles, width, height)
