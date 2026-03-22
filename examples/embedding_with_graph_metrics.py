"""
Embedding with graph-based descriptors
=======================================

This script demonstrates how to build PFASGroups embeddings with:

  * ``effective_graph_resistance`` (EGR) — the best single graph metric for
    group discrimination (benchmark preset ``'best'``);
  * ``n_spacer`` — telomer CH₂ linker chain length;
  * ``ring_size`` — smallest ring size overlapping the matched component;
  * molecule-wide descriptors via ``molecule_metrics``;
  * multi-halogen (4 × N_G) embeddings with ``halogens=None``.

Run from the repository root::

    python examples/embedding_with_graph_metrics.py
"""

from __future__ import annotations

import numpy as np
from HalogenGroups import parse_smiles
from PFASGroups import EMBEDDING_PRESETS

# ---------------------------------------------------------------------------
# Example molecules
# ---------------------------------------------------------------------------

SMILES = [
    # ── Linear perfluoroalkyl ───────────────────────────────────────────────
    "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",   # PFOA
    "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                          # PFBA
    "C(C(C(C(C(C(C(C(F)(F)S(=O)(=O)O)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F",  # PFOS (C8)
    # ── Fluorotelomers (CH₂ spacers) ───────────────────────────────────────
    "OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                             # 4:2 FTOH (2 CH₂)
    "OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",              # 6:2 FTOH (2 CH₂)
    "OC(=O)CC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",         # 6:2 FTCA (1 CH₂)
    "OC(=O)CCC(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                       # 4:2 FTCA (2 CH₂)
    # ── Cyclic / aryl ──────────────────────────────────────────────────────
    "Fc1c(F)c(F)c(F)c(F)c1F",                                        # perfluorobenzene
    "FC1(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C1(F)F",                    # perfluorocyclohexane
    "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C1=CNN=C1",  # PFOAZ (pyrazole)
    # ── Polychlorinated ────────────────────────────────────────────────────
    "ClC(Cl)(Cl)C(Cl)(Cl)Cl",                                        # hexachloroethane
]

LABELS = [
    "PFOA", "PFBA", "PFOS",
    "4:2 FTOH", "6:2 FTOH", "6:2 FTCA", "4:2 FTCA",
    "perf-benzene", "perf-cyclohexane", "PFOAZ",
    "HCA (Cl)",
]


def section(title: str) -> None:
    print(f"\n{'─' * 60}")
    print(f"  {title}")
    print(f"{'─' * 60}")


# ── 1. Parse all molecules ──────────────────────────────────────────────────
section("1. Parse molecules (all halogens)")
results = parse_smiles(SMILES)
print(f"Parsed {len(results)} molecules.")


# ── 2. Default binary embedding (fluorine only) ─────────────────────────────
section("2. Binary embedding  — fluorine only (default)")
arr_bin = results.to_array()
cols_bin = results.column_names()
print(f"Shape: {arr_bin.shape}   (n_molecules × n_groups)")
print(f"Non-zero entries per molecule: {(arr_bin > 0).sum(axis=1).tolist()}")


# ── 3. Preset 'best' ─────────────────────────────────────────────────────────
section("3. Preset 'best'  — binary + effective_graph_resistance")
print("  Description:", EMBEDDING_PRESETS['best']['description'])
arr_best = results.to_array(preset='best')
cols_best = results.column_names(preset='best')
print(f"  Shape: {arr_best.shape}   (= n_groups × 2 metric blocks)")


# ── 4. EGR embedding via component_metrics ───────────────────────────────────
section("4. EGR via component_metrics")
arr_egr = results.to_array(
    component_metrics=['binary', 'effective_graph_resistance']
)
# Note: identical to preset='best' above when no molecule_metrics are added
print(f"  Shape: {arr_egr.shape}")


# ── 5. n_spacer — telomer CH₂ linker length ──────────────────────────────────
section("5. n_spacer  — telomer spacer length per group")
arr_ns = results.to_array(component_metrics=['n_spacer'])
cols_ns = results.column_names(component_metrics=['n_spacer'])

# Find columns that are non-zero for at least one molecule
nonzero_cols = np.where((arr_ns > 0).any(axis=0))[0]
print(f"  Groups with non-zero n_spacer:")
for idx in nonzero_cols:
    values = arr_ns[:, idx]
    nz = [(LABELS[i], float(values[i])) for i in range(len(LABELS)) if values[i] > 0]
    print(f"    {cols_ns[idx]}: {nz}")


# ── 6. ring_size — smallest ring overlapping the component ───────────────────
section("6. ring_size  — cyclic/aryl ring size per group")
arr_rs = results.to_array(component_metrics=['ring_size'])
cols_rs = results.column_names(component_metrics=['ring_size'])

nonzero_cols = np.where((arr_rs > 0).any(axis=0))[0]
print(f"  Groups with non-zero ring_size:")
for idx in nonzero_cols:
    values = arr_rs[:, idx]
    nz = [(LABELS[i], float(values[i])) for i in range(len(LABELS)) if values[i] > 0]
    print(f"    {cols_rs[idx]}: {nz}")


# ── 7. Combined embedding: EGR + n_spacer + ring_size ────────────────────────
section("7. Combined EGR + n_spacer + ring_size + molecule metrics")
MOL_METRICS = [
    'n_components', 'max_size', 'mean_branching', 'max_branching',
    'mean_component_fraction', 'max_component_fraction',
]
arr_full = results.to_array(
    component_metrics=['effective_graph_resistance', 'n_spacer', 'ring_size'],
    molecule_metrics=MOL_METRICS,
)
cols_full = results.column_names(
    component_metrics=['effective_graph_resistance', 'n_spacer', 'ring_size'],
    molecule_metrics=MOL_METRICS,
)
n_groups = arr_bin.shape[1]                  # N_G from binary embedding above
print(f"  component_metrics blocks : 3  ({n_groups} groups each)")
print(f"  molecule_metrics columns : {len(MOL_METRICS)}")
print(f"  Total shape              : {arr_full.shape}")
print(f"  Last columns (mol-wide)  : {cols_full[-len(MOL_METRICS):]}")


# ── 8. Multi-halogen: 4 × N_G columns (F + Cl + Br + I) ─────────────────────
section("8. Multi-halogen embedding  (one parse per halogen, then hstack)")
# Each halogen contributes one block of N_G × n_metrics columns.
# The simplest approach is to parse once per halogen and concatenate.
arrs_4x, cols_4x_all = [], []
for hal in ['F', 'Cl', 'Br', 'I']:
    r_hal = parse_smiles(SMILES, halogens=hal)
    arrs_4x.append(r_hal.to_array(
        component_metrics=['effective_graph_resistance']
    ))
arr_4x = np.hstack(arrs_4x)
print(f"  Shape: {arr_4x.shape}   (= 4 halogens × N_G × 1 metric block)")
print(f"  Formula: 4 × {arrs_4x[0].shape[1]} = {arr_4x.shape[1]} columns")


# ── 9. Single-molecule access via PFASEmbedding ──────────────────────────────
section("9. Single-molecule PFASEmbedding")
# parse_smiles returns a PFASEmbeddingSet; each element is a PFASEmbedding
fp = results[0]                          # PFASEmbedding (dict-like) for first mol
arr_1d = fp.to_array(preset='best')
print(f"  SMILES            : {fp.smiles}")
print(f"  Shape (1-D, 1 mol): {arr_1d.shape}")

# PFASEmbeddingSet.to_array() stacks all molecules into a 2-D matrix
arr_2d = results.to_array(preset='best')
print(f"  Shape (2-D, {len(results)} mols): {arr_2d.shape}")


# ── 10. Accessing raw n_spacer / ring_size per component ────────────────────
section("10. Raw per-component values (parse_groups_in_mol)")
from PFASGroups.parser import parse_groups_in_mol
from rdkit import Chem

for smi, label in zip(SMILES[3:7], LABELS[3:7]):  # telomers only
    mol = Chem.MolFromSmiles(smi)
    matches, *_ = parse_groups_in_mol(mol)
    for grp, _idx, _smarts, comps in matches:
        ns_vals = [c.get('n_spacer', 0) for c in comps]
        rs_vals = [c.get('ring_size', 0) for c in comps]
        if any(ns_vals) or any(rs_vals):
            print(f"  {label:12s} | group: {grp.name:35s} | "
                  f"n_spacer={ns_vals}  ring_size={rs_vals}")


print("\nDone.")
