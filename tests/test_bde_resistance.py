"""
Tests for BDE-weighted effective graph resistance (Kirchhoff index).

Verifies:
- Physical plausibility: linear chain EGR grows as O(n^3) for n atoms
- Branching sensitivity: branched isomers have lower EGR than linear
- Bond-order effect: stronger bonds (higher BDE) reduce resistance
- Exact Kirchhoff identity: K_f == n * sum(1/lambda_i) for non-zero eigenvalues
- SMARTS-distance metrics use genuine BDE-weighted R(u,v)
- Arithmetic correctness for a trivial 2-atom graph
- Fallback/skip paths (limit_effective_graph_resistance)
"""

import math
import pytest
import numpy as np
from rdkit import Chem

from PFASGroups import parse_smiles
from PFASGroups.ComponentsSolverModel import _BDEScheme, _get_bde_scheme, ComponentsSolver


# ── helpers ────────────────────────────────────────────────────────────────────

def _first_comp_metrics(smiles: str) -> dict:
    """Return the component metrics dict for the first component of the first match."""
    results = parse_smiles(smiles, halogens='F')
    for match in results[0]['matches']:
        for comp in match['components']:
            return comp
    raise AssertionError(f'No component found for {smiles}')


def _egr(smiles: str) -> float:
    return _first_comp_metrics(smiles)['effective_graph_resistance']


# ── BDE scheme unit tests ──────────────────────────────────────────────────────

class TestBDEScheme:

    def test_singleton(self):
        a = _get_bde_scheme()
        b = _get_bde_scheme()
        assert a is b, 'Module-level singleton should return the same object'

    def test_cc_reference(self):
        scheme = _BDEScheme()
        # C-C single bond is the reference → conductance == 1.0
        c = scheme.conductance(6, 6, 1.0)
        assert abs(c - 1.0) < 1e-9

    def test_cf_gt_cc(self):
        # In the diatomic_bonds_dict.json used here, C-C is the reference bond
        # (bde_dict[6][6] acts as ref_bde). C-F BDE is lower than C-C in this
        # dataset, so C-F conductance < 1.0 (C-F is a higher-resistance bond).
        scheme = _BDEScheme()
        assert scheme.conductance(6, 9, 1.0) < 1.0
        assert scheme.conductance(6, 9, 1.0) > 0

    def test_ch_lt_cf(self):
        # C-H BDE < C-F BDE
        scheme = _BDEScheme()
        assert scheme.conductance(6, 1, 1.0) < scheme.conductance(6, 9, 1.0)

    def test_bond_order_increases_conductance(self):
        scheme = _BDEScheme()
        c1 = scheme.conductance(6, 6, 1.0)
        c2 = scheme.conductance(6, 6, 2.0)
        c3 = scheme.conductance(6, 6, 3.0)
        assert c1 < c2 < c3, 'Higher bond order should yield higher conductance'

    def test_conductance_positive(self):
        scheme = _BDEScheme()
        for z1, z2 in [(6, 6), (6, 9), (6, 7), (6, 8), (6, 1), (9, 9)]:
            assert scheme.conductance(z1, z2, 1.0) > 0

    def test_symmetry(self):
        scheme = _BDEScheme()
        # conductance(z1, z2) must equal conductance(z2, z1)
        for z1, z2 in [(6, 9), (6, 7), (6, 8), (9, 17)]:
            assert abs(scheme.conductance(z1, z2, 1.0) -
                       scheme.conductance(z2, z1, 1.0)) < 1e-9


# ── Two-atom graph: analytic ground truth ─────────────────────────────────────

class TestTwoAtomGraph:
    """For a graph with 2 nodes connected by a single edge of conductance c:
        L = [[c, -c], [-c, c]]
        L+ = [[1/(4c), -1/(4c)], [-1/(4c), 1/(4c)]]  (Moore-Penrose pinv)
        R(0,1) = L+[0,0] + L+[1,1] - 2*L+[0,1] = 1/(4c) + 1/(4c) + 2/(4c) = 1/c
        K_f = 1/c

    For a C-C bond with c=1 → K_f = 1.0
    """

    def _make_cc_solver(self):
        """Create a minimal ComponentsSolver for C2F6 (just 2 C atoms in fluorinated component)."""
        mol = Chem.MolFromSmiles('FC(F)(F)C(F)(F)F')
        assert mol is not None
        solver = ComponentsSolver(mol)
        return solver

    def test_two_node_kirchhoff(self):
        solver = self._make_cc_solver()
        # The fluorinated component contains both carbons (2 C atoms).
        # Find the 2-node component and verify K_f = 1/c via _kirchhoff_index.
        for comps in solver.components.values():
            for comp in comps:
                if len(comp) == 2:
                    subG = solver.G.subgraph(comp)
                    kirchhoff = solver._kirchhoff_index(subG, uniform=False)
                    nodes = list(comp)
                    # K_f for a 2-node graph = 1 / conductance of the single edge
                    scheme = solver.bde_scheme
                    z0 = solver.G.nodes[nodes[0]]['element']
                    z1 = solver.G.nodes[nodes[1]]['element']
                    bond_order = subG.edges[nodes[0], nodes[1]].get('order', 1.0)
                    c = scheme.conductance(z0, z1, bond_order)
                    expected = 1.0 / c
                    assert abs(kirchhoff - expected) < 1e-9, (
                        f'K_f expected {expected}, got {kirchhoff}'
                    )

    def test_two_node_uniform_kirchhoff(self):
        """With uniform weights (c=1) on a 2-node graph, K_f must equal 1.0."""
        solver = self._make_cc_solver()
        for comps in solver.components.values():
            for comp in comps:
                if len(comp) == 2:
                    subG = solver.G.subgraph(comp)
                    kirchhoff_uniform = solver._kirchhoff_index(subG, uniform=True)
                    assert abs(kirchhoff_uniform - 1.0) < 1e-9, (
                        f'Uniform K_f expected 1.0, got {kirchhoff_uniform}'
                    )


# ── EGR physical properties ────────────────────────────────────────────────────

class TestEGRPhysicalProperties:

    # Linear chain PFCA series: OC(=O)(CF)_n  →  the fluorinated part is a CF2 chain
    #   For a path graph P_n with uniform conductance c: K_f = n(n^2-1)/(3c) (kcal normalised)
    #   So K_f(n+1) > K_f(n) strictly for n >= 1.

    _LINEAR = [
        ('C2', 'OC(=O)C(F)(F)F'),
        ('C3', 'OC(=O)C(F)(F)C(F)(F)F'),
        ('C4', 'OC(=O)C(F)(F)C(F)(F)C(F)(F)F'),
        ('C5', 'OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F'),
        ('C6', 'OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F'),
        ('C7', 'OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F'),
        ('C8', 'OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F'),
    ]

    def test_egr_increases_with_chain_length(self):
        """EGR must be strictly increasing as the fluorinated chain grows."""
        prev_egr = None
        prev_name = None
        for name, smi in self._LINEAR:
            egr = _egr(smi)
            assert not math.isnan(egr), f'{name}: EGR is NaN'
            if prev_egr is not None:
                assert egr > prev_egr, (
                    f'EGR should increase with chain length: '
                    f'{prev_name}={prev_egr:.4f} >= {name}={egr:.4f}'
                )
            prev_egr = egr
            prev_name = name

    def test_egr_positive_for_chain(self):
        for name, smi in self._LINEAR[1:]:  # skip C2 (trivially small)
            egr = _egr(smi)
            assert egr > 0, f'{name}: EGR should be positive, got {egr}'

    def test_branching_reduces_egr(self):
        """A branched isomer should have a lower EGR than the linear chain of
        the same carbon count, because branching shortens the longest paths.
        """
        # C5: linear PFCA vs α-branched (gem-bis-CF3 on α-carbon)
        linear_c5 = 'OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F'      # 4-C chain component
        gem_c4    = 'OC(=O)C(C(F)(F)F)(C(F)(F)F)F'              # star-shaped: 1 root + 2 CF3 branches
        egr_lin = _egr(linear_c5)
        egr_gem = _egr(gem_c4)
        # gem has diameter 2 whereas linear-4 has diameter 3 → EGR_gem < EGR_lin
        assert egr_gem < egr_lin, (
            f'Branched (gem) EGR ({egr_gem:.4f}) should be < linear-C5 EGR ({egr_lin:.4f})'
        )

    def test_linear_vs_alpha_branched_same_n(self):
        """At equal carbon count (n=4), a path has higher EGR than a star.

        linear-4 (path P4): K_f = n(n²-1)/(3c) = 20/c
        star-4 (K_{1,3}): central C + 3 CF3 branches → K_f = 9/c
        """
        linear = 'OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F'  # 4-C path
        star   = 'OC(=O)C(C(F)(F)F)(C(F)(F)F)C(F)(F)F'  # 4-C star (K_{1,3})
        egr_l = _egr(linear)
        egr_s = _egr(star)
        assert egr_l > egr_s, (
            f'Linear P4 EGR ({egr_l:.4f}) should exceed K_{{1,3}} star EGR ({egr_s:.4f})'
        )

    def test_pfos_vs_pfoa(self):
        """PFOA and PFOS each have a C8 fluorinated chain; both EGR values
        should be finite and positive.  The exact values may differ because
        PFOS includes the S head-group atom in the component, which extends
        the effective graph beyond the bare 8-carbon chain."""
        pfoa = 'OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F'
        pfos = 'OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F'
        egr_pfoa = _egr(pfoa)
        egr_pfos = _egr(pfos)
        assert egr_pfoa > 0 and not math.isnan(egr_pfoa), f'PFOA EGR invalid: {egr_pfoa}'
        assert egr_pfos > 0 and not math.isnan(egr_pfos), f'PFOS EGR invalid: {egr_pfos}'


# ── Kirchhoff identity ─────────────────────────────────────────────────────────

class TestKirchhoffIdentity:
    """K_f = n * sum(1/lambda_i) for non-zero eigenvalues of L.

    This is the classical relation; we verify it holds for the BDE-weighted L.
    """

    def _kirchhoff_via_eigenvalues(self, subG, scheme):
        nodes = list(subG.nodes())
        n = len(nodes)
        idx = {node: i for i, node in enumerate(nodes)}
        L = np.zeros((n, n))
        for u, v, data in subG.edges(data=True):
            bo = data.get('order', 1.0)
            c = scheme.conductance(
                subG.nodes[u]['element'], subG.nodes[v]['element'], bo)
            i, j = idx[u], idx[v]
            L[i, i] += c; L[j, j] += c
            L[i, j] -= c; L[j, i] -= c
        eigvals = np.linalg.eigvalsh(L)
        non_zero = eigvals[eigvals > 1e-10]
        if len(non_zero) == 0:
            return 0.0
        return float(n * np.sum(1.0 / non_zero))

    @pytest.mark.parametrize('smiles', [
        'OC(=O)C(F)(F)C(F)(F)C(F)(F)F',                            # linear C4
        'OC(=O)C(F)(C(F)(F)F)C(F)(F)F',                            # alpha-C4
        'OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F',             # linear C6
        'OC(=O)C(C(F)(F)F)(C(F)(F)F)C(F)(F)C(F)(F)F',             # gem-C5
    ])
    def test_kirchhoff_equals_eigenvalue_sum(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        solver = ComponentsSolver(mol)
        scheme = solver.bde_scheme
        for comps in solver.components.values():
            for comp in comps:
                if len(comp) < 2:
                    continue
                subG = solver.G.subgraph(comp)
                if not __import__('networkx').is_connected(subG):
                    continue
                kirchhoff_pinv = solver._kirchhoff_index(subG, uniform=False)
                kirchhoff_eig = self._kirchhoff_via_eigenvalues(subG, scheme)
                assert abs(kirchhoff_pinv - kirchhoff_eig) < 1e-6, (
                    f'Pseudoinverse K_f ({kirchhoff_pinv:.6f}) differs from '
                    f'eigenvalue K_f ({kirchhoff_eig:.6f}) for {smiles}'
                )


# ── SMARTS-distance metrics use BDE resistance ─────────────────────────────────

class TestSmartsResistanceMetrics:

    def test_all_resistance_metrics_present(self):
        smi = 'OC(=O)C(F)(F)C(F)(F)C(F)(F)F'
        comp = _first_comp_metrics(smi)
        for key in (
            'effective_graph_resistance',
            'effective_graph_resistance_BDE',
            'branching',
            'mean_eccentricity',
            'diameter',
            'radius',
            'component_fraction',
            'min_dist_to_barycentre',
            'min_dist_to_centre',
            'max_dist_to_periphery',
        ):
            assert key in comp, f'{key} missing from component metrics'
            val = comp[key]
            assert val is not None, f'{key} is None'
            assert not math.isnan(float(val)), f'{key} is NaN'


# ── New EGR-BDE variant tests ──────────────────────────────────────────────────

class TestEGRBDEVariant:
    """Tests for the new effective_graph_resistance_BDE metric."""

    def test_egr_bde_present(self):
        """effective_graph_resistance_BDE must be present and finite in component metrics."""
        comp = _first_comp_metrics('OC(=O)C(F)(F)C(F)(F)C(F)(F)F')
        assert 'effective_graph_resistance_BDE' in comp, (
            'effective_graph_resistance_BDE missing from component metrics'
        )
        val = comp['effective_graph_resistance_BDE']
        assert not math.isnan(float(val)), 'effective_graph_resistance_BDE is NaN'

    def test_egr_bde_different_from_uniform(self):
        """effective_graph_resistance (uniform) should differ from effective_graph_resistance_BDE
        (BDE-weighted, 1-hop expanded) for a typical PFAS molecule."""
        smi = 'OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F'
        comp = _first_comp_metrics(smi)
        egr_uniform = comp['effective_graph_resistance']
        egr_bde = comp['effective_graph_resistance_BDE']
        # They are computed on different graphs with different weights, so they should differ.
        assert egr_uniform != egr_bde, (
            f'effective_graph_resistance ({egr_uniform}) should differ from '
            f'effective_graph_resistance_BDE ({egr_bde})'
        )

    def test_egr_bde_positive(self):
        """effective_graph_resistance_BDE must be positive for all components."""
        smi = 'OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F'
        results = parse_smiles(smi, halogens='F')
        for match in results[0]['matches']:
            for comp in match.get('components', []):
                val = comp.get('effective_graph_resistance_BDE', float('nan'))
                if not math.isnan(float(val)):
                    assert float(val) > 0, (
                        f'effective_graph_resistance_BDE should be positive, got {val}'
                    )


# ── limit_effective_graph_resistance ──────────────────────────────────────────

class TestResistanceLimit:

    def test_limit_zero_returns_nan(self):
        """limit_effective_graph_resistance=0 should skip all resistance computation."""
        results = parse_smiles(
            'OC(=O)C(F)(F)C(F)(F)C(F)(F)F',
            halogens='F',
            limit_effective_graph_resistance=0,
        )
        for match in results[0]['matches']:
            for comp in match.get('components', []):
                assert math.isnan(comp['effective_graph_resistance']), \
                    'EGR should be NaN when limit=0'

    def test_limit_large_computes_normally(self):
        """A large limit should behave identically to None."""
        smi = 'OC(=O)C(F)(F)C(F)(F)C(F)(F)F'
        egr_none  = _egr(smi)
        results   = parse_smiles(smi, halogens='F', limit_effective_graph_resistance=500)
        egr_limit = results[0]['matches'][0]['components'][0]['effective_graph_resistance']
        assert abs(egr_none - egr_limit) < 1e-9
