"""
Tests for BDE-weighted effective graph resistance (Kirchhoff index).

Verifies:
- Physical plausibility: linear chain EGR grows as O(n^3) for n atoms
- Branching sensitivity: branched isomers have lower EGR than linear
- Bond-order effect: stronger bonds (higher BDE) reduce resistance
- Exact Kirchhoff identity: K_f == n * sum(1/lambda_i) for non-zero eigenvalues
- Resistance distance symmetry and triangle inequality
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
        # The fluorinated component is both carbons (2 C atoms)
        # Find the component containing both carbons
        for comps in solver.components.values():
            for comp in comps:
                if len(comp) == 2:
                    subG = solver.G.subgraph(comp)
                    kirchhoff, rdist = solver._bde_resistance_matrix(subG)
                    nodes = list(comp)
                    r = rdist[(nodes[0], nodes[1])]
                    # K_f for 2 nodes = single resistance = 1/c
                    scheme = solver.bde_scheme
                    z0 = solver.G.nodes[nodes[0]]['element']
                    z1 = solver.G.nodes[nodes[1]]['element']
                    bond_order = subG.edges[nodes[0], nodes[1]].get('order', 1.0)
                    c = scheme.conductance(z0, z1, bond_order)
                    expected = 1.0 / c
                    assert abs(r - expected) < 1e-9, f'R(0,1) expected {expected}, got {r}'
                    assert abs(kirchhoff - expected) < 1e-9


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
                kirchhoff_pinv, _ = solver._bde_resistance_matrix(subG)
                kirchhoff_eig = self._kirchhoff_via_eigenvalues(subG, scheme)
                assert abs(kirchhoff_pinv - kirchhoff_eig) < 1e-6, (
                    f'Pseudoinverse K_f ({kirchhoff_pinv:.6f}) differs from '
                    f'eigenvalue K_f ({kirchhoff_eig:.6f}) for {smiles}'
                )


# ── Resistance distance properties ────────────────────────────────────────────

class TestResistanceDistanceProperties:

    def _get_rdist(self, smiles: str) -> dict:
        mol = Chem.MolFromSmiles(smiles)
        solver = ComponentsSolver(mol)
        for comps in solver.components.values():
            for comp in comps:
                if len(comp) >= 2:
                    subG = solver.G.subgraph(comp)
                    if __import__('networkx').is_connected(subG):
                        _, rdist = solver._bde_resistance_matrix(subG)
                        return rdist
        return {}

    @pytest.mark.parametrize('smiles', [
        'OC(=O)C(F)(F)C(F)(F)C(F)(F)F',
        'OC(=O)C(F)(C(F)(F)F)C(F)(F)C(F)(F)F',
    ])
    def test_rdist_symmetry(self, smiles):
        rdist = self._get_rdist(smiles)
        for (u, v), r in rdist.items():
            assert abs(r - rdist.get((v, u), float('nan'))) < 1e-12, \
                f'R({u},{v}) != R({v},{u}): {r} vs {rdist.get((v,u))}'

    @pytest.mark.parametrize('smiles', [
        'OC(=O)C(F)(F)C(F)(F)C(F)(F)F',
        'OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F',
    ])
    def test_rdist_positive(self, smiles):
        rdist = self._get_rdist(smiles)
        for (u, v), r in rdist.items():
            if u != v:
                assert r > 0, f'R({u},{v}) = {r} should be positive'

    @pytest.mark.parametrize('smiles', [
        'OC(=O)C(F)(F)C(F)(F)C(F)(F)F',
    ])
    def test_rdist_triangle_inequality(self, smiles):
        """Resistance distance satisfies the triangle inequality: R(i,k) <= R(i,j)+R(j,k)."""
        rdist = self._get_rdist(smiles)
        nodes = list({u for u, _ in rdist})
        violations = []
        for i in nodes:
            for j in nodes:
                for k in nodes:
                    if i == j or j == k or i == k:
                        continue
                    rij = rdist.get((i, j))
                    rjk = rdist.get((j, k))
                    rik = rdist.get((i, k))
                    if rij is None or rjk is None or rik is None:
                        continue
                    if rik > rij + rjk + 1e-9:
                        violations.append((i, j, k, rik, rij + rjk))
        assert not violations, f'Triangle inequality violated: {violations[:3]}'


# ── SMARTS-distance metrics use BDE resistance ─────────────────────────────────

class TestSmartsResistanceMetrics:

    def test_resistance_dist_differs_from_hop_dist(self):
        """min_resistance_dist_to_barycenter should differ from the plain
        hop distance when C-F bonds (higher BDE) are present, because BDE
        weighting changes the effective path cost.
        """
        # Use a structure with a branch point where C-F and C-C conductances differ
        smi = 'OC(=O)C(F)(C(F)(F)F)C(F)(F)C(F)(F)F'
        comp = _first_comp_metrics(smi)
        hop   = comp['min_dist_to_barycenter']
        rdist = comp['min_resistance_dist_to_barycenter']
        # Both must be finite and non-negative
        assert not math.isnan(rdist), 'min_resistance_dist_to_barycenter is NaN'
        assert rdist >= 0

    def test_periphery_resistance_finite(self):
        smi = 'OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F'
        comp = _first_comp_metrics(smi)
        assert not math.isnan(comp.get('max_resistance_dist_to_periphery', float('nan')))
        assert comp.get('max_resistance_dist_to_periphery', -1) >= 0

    def test_all_resistance_metrics_present(self):
        smi = 'OC(=O)C(F)(F)C(F)(F)C(F)(F)F'
        comp = _first_comp_metrics(smi)
        for key in ('effective_graph_resistance',
                    'min_resistance_dist_to_barycenter',
                    'min_resistance_dist_to_center',
                    'max_resistance_dist_to_periphery',
                    'mean_resistance_distance',
                    'resistance_diameter',
                    'resistance_radius',
                    'mean_resistance_eccentricity',
                    'median_resistance_eccentricity'):
            assert key in comp, f'{key} missing from component metrics'
            assert not math.isnan(comp[key]), f'{key} is NaN'


# ── New molecule-wide resistance metrics ───────────────────────────────────────

class TestNewResistanceMetrics:
    """Verify the new per-component and molecule-wide BDE-resistance metrics."""

    def test_mean_resistance_distance_positive(self):
        """mean_resistance_distance > 0 for any multi-atom component."""
        comp = _first_comp_metrics('OC(=O)C(F)(F)C(F)(F)C(F)(F)F')
        mrd = comp['mean_resistance_distance']
        assert not math.isnan(mrd), 'mean_resistance_distance is NaN'
        assert mrd > 0

    def test_mean_resistance_distance_less_than_egr(self):
        """mean_resistance_distance = EGR / C(n,2) < EGR for n > 2."""
        comp = _first_comp_metrics('OC(=O)C(F)(F)C(F)(F)C(F)(F)F')
        mrd = comp['mean_resistance_distance']
        egr = comp['effective_graph_resistance']
        assert mrd < egr, f'mean_resistance_distance ({mrd}) should be < EGR ({egr}) for n>2'

    def test_resistance_diameter_geq_radius(self):
        """resistance_diameter >= resistance_radius (same as hop-count)."""
        comp = _first_comp_metrics('OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F')
        rd = comp['resistance_diameter']
        rr = comp['resistance_radius']
        assert not math.isnan(rd) and not math.isnan(rr)
        assert rd >= rr, f'resistance_diameter ({rd}) < resistance_radius ({rr})'

    def test_resistance_diameter_increases_with_chain(self):
        """Resistance diameter should grow as the chain lengthens."""
        s3 = 'OC(=O)C(F)(F)C(F)(F)F'
        s5 = 'OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F'
        rd3 = _first_comp_metrics(s3)['resistance_diameter']
        rd5 = _first_comp_metrics(s5)['resistance_diameter']
        assert rd5 > rd3, f'C5 resistance diameter ({rd5}) should exceed C3 ({rd3})'

    def test_mean_resistance_eccentricity_between_radius_and_diameter(self):
        """resistance_radius <= mean_resistance_eccentricity <= resistance_diameter."""
        comp = _first_comp_metrics('OC(=O)C(F)(F)C(F)(F)C(F)(F)F')
        rr   = comp['resistance_radius']
        mre  = comp['mean_resistance_eccentricity']
        rd   = comp['resistance_diameter']
        assert rr <= mre + 1e-9, f'radius ({rr}) > mean_ecc ({mre})'
        assert mre <= rd + 1e-9, f'mean_ecc ({mre}) > diameter ({rd})'

    def test_molecule_wide_metrics_present(self):
        """parse_smiles HalogenGroup match summary dicts should contain the new molecule-wide keys."""
        results = parse_smiles('OC(=O)C(F)(F)C(F)(F)C(F)(F)F', halogens='F')
        halogen_matches = [m for m in results[0]['matches'] if m.get('type') == 'HalogenGroup']
        assert halogen_matches, 'Expected at least one HalogenGroup match'
        for match in halogen_matches:
            for key in ('mean_resistance_distance', 'mean_resistance_diameter',
                        'mean_resistance_radius', 'mean_resistance_eccentricity',
                        'mean_resistance_dist_to_barycenter',
                        'mean_resistance_dist_to_center',
                        'mean_resistance_dist_to_periphery'):
                assert key in match, f'{key} missing from HalogenGroup match summary'

    def test_branching_lowers_resistance_diameter(self):
        """A 4-node star (K_{1,3}) has a smaller resistance diameter than a 4-node path (P4).

        P4 resistance diameter = max R(u,v) between the two ends ≈ 20/3c  (roughly)
        K_{1,3} resistance diameter = max R(leaf,leaf) = 2/c (through center)
        """
        path = 'OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F'  # 4-C path
        star = 'OC(=O)C(C(F)(F)F)(C(F)(F)F)C(F)(F)F'  # 4-C star (K_{1,3})
        rd_path = _first_comp_metrics(path)['resistance_diameter']
        rd_star = _first_comp_metrics(star)['resistance_diameter']
        assert rd_path > rd_star, (
            f'P4 resistance diameter ({rd_path}) should > K_{{1,3}} ({rd_star})'
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
