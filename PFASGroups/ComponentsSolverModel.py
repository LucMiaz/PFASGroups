import os
import json
import math
import warnings
import numpy as np
import networkx as nx
from rdkit import Chem
from .core import add_componentSmarts, mol_to_nx

# ── BDE helpers (self-contained, no dependency on molecular_quantum_graph) ────

def _bde_keys_to_int(x):
    return {int(k) if isinstance(k, str) and k.isdigit() else k: v
            for k, v in x.items()}


def _load_bde_dict_local():
    """Load diatomic BDE dict (kcal/mol) from PFASGroups' own data folder."""
    _path = os.path.join(os.path.dirname(__file__), 'data', 'diatomic_bonds_dict.json')
    try:
        with open(_path) as fh:
            raw = json.load(fh, object_hook=_bde_keys_to_int)
        return raw
    except Exception as exc:
        warnings.warn(f'PFASGroups: could not load BDE dict ({exc}); '
                      'falling back to uniform resistance.')
        return None


# ── Hard-coded bond-order scaling model ──────────────────────────────────────
# Best model selected by Psi4 B3LYP/6-31G* calibration on 134 diatomic
# molecules (see molecular_quantum_graph/bde_computation/bond_order_calibration/
# for the full analysis).  Model: poly2  f(n) = 1 + a*(n-1) + b*(n-1)^2
_BOND_ORDER_MODEL_NAME  = 'poly2'
_BOND_ORDER_MODEL_PARAMS = {'a': 1.2650122708517233, 'b': -0.3142013833397031}


def _bond_order_factor(n: float, model_name: str, params: dict) -> float:
    """Evaluate the bond-order scaling factor f(n), f(1)=1."""
    x = n - 1.0
    if model_name == 'linear':
        val = 1.0 + params.get('alpha', 0.3) * x
    elif model_name == 'power':
        val = n ** params.get('beta', 0.6)
    elif model_name == 'log':
        val = 1.0 + params.get('a', 1.0) * math.log(max(n, 1e-10))
    elif model_name == 'poly2':
        val = 1.0 + params.get('a', 0.3) * x + params.get('b', 0.0) * x ** 2
    elif model_name == 'poly3':
        val = (1.0
               + params.get('a', 0.3) * x
               + params.get('b', 0.0) * x ** 2
               + params.get('c', 0.0) * x ** 3)
    else:
        val = 1.0 + 0.3 * x
    return max(val, 1e-6)


class _BDEScheme:
    """Lightweight BDE weighting scheme bundled with PFASGroups."""

    def __init__(self):
        self.bde_dict = _load_bde_dict_local()
        self._model_name   = _BOND_ORDER_MODEL_NAME
        self._model_params = _BOND_ORDER_MODEL_PARAMS
        # C-C single-bond BDE as normalisation reference (kcal/mol)
        self.ref_bde = (
            self.bde_dict.get(6, {}).get(6, 83.1)
            if self.bde_dict else 83.1
        )

    def conductance(self, z1: int, z2: int, bond_order: float = 1.0) -> float:
        """Return BDE conductance = BDE(z1,z2,order) / ref_bde.

        Higher BDE → stronger bond → higher conductance → shorter resistance path.
        Returns 1.0 (uniform) if BDE data is unavailable.
        """
        if self.bde_dict is None:
            return 1.0
        try:
            base = self.bde_dict[z1][z2]
        except KeyError:
            try:
                base = self.bde_dict[z2][z1]
            except KeyError:
                b1 = self.bde_dict.get(z1, {}).get(6, 80.0)
                b2 = self.bde_dict.get(z2, {}).get(6, 80.0)
                base = (b1 + b2) / 2.0
        bde = base * _bond_order_factor(bond_order, self._model_name, self._model_params)
        return bde / self.ref_bde


# Module-level singleton — loaded once per process
_BDE_SCHEME: '_BDEScheme | None' = None


def _get_bde_scheme() -> _BDEScheme:
    global _BDE_SCHEME
    if _BDE_SCHEME is None:
        _BDE_SCHEME = _BDEScheme()
    return _BDE_SCHEME

class ComponentsSolver:
    """Class to hold components information with comprehensive graph metrics."""
    @add_componentSmarts()
    def __init__(self, mol, **kwargs):
        self.componentSmartss = kwargs.get('componentSmartss')
        self.mol = mol
        self.mol_size = mol.GetNumAtoms()  # Total atoms in molecule for fraction calculation
        self.total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')  # Total carbon atoms

        # Per-halogen density metrics (F, Cl, Br, I)
        # Naming convention:
        #   total_{sym}s        : total count of that halogen (total_fluorines, total_chlorines, …)
        #   per{X}ination_density: halogen count / heavy-atom count (perfluorination_density, …)
        #   c{sym}2_count       : carbons bearing ≥2 of that halogen (cf2_count, ccl2_count, …)
        #   c{sym}2_density     : c{sym}2_count / total_carbons
        _HALOGEN_META = [
            ('F',  'fluorines',  'perfluorination_density',   'cf2_count',  'cf2_density'),
            ('Cl', 'chlorines',  'perchlorination_density',  'ccl2_count', 'ccl2_density'),
            ('Br', 'bromines',   'perbromination_density',   'cbr2_count', 'cbr2_density'),
            ('I',  'iodines',    'periodination_density',    'ci2_count',  'ci2_density'),
        ]
        for halogen_sym, total_attr, per_attr, cx2_count_attr, cx2_density_attr in _HALOGEN_META:
            total_hal = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == halogen_sym)
            setattr(self, 'total_' + total_attr, total_hal)
            setattr(self, per_attr, total_hal / self.mol_size if self.mol_size > 0 else 0.0)
            cx2 = sum(
                1 for atom in mol.GetAtoms()
                if atom.GetSymbol() == 'C' and
                sum(1 for nb in atom.GetNeighbors() if nb.GetSymbol() == halogen_sym) >= 2
            )
            setattr(self, cx2_count_attr, cx2)
            setattr(self, cx2_density_attr, cx2 / self.total_carbons if self.total_carbons > 0 else 0.0)
        self.G = mol_to_nx(mol)
        self.bde_scheme = _get_bde_scheme()
        self.limit_effective_graph_resistance = kwargs.get('limit_effective_graph_resistance', None)
        self.skip_component_metrics = not kwargs.get('compute_component_metrics', True)
        self.total_branching = self._compute_total_branching()
        self.components = self.get_fluorinated_subgraph()
        self.extended_components = {k:{0:v} for k,v in self.components.items()}
        # Mapping from (pathType, max_dist, extended_component_index) -> original_component_index
        self.component_to_original_index = {}
        self.levels = {0}
        # Cache for component metrics
        self._component_metrics_cache = {}
        # Precompute full component sizes (including all attached atoms: H, F, Cl, Br, I)
        self.component_full_sizes = {}
        for path_type, components_list in self.components.items():
            if path_type not in self.component_full_sizes:
                self.component_full_sizes[path_type] = {}
            for i, comp in enumerate(components_list):
                full_comp = self.get_full_component_atoms(comp)
                self.component_full_sizes[path_type][i] = len(full_comp)
        # Compute and store metrics for all components on creation (if enabled)
        if not self.skip_component_metrics:
            self._precompute_component_metrics()
        # Initialize mapping for level 0 (original components)
        self._init_component_mapping()

    def __enter__(self):
        return self
    def __exit__(self, exc_type, exc_value, traceback):
        self.components = None
        self.extended_components = None
        self.mol = None
        self.G = None
        self._component_metrics_cache = None

    def __len__(self):
        return len(self.components)

    def get(self, pathType, max_dist=0, default = []):
        if max_dist not in self.levels:
            self.extend_components(max_dist)
        return self.extended_components.get(pathType, {}).get(max_dist, default)

    def max_size(self):
        return max([len(x) for x in self.components]) if len(self.components)>0 else 0

    def sizes(self):
        return [len(x) for x in self.components]

    def _connected_components(self, subset):
        """Find connected components in a molecule."""
        G = self.G.subgraph(subset)
        components = list(nx.connected_components(G))
        return components

    @staticmethod
    def _extract_component_smarts(entry):
        if entry is None:
            return None
        if isinstance(entry, (list, tuple)):
            return entry[0] if len(entry) > 0 else None
        if isinstance(entry, dict):
            return entry.get('smarts', entry.get('component', entry.get('chain')))
        return entry

    def get_fluorinated_subgraph(self, **kwargs):
        """Get the fluorinated indices by connected components of a molecule based on path SMARTS."""
        subsets = {}
        for pathName, d in self.componentSmartss.items():
            path_smarts = self._extract_component_smarts(d)
            if path_smarts is None:
                subsets[pathName] = []
                continue
            if isinstance(path_smarts, str):
                path_smarts = Chem.MolFromSmarts(path_smarts)
                if path_smarts is None:
                    subsets[pathName] = []
                    continue
                path_smarts.UpdatePropertyCache()
                Chem.GetSymmSSSR(path_smarts)
                path_smarts.GetRingInfo().NumRings()
            matches = self.mol.GetSubstructMatches(path_smarts)
            subset = [y for x in matches for y in x]
            if len(subset)==0:
                subsets[pathName]= []
                continue
            components = self._connected_components(subset)
            subsets[pathName] = components
        return subsets

    def _compute_total_branching(self):
        non_hf_atoms = [
            atom.GetIdx()
            for atom in self.mol.GetAtoms()
            if atom.GetSymbol() not in ['H', 'F','Cl','Br','I']
        ]
        if not non_hf_atoms:
            return 0.0
        return calculate_branching(self.G, non_hf_atoms)

    def get_full_component_atoms(self, component):
        """Get all atoms in a component including those attached (H, F, halogens).

        Parameters
        ----------
        component : set
            Set of atom indices representing the carbon backbone

        Returns
        -------
        set
            Set of all atom indices including the component and all directly attached atoms
        """
        full_component = set(component)
        # Add all neighbors of component atoms that are not already in the component
        # This includes H, F, and other halogens attached to the carbon backbone
        for atom_idx in component:
            for neighbor_idx in self.G.neighbors(atom_idx):
                if self.mol.GetAtomWithIdx(neighbor_idx).GetSymbol() in ['H', 'F', 'Cl', 'Br', 'I']:
                    full_component.add(neighbor_idx)
        return full_component

    def get_total_components_fraction(self, matched_components_list):
        """Calculate the fraction of carbon atoms in the molecule covered by the union of all components.

        Parameters
        ----------
        matched_components_list : list of dict
            List of matched component dictionaries, each with 'component' and 'smarts_matches' keys

        Returns
        -------
        float
            Fraction of carbon atoms covered by the union of all components (0.0 to 1.0)
        """
        if len(matched_components_list) == 0 or self.total_carbons == 0:
            return 0.0

        # Union all carbon atoms from components (augmented components already include
        # SMARTS-match and linker atoms, so just iterate over 'component' atom sets).
        union_carbon_atoms = set()
        for comp_dict in matched_components_list:
            component = set(comp_dict.get('component', []))
            smarts_matches = comp_dict.get('smarts_matches')

            # Add carbon atoms from component
            for atom_idx in component:
                if self.mol.GetAtomWithIdx(atom_idx).GetSymbol() == 'C':
                    union_carbon_atoms.add(atom_idx)

            # smarts_matches are intentionally excluded: the augmented component already
            # contains those atoms via get_augmented_component. Adding them separately
            # would double-count for OECD groups and overcount for telomers.

        # Total coverage = union of C atoms from all augmented components.
        # smarts_extra_atoms is intentionally not added here: the augmented components
        # already incorporate SMARTS-match and linker atoms, so no additive correction
        # is needed and it would push telomers above 1.0.
        total_fraction = len(union_carbon_atoms) / self.total_carbons
        return total_fraction

    def _precompute_component_metrics(self):
        """Precompute metrics for all initial components."""
        for path_type, components_list in self.components.items():
            for comp in components_list:
                # This will populate the cache
                self.compute_component_metrics(comp)

    def _init_component_mapping(self):
        """Initialize mapping for original components (max_dist=0)."""
        for pathType in self.components.keys():
            for i in range(len(self.components[pathType])):
                self.component_to_original_index[(pathType, 0, i)] = i

    def extend_components(self, max_dist):
        """Extend a component in a graph by a maximum distance.
        This is used to match functional groups that are not directly connected to the component. Different components that overlap are not merged.
        """
        if max_dist>0:
            for pathType, components in self.components.items():
                extended_components = []
                for i, component in enumerate(components):
                    extended = component.copy()
                    for node in component:
                        lengths = nx.single_source_shortest_path_length(self.G, node, cutoff=max_dist)
                        extended.update([n for n,d in lengths.items() if d<=max_dist])
                    extended_components.append(extended)
                    # Map this extended component back to its original component
                    self.component_to_original_index[(pathType, max_dist, i)] = i
                self.extended_components.setdefault(pathType, {})[max_dist] = extended_components
            self.levels.add(max_dist)
    def shortest_path_to_component(self, atom, component):
        """Get shortest paths from SMARTS matches to original component within extended component.

        Parameters
        ----------
        atom_index : int
        component : set

        Returns
        -------
        dict
            Mapping from SMARTS match atom index to shortest path list to original component
        """
        path = nx.shortest_path(self.G, atom, list(component)[0])
        path = set(path).difference(component)
        if len(path)==0:
            return None, []
        return [x for x in path if x!=atom]
    def get_augmented_component(self, pathType, max_dist, component_index, smarts_matches, linker_smarts=None):
        """Get original component augmented with connecting atoms to SMARTS matches.

        Parameters
        ----------
        pathType : str
            Type of component path (e.g., 'Perfluoroalkyl')
        max_dist : int
            Distance used for extension
        component_index : int
            Index of the component in the extended components list
        smarts_matches : set
            Set of atom indices that matched the SMARTS pattern
        linker_smarts : Chem.Mol, optional
            Compiled SMARTS pattern for validating linker atoms.
            If provided, only paths where intermediate atoms match this pattern are accepted.

        Returns
        -------
        set
            Original component augmented with shortest path atoms connecting SMARTS matches
        """
        if max_dist == 0:
            # No augmentation needed, return original component
            return self.components[pathType][component_index]

        # Get original and extended components
        orig_index = self.component_to_original_index.get((pathType, max_dist, component_index), component_index)
        orig_comp = self.components[pathType][orig_index]
        ext_comp = self.extended_components[pathType][max_dist][component_index]

        # Start with original component
        augmented = set(orig_comp)
        linker_matches = set([y for x in self.mol.GetSubstructMatches(linker_smarts) for y in x]) if linker_smarts is not None else []
        # Add shortest paths from SMARTS matches to original component
        for smarts_atom in smarts_matches:
            if smarts_atom in ext_comp and smarts_atom not in orig_comp:
                try:
                    linker_atoms = self.shortest_path_to_component(smarts_atom, orig_comp)


                    # Validate linker atoms only if there are intermediate atoms
                    # Direct connections (no linker) are accepted without validation
                    if linker_smarts is not None and len(linker_atoms) > 0:
                        # Validate the intermediate linker atoms
                        if not all(atom in linker_matches for atom in linker_atoms):
                            continue  # Skip this SMARTS atom if linker validation fails
                except nx.NetworkXNoPath:
                    continue  # Skip this SMARTS atom if no path exists
                # Add the complete path: SMARTS atom + linker atoms + component border atom
                augmented.update([smarts_atom] + linker_atoms)
            elif smarts_atom in orig_comp:
                # SMARTS atom already in original component
                augmented.add(smarts_atom)

        # Verify that augmented component still contains the SMARTS matches
        # Count how many SMARTS atoms ended up in the augmented component
        smarts_in_augmented = sum(1 for atom in smarts_matches if atom in augmented)

        # If no SMARTS atoms made it into the augmented component, reject it
        # This happens when all SMARTS atoms failed linker validation
        if smarts_in_augmented == 0 and len(smarts_matches) > 0:
            return []

        return augmented

    def _kirchhoff_index_pinv(self, subG, uniform: bool = False) -> float:
        """Kirchhoff index via Laplacian pseudoinverse. Fast for n < 30."""
        nodes = list(subG.nodes())
        n = len(nodes)
        idx = {node: i for i, node in enumerate(nodes)}

        L = np.zeros((n, n), dtype=float)
        for u, v, data in subG.edges(data=True):
            if uniform:
                c = 1.0
            else:
                bond_order = data.get('order', 1.0)
                z_u = subG.nodes[u].get('element', 6)
                z_v = subG.nodes[v].get('element', 6)
                c = self.bde_scheme.conductance(z_u, z_v, bond_order)
            i, j = idx[u], idx[v]
            L[i, i] += c
            L[j, j] += c
            L[i, j] -= c
            L[j, i] -= c

        Lp = np.linalg.pinv(L)
        diag = np.diag(Lp)
        kirchhoff = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                kirchhoff += max(diag[i] + diag[j] - 2.0 * Lp[i, j], 0.0)
        return kirchhoff

    def _kirchhoff_index(self, subG, uniform: bool = False) -> float:
        """Kirchhoff (effective graph resistance) index via Laplacian eigenspectrum.

        Equivalent to nx.effective_graph_resistance (Theorem 2.2, Ellens et al.
        2011) but without the internal G.copy(). Differences from nx:
          - Uses np.linalg.eigvalsh (symmetric-aware, returns sorted reals)
            instead of np.linalg.eigvals; results are identical for symmetric L.
          - Builds the weighted Laplacian directly from edge data rather than
            relying on graph weight attributes, enabling the BDE-weighted variant.

        Parameters
        ----------
        subG : networkx.Graph
            Subgraph to operate on (original component or 1-hop expanded).
        uniform : bool
            True  → all edge conductances = 1 (topological / unweighted).
            False → BDE-calibrated conductances (bond-strength weighted).

        Returns
        -------
        float
            Sum of all pairwise effective resistance distances.
        """
        n = subG.number_of_nodes()
        nodes = list(subG.nodes())
        idx = {node: i for i, node in enumerate(nodes)}

        L = np.zeros((n, n), dtype=float)
        for u, v, data in subG.edges(data=True):
            if uniform:
                c = 1.0
            else:
                bond_order = data.get('order', 1.0)
                z_u = subG.nodes[u].get('element', 6)
                z_v = subG.nodes[v].get('element', 6)
                c = self.bde_scheme.conductance(z_u, z_v, bond_order)
            i, j = idx[u], idx[v]
            L[i, i] += c
            L[j, j] += c
            L[i, j] -= c
            L[j, i] -= c

        # Eigenvalues only; skip the zero eigenvalue (index 0)
        mu = np.sort(np.linalg.eigvalsh(L))
        return float(np.sum(1.0 / mu[1:]) * n)

    def compute_component_metrics(self, component):
        """Compute comprehensive graph metrics for a component.

        Parameters
        ----------
        component : set or frozenset
            Set of atom indices in the component

        Returns
        -------
        dict
            Dictionary with graph metrics including:
            - diameter: maximum eccentricity
            - radius: minimum eccentricity
            - eccentricity_values: dict mapping node to its eccentricity
            - center: nodes with minimum eccentricity
            - periphery: nodes with maximum eccentricity
            - barycenter: nodes minimizing sum of distances
            - effective_graph_resistance: BDE-weighted Kirchhoff index
            - _rdist: internal BDE-weighted pairwise resistance distance dict
        """
        # If metrics computation is disabled, return minimal metrics (size only)
        if self.skip_component_metrics:
            return {
                'size': len(component),
                'diameter': float('nan'),
                'radius': float('nan'),
                'eccentricity_values': {},
                'center': [],
                'periphery': [],
                'barycenter': [],
                'effective_graph_resistance': float('nan'),
                'effective_graph_resistance_BDE': float('nan'),
            }

        cache_key = frozenset(component)
        if cache_key in self._component_metrics_cache:
            return self._component_metrics_cache[cache_key]

        if len(component) <= 1:
            metrics = {
                'diameter': 0,
                'radius': 0,
                'eccentricity_values': {list(component)[0]: 0} if len(component) == 1 else {},
                'center': list(component),
                'periphery': list(component),
                'barycenter': list(component),
                'effective_graph_resistance': 0.0,
                'effective_graph_resistance_BDE': 0.0,
            }
            self._component_metrics_cache[cache_key] = metrics
            return metrics

        # Create subgraph for this component
        subG = self.G.subgraph(component)

        # Check if connected
        if not nx.is_connected(subG):
            metrics = {
                'diameter': float('inf'),
                'radius': 0,
                'eccentricity_values': {},
                'center': [],
                'periphery': [],
                'barycenter': [],
                'effective_graph_resistance': float('inf'),
                'effective_graph_resistance_BDE': float('inf'),
            }
            self._component_metrics_cache[cache_key] = metrics
            return metrics

        try:
            # Compute eccentricity for each node
            eccentricity_values = nx.eccentricity(subG)

            # Diameter and radius
            diameter = nx.diameter(subG)
            radius = nx.radius(subG)

            # Center and periphery
            center = nx.center(subG)
            periphery = nx.periphery(subG)

            # Barycenter: nodes that minimize total distance to all other nodes
            # total_distances = {}
            # for node in subG.nodes():
            #     lengths = nx.single_source_shortest_path_length(subG, node)
            #     total_distances[node] = sum(lengths.values())

            # min_total_dist = min(total_distances.values())
            # barycenter = [node for node, dist in total_distances.items() if dist == min_total_dist]
            barycenter = nx.barycenter(subG)

            # Effective graph resistance (Kirchhoff index) — two variants:
            #   uniform:  topological (edge weights = 1), original C-skeleton component
            #   BDE:      bond-strength-weighted, component expanded 1 hop to include F/H/Cl/Br…
            try:
                should_compute_resistance = (
                    self.limit_effective_graph_resistance is None or
                    (self.limit_effective_graph_resistance > 0
                     and len(component) < self.limit_effective_graph_resistance)
                )

                if should_compute_resistance:
                    effective_graph_resistance = self._kirchhoff_index(subG, uniform=True)
                    # 1-hop expansion: add all neighbours of every component node
                    expanded = set(component)
                    for node in list(component):
                        expanded.update(self.G.neighbors(node))
                    subG_exp = self.G.subgraph(expanded)
                    effective_graph_resistance_BDE = self._kirchhoff_index(subG_exp, uniform=False)
                else:
                    effective_graph_resistance = float('nan')
                    effective_graph_resistance_BDE = float('nan')
            except Exception:
                effective_graph_resistance = float('nan')
                effective_graph_resistance_BDE = float('nan')

            metrics = {
                'diameter': diameter,
                'radius': radius,
                'eccentricity_values': eccentricity_values,
                'center': center,
                'periphery': periphery,
                'barycenter': barycenter,
                'effective_graph_resistance': effective_graph_resistance,
                'effective_graph_resistance_BDE': effective_graph_resistance_BDE,
            }

        except Exception as e:
            # Fallback for any computation errors
            metrics = {
                'diameter': float('nan'),
                'radius': float('nan'),
                'eccentricity_values': {},
                'center': [],
                'periphery': [],
                'barycenter': [],
                'effective_graph_resistance': float('nan'),
                'effective_graph_resistance_BDE': float('nan'),
            }

        self._component_metrics_cache[cache_key] = metrics
        return metrics

    def compute_smarts_component_metrics(self, component, smarts_matches):
        """Compute metrics relating SMARTS matches to component structural features.

        Parameters
        ----------
        component : set
            Set of atom indices in the component
        smarts_matches : set
            Set of atom indices matching the SMARTS pattern

        Returns
        -------
        dict
            Dictionary with SMARTS-specific metrics, or None if no smarts_matches
        """
        if smarts_matches is None or len(smarts_matches) == 0:
            return None

        comp_metrics = self.compute_component_metrics(component)
        smarts_in_comp = smarts_matches.intersection(component)

        if len(smarts_in_comp) == 0 or len(component) <= 1:
            return {
                'min_dist_to_barycenter': 0,
                'min_dist_to_center': 0,
                'max_dist_to_periphery': 0,
            }

        subG = self.G.subgraph(component)

        min_dist_to_barycenter = float('inf')
        min_dist_to_center = float('inf')
        max_dist_to_periphery = 0

        try:
            for smarts_node in smarts_in_comp:
                if smarts_node not in subG:
                    continue

                for bc_node in comp_metrics['barycenter']:
                    try:
                        dist = nx.shortest_path_length(subG, smarts_node, bc_node)
                        min_dist_to_barycenter = min(min_dist_to_barycenter, dist)
                    except Exception:
                        pass

                for center_node in comp_metrics['center']:
                    try:
                        dist = nx.shortest_path_length(subG, smarts_node, center_node)
                        min_dist_to_center = min(min_dist_to_center, dist)
                    except Exception:
                        pass

                for periph_node in comp_metrics['periphery']:
                    try:
                        dist = nx.shortest_path_length(subG, smarts_node, periph_node)
                        max_dist_to_periphery = max(max_dist_to_periphery, dist)
                    except Exception:
                        pass

            if min_dist_to_barycenter == float('inf'):
                min_dist_to_barycenter = 0
            if min_dist_to_center == float('inf'):
                min_dist_to_center = 0

        except Exception:
            min_dist_to_barycenter = 0
            min_dist_to_center = 0
            max_dist_to_periphery = 0

        return {
            'min_dist_to_barycenter': min_dist_to_barycenter,
            'min_dist_to_center': min_dist_to_center,
            'max_dist_to_periphery': max_dist_to_periphery,
        }

    def get_matched_component_dict(self, component, smarts_matches=None, smarts_type='unknown', pfas_group=None, comp_id = None):
        """Get a dictionary with all metrics for a matched component.

        Parameters
        ----------
        component : set
            Set of atom indices in the component
        smarts_matches : set or None
            Set of atom indices matching the SMARTS pattern (None if no SMARTS)
        smarts_type : str
            Type identifier for the SMARTS pattern
        pfas_group : HalogenGroup or None
            HalogenGroup object with precomputed SMARTS atom counts

        Returns
        -------
        dict
            Complete dictionary with all component metrics
        """
        # Basic branching metric
        smarts_set = smarts_matches if smarts_matches is not None else set()
        basic_metrics = calculate_component_metrics(self.G, component, smarts_set)

        # Comprehensive graph metrics (cached)
        comp_metrics = self.compute_component_metrics(component)

        # SMARTS-specific metrics (computed on the fly, None if no smarts)
        smarts_metrics = self.compute_smarts_component_metrics(component, smarts_matches)

        # Get precomputed SMARTS extra atoms count if pfas_group is available
        smarts_extra_atoms = 0
        if pfas_group is not None and smarts_matches is not None and len(smarts_matches) > 0:
            if pfas_group.smarts_extra_atoms is not None:
                # Sum the extra atoms from all SMARTS patterns
                # For groups with multiple matches, we count each match
                smarts_extra_atoms =  pfas_group.component_specific_extra_atoms[comp_id] if comp_id is not None else sum(pfas_group.smarts_extra_atoms) * len(smarts_matches)

        # Calculate mean and median eccentricity from eccentricity_values
        eccentricity_values = comp_metrics.get('eccentricity_values', {})
        if len(eccentricity_values) > 0:
            ecc_list = list(eccentricity_values.values())
            mean_eccentricity = sum(ecc_list) / len(ecc_list)
            sorted_ecc = sorted(ecc_list)
            n = len(sorted_ecc)
            if n % 2 == 0:
                median_eccentricity = (sorted_ecc[n//2 - 1] + sorted_ecc[n//2]) / 2.0
            else:
                median_eccentricity = sorted_ecc[n//2]
        else:
            mean_eccentricity = 0.0
            median_eccentricity = 0.0

        # Calculate component fraction based on carbon atoms only
        # Count carbon atoms in component
        component_carbons = sum(1 for atom_idx in component if self.mol.GetAtomWithIdx(atom_idx).GetSymbol() == 'C')

        # Component-level halogen density metrics (F, Cl, Br, I)
        component_carbons_indices = [idx for idx in component
                                     if self.mol.GetAtomWithIdx(idx).GetSymbol() == 'C']
        _n_comp_c = len(component_carbons_indices)
        _COMP_HALOGEN_META = [
            ('F',  'component_f_count',  'component_perfluorination_density',  'component_cf2_count',  'component_cf2_density'),
            ('Cl', 'component_cl_count', 'component_perchlorination_density', 'component_ccl2_count', 'component_ccl2_density'),
            ('Br', 'component_br_count', 'component_perbromination_density',  'component_cbr2_count', 'component_cbr2_density'),
            ('I',  'component_i_count',  'component_periodination_density',   'component_ci2_count',  'component_ci2_density'),
        ]
        _comp_hal_vals = {}
        for halogen_sym, cnt_key, per_key, cx2_cnt_key, cx2_den_key in _COMP_HALOGEN_META:
            hal_count = sum(
                1 for atom_idx in component_carbons_indices
                for nb in self.mol.GetAtomWithIdx(atom_idx).GetNeighbors()
                if nb.GetSymbol() == halogen_sym
            )
            cx2_count = sum(
                1 for atom_idx in component_carbons_indices
                if sum(1 for nb in self.mol.GetAtomWithIdx(atom_idx).GetNeighbors()
                       if nb.GetSymbol() == halogen_sym) >= 2
            )
            _comp_hal_vals[cnt_key]    = hal_count
            _comp_hal_vals[per_key]    = hal_count / _n_comp_c if _n_comp_c > 0 else 0.0
            _comp_hal_vals[cx2_cnt_key] = cx2_count
            _comp_hal_vals[cx2_den_key] = cx2_count / _n_comp_c if _n_comp_c > 0 else 0.0
        # Convenience aliases for backward compatibility
        component_f_count = _comp_hal_vals['component_f_count']
        component_cf2_count = _comp_hal_vals['component_cf2_count']

        # component_fraction: ratio of carbon atoms in this component to all carbon atoms in the
        # molecule, counting only C atoms actually present in the (augmented) component atom set.
        # The augmented component already includes linker atoms and SMARTS-match atoms added by
        # get_augmented_component, so no further additive correction is needed.
        # Excluded: O, F, and other heteroatoms that may appear in telomer linker paths.
        component_fraction = component_carbons / self.total_carbons if self.total_carbons > 0 else 0.0

        # --- n_spacer (telomer CH₂ linker chain length) -------------------------
        # Only meaningful for groups that have a linker_smarts (telomers).
        # Formula: the augmented component includes the orig pfluorinated component
        # PLUS the CH₂ linker atoms PLUS the SMARTS-match atom (the functional group
        # C adjacent to the linker).  So n_spacer = |augmented - orig_comp| + 1.
        # For n=1: the SMARTS-match atom IS already in orig_comp, so ΔΔ=0 → spacer=1.
        n_spacer = 0
        if pfas_group is not None and getattr(pfas_group, 'linker_smarts', None) is not None and comp_id is not None:
            max_dist_for_spacer = getattr(pfas_group, 'max_dist_from_comp', 0)
            orig_idx = self.component_to_original_index.get(
                (smarts_type, max_dist_for_spacer, comp_id), comp_id
            )
            if smarts_type in self.components and orig_idx < len(self.components[smarts_type]):
                orig_comp_set = set(self.components[smarts_type][orig_idx])
                n_spacer = len(set(component) - orig_comp_set) + 1

        # --- ring_size (smallest ring overlapping with component) ---------------
        ring_size = 0
        component_ring_atoms = [
            idx for idx in component if self.mol.GetAtomWithIdx(idx).IsInRing()
        ]
        if component_ring_atoms:
            ring_info = self.mol.GetRingInfo()
            comp_set = set(component)
            for ring in sorted(ring_info.AtomRings(), key=len):
                if comp_set & set(ring):
                    ring_size = len(ring)
                    break

        result = {
            'component': sorted(list(component)),
            'size': len(component),
            'component_fraction': component_fraction,  # C atoms in augmented component / total C atoms in molecule (always ≤ 1.0)
            'smarts_matches': sorted(list(smarts_matches)) if smarts_matches is not None else None,  # Store for union calculation
            'smarts_extra_atoms': smarts_extra_atoms,  # Extra carbons from functional group
            'SMARTS': smarts_type,
            # Telomer spacer and cyclic ring metrics
            'n_spacer': n_spacer,
            'ring_size': ring_size,
            # Basic metrics
            'branching': basic_metrics['branching'],
            'branching_ratio_to_molecule': basic_metrics['branching'] / self.total_branching if self.total_branching > 0 else 0.0,
            'total_branching': self.total_branching,
            'smarts_centrality': basic_metrics['smarts_centrality'],
            # Graph structure metrics
            'diameter': comp_metrics.get('diameter', float('nan')),
            'radius': comp_metrics.get('radius', float('nan')),
            'effective_graph_resistance': comp_metrics.get('effective_graph_resistance', float('nan')),
            'effective_graph_resistance_BDE': comp_metrics.get('effective_graph_resistance_BDE', float('nan')),
            'eccentricity_values': comp_metrics.get('eccentricity_values', {}),
            'mean_eccentricity': mean_eccentricity,
            'median_eccentricity': median_eccentricity,
            'center': comp_metrics.get('center', []),
            'periphery': comp_metrics.get('periphery', []),
            'barycenter': comp_metrics.get('barycenter', []),
            # Distance metrics (with defaults)
            'min_dist_to_barycenter': 0,
            'min_dist_to_center': 0,
            'max_dist_to_periphery': 0,
            # Halogen density metrics — component-level (F, Cl, Br, I)
            **_comp_hal_vals,
            # Global molecule-level halogen density (for context)
            'molecule_perfluorination_density':   self.perfluorination_density,
            'molecule_cf2_density':               self.cf2_density,
            'molecule_perchlorination_density':   self.perchlorination_density,
            'molecule_ccl2_density':              self.ccl2_density,
            'molecule_perbromination_density':    self.perbromination_density,
            'molecule_cbr2_density':              self.cbr2_density,
            'molecule_periodination_density':     self.periodination_density,
            'molecule_ci2_density':               self.ci2_density,
        }

        # Override distance metrics with SMARTS-specific values if available
        if smarts_metrics is not None:
            result.update({
                'min_dist_to_barycenter': smarts_metrics.get('min_dist_to_barycenter', 0),
                'min_dist_to_center': smarts_metrics.get('min_dist_to_center', 0),
                'max_dist_to_periphery': smarts_metrics.get('max_dist_to_periphery', 0),
            })

        return result


def calculate_branching(G, subset):
    # Calculate branching: measure of branching vs linearity
    # For linear chains: branching → 1.0
    # For highly branched: branching → 0.0
    try:
        subG = G.subgraph(subset)
        # Count branch points (degree > 2)
        branch_points = sum(1 for node in subG.nodes() if subG.degree(node) > 2)
        # Normalize by component size
        branching = 1.0 - (branch_points / max(1, len(subset) - 2))  # -2 to account for endpoints
        branching = max(0.0, min(1.0, branching))  # Clamp to [0, 1]
    except:
        branching = 0.0
    return branching

def calculate_component_metrics(G, component, smarts_matches):
    """Calculate branching and centrality metrics for a component.

    Parameters
    ----------
    G : networkx.Graph
        Molecular graph
    component : set
        Set of atom indices in the component
    smarts_matches : set
        Set of atom indices matching the SMARTS pattern

    Returns
    -------
    dict
        Dictionary with 'branching' (float) and 'smarts_centrality' (float)
    """
    if len(component) <= 1:
        return {'branching': 0.0, 'smarts_centrality': 1.0}

    # Create subgraph for this component
    subG = G.subgraph(component)

    # Calculate branching: fraction of nodes with degree > 2
    branching = calculate_branching(subG, component)

    # Calculate SMARTS centrality: how central the matched atoms are
    smarts_in_component = smarts_matches.intersection(component)
    if len(smarts_in_component) == 0:
        smarts_centrality = 0.0
    else:
        try:
            # Calculate average shortest path distance from SMARTS matches to all other nodes
            total_distance = 0
            count = 0
            for smarts_node in smarts_in_component:
                if smarts_node in subG:
                    lengths = nx.single_source_shortest_path_length(subG, smarts_node)
                    for node, dist in lengths.items():
                        if node != smarts_node:
                            total_distance += dist
                            count += 1

            if count > 0:
                avg_distance = total_distance / count
                # Calculate maximum possible average distance (for peripheral node)
                # For a linear chain of n nodes, max avg distance is ~n/3
                max_possible_distance = len(component) / 3.0
                # Centrality: 1.0 = central, 0.0 = peripheral
                smarts_centrality = 1.0 - min(1.0, avg_distance / max(1.0, max_possible_distance))
            else:
                smarts_centrality = 1.0
        except:
            smarts_centrality = 0.5

    return {'branching': branching, 'smarts_centrality': smarts_centrality}
