import networkx as nx
from rdkit import Chem
from .core import add_componentSmarts, mol_to_nx

class ComponentsSolver:
    """Class to hold components information with comprehensive graph metrics."""
    @add_componentSmarts()
    def __init__(self, mol, **kwargs):
        self.componentSmartss = kwargs.get('componentSmartss')
        self.mol = mol
        self.mol_size = mol.GetNumAtoms()  # Total atoms in molecule for fraction calculation
        self.total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')  # Total carbon atoms
        self.G = mol_to_nx(mol)
        self.limit_effective_graph_resistance = kwargs.get('limit_effective_graph_resistance',None)
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
            if atom.GetSymbol() not in ['H', 'F']
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

        # Union all carbon atoms from components and SMARTS matches
        union_carbon_atoms = set()
        for comp_dict in matched_components_list:
            component = set(comp_dict.get('component', []))
            smarts_matches = comp_dict.get('smarts_matches')

            # Add carbon atoms from component
            for atom_idx in component:
                if self.mol.GetAtomWithIdx(atom_idx).GetSymbol() == 'C':
                    union_carbon_atoms.add(atom_idx)

            # Add carbon atoms from SMARTS matches
            if smarts_matches is not None:
                for atom_idx in smarts_matches:
                    if self.mol.GetAtomWithIdx(atom_idx).GetSymbol() == 'C':
                        union_carbon_atoms.add(atom_idx)

        # Add extra carbons from functional groups
        # This accounts for carbons in functional groups that are not matched by SMARTS
        # For example, in -COOH, if SMARTS matches the adjacent carbon, we need to add
        # the carbonyl carbon from smarts_extra_atoms
        extra_carbons_count = 0
        for comp_dict in matched_components_list:
            # Get smarts_extra_atoms from pfas_group if available
            if 'smarts_extra_atoms' in comp_dict:
                extra_carbons_count += comp_dict['smarts_extra_atoms']

        # Total = union of matched carbons + extra functional group carbons
        # Cap at total_carbons to avoid exceeding 1.0
        total_carbon_count = min(len(union_carbon_atoms) + extra_carbons_count, self.total_carbons)
        total_fraction = total_carbon_count / self.total_carbons
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
            - effective_graph_resistance: sum of resistance distances
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
                'effective_graph_resistance': float('nan')
            }

        # Use frozenset for caching
        comp_key = frozenset(component)
        if comp_key in self._component_metrics_cache:
            return self._component_metrics_cache[comp_key]

        if len(component) <= 1:
            metrics = {
                'diameter': 0,
                'radius': 0,
                'eccentricity_values': {list(component)[0]: 0} if len(component) == 1 else {},
                'center': list(component),
                'periphery': list(component),
                'barycenter': list(component),
                'effective_graph_resistance': 0.0
            }
            self._component_metrics_cache[comp_key] = metrics
            return metrics

        # Create subgraph for this component
        subG = self.G.subgraph(component)

        # Check if connected
        if not nx.is_connected(subG):
            # For disconnected components, compute metrics separately
            metrics = {
                'diameter': float('inf'),
                'radius': 0,
                'eccentricity_values': {},
                'center': [],
                'periphery': [],
                'barycenter': [],
                'effective_graph_resistance': float('inf')
            }
            self._component_metrics_cache[comp_key] = metrics
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
            total_distances = {}
            for node in subG.nodes():
                lengths = nx.single_source_shortest_path_length(subG, node)
                total_distances[node] = sum(lengths.values())

            min_total_dist = min(total_distances.values())
            barycenter = [node for node, dist in total_distances.items() if dist == min_total_dist]

            # Effective graph resistance (requires matrix operations)
            try:
                # Compute resistance distance based on limit setting
                # None = compute for all, int > 0 = compute if component size < limit, 0 = skip all
                should_compute_resistance = (
                    self.limit_effective_graph_resistance is None or
                    (self.limit_effective_graph_resistance > 0 and len(component) < self.limit_effective_graph_resistance)
                )

                if should_compute_resistance:
                    resistance_sum = 0.0
                    nodes = list(subG.nodes())
                    for i, u in enumerate(nodes):
                        for v in nodes[i+1:]:
                            try:
                                # Resistance distance approximation using shortest path
                                # For more accurate computation, would need Laplacian pseudoinverse
                                sp_length = nx.shortest_path_length(subG, u, v)
                                resistance_sum += sp_length
                            except:
                                resistance_sum += float('inf')
                    effective_graph_resistance = resistance_sum
                else:
                    effective_graph_resistance = float('nan')
            except:
                effective_graph_resistance = float('nan')

            metrics = {
                'diameter': diameter,
                'radius': radius,
                'eccentricity_values': eccentricity_values,
                'center': center,
                'periphery': periphery,
                'barycenter': barycenter,
                'effective_graph_resistance': effective_graph_resistance
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
                'effective_graph_resistance': float('nan')
            }

        self._component_metrics_cache[comp_key] = metrics
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
                'min_resistance_dist_to_barycenter': 0.0,
                'min_dist_to_center': 0,
                'min_resistance_dist_to_center': 0.0,
                'max_dist_to_periphery': 0,
                'max_resistance_dist_to_periphery': 0.0
            }

        subG = self.G.subgraph(component)

        # Initialize metrics
        min_dist_to_barycenter = float('inf')
        min_resistance_to_barycenter = float('inf')
        min_dist_to_center = float('inf')
        min_resistance_to_center = float('inf')
        max_dist_to_periphery = 0
        max_resistance_to_periphery = 0.0

        try:
            # Compute distances from SMARTS matches to structural features
            for smarts_node in smarts_in_comp:
                if smarts_node not in subG:
                    continue

                # Distances to barycenter nodes
                for bc_node in comp_metrics['barycenter']:
                    try:
                        dist = nx.shortest_path_length(subG, smarts_node, bc_node)
                        min_dist_to_barycenter = min(min_dist_to_barycenter, dist)
                        # Resistance distance approximation
                        min_resistance_to_barycenter = min(min_resistance_to_barycenter, float(dist))
                    except:
                        pass

                # Distances to center nodes
                for center_node in comp_metrics['center']:
                    try:
                        dist = nx.shortest_path_length(subG, smarts_node, center_node)
                        min_dist_to_center = min(min_dist_to_center, dist)
                        min_resistance_to_center = min(min_resistance_to_center, float(dist))
                    except:
                        pass

                # Distances to periphery nodes
                for periph_node in comp_metrics['periphery']:
                    try:
                        dist = nx.shortest_path_length(subG, smarts_node, periph_node)
                        max_dist_to_periphery = max(max_dist_to_periphery, dist)
                        max_resistance_to_periphery = max(max_resistance_to_periphery, float(dist))
                    except:
                        pass

            # Handle cases where no valid distances were found
            if min_dist_to_barycenter == float('inf'):
                min_dist_to_barycenter = 0
                min_resistance_to_barycenter = 0.0
            if min_dist_to_center == float('inf'):
                min_dist_to_center = 0
                min_resistance_to_center = 0.0

        except Exception as e:
            # Fallback values
            min_dist_to_barycenter = 0
            min_resistance_to_barycenter = 0.0
            min_dist_to_center = 0
            min_resistance_to_center = 0.0
            max_dist_to_periphery = 0
            max_resistance_to_periphery = 0.0

        return {
            'min_dist_to_barycenter': min_dist_to_barycenter,
            'min_resistance_dist_to_barycenter': min_resistance_to_barycenter,
            'min_dist_to_center': min_dist_to_center,
            'min_resistance_dist_to_center': min_resistance_to_center,
            'max_dist_to_periphery': max_dist_to_periphery,
            'max_resistance_dist_to_periphery': max_resistance_to_periphery
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

        # Count carbon atoms in SMARTS matches that are NOT already in the component
        smarts_carbons_not_in_component = 0
        if smarts_matches is not None:
            for atom_idx in smarts_matches:
                if self.mol.GetAtomWithIdx(atom_idx).GetSymbol() == 'C' and atom_idx not in component:
                    smarts_carbons_not_in_component += 1

        # smarts_extra_atoms represents additional carbons in the functional group beyond the matched carbon
        # For each SMARTS pattern, extra_atoms indicates carbons beyond the primary matched atom
        # These are summed across all SMARTS matches to get the total extra carbons
        # So we always add smarts_extra_atoms regardless of whether matched carbon is in component

        total_carbons_in_component = component_carbons + smarts_carbons_not_in_component + smarts_extra_atoms

        component_fraction = total_carbons_in_component / self.total_carbons if self.total_carbons > 0 else 0.0

        result = {
            'component': sorted(list(component)),
            'size': len(component),
            'component_fraction': component_fraction,  # Fraction of molecule covered by component
            'smarts_matches': sorted(list(smarts_matches)) if smarts_matches is not None else None,  # Store for union calculation
            'smarts_extra_atoms': smarts_extra_atoms,  # Extra carbons from functional group
            'SMARTS': smarts_type,
            # Basic metrics
            'branching': basic_metrics['branching'],
            'branching_ratio_to_molecule': basic_metrics['branching'] / self.total_branching if self.total_branching > 0 else 0.0,
            'total_branching': self.total_branching,
            'smarts_centrality': basic_metrics['smarts_centrality'],
            # Graph structure metrics
            'diameter': comp_metrics.get('diameter', float('nan')),
            'radius': comp_metrics.get('radius', float('nan')),
            'effective_graph_resistance': comp_metrics.get('effective_graph_resistance', float('nan')),
            'eccentricity_values': comp_metrics.get('eccentricity_values', {}),
            'mean_eccentricity': mean_eccentricity,
            'median_eccentricity': median_eccentricity,
            'center': comp_metrics.get('center', []),
            'periphery': comp_metrics.get('periphery', []),
            'barycenter': comp_metrics.get('barycenter', []),
            # Distance metrics (with defaults)
            'min_dist_to_barycenter': 0,
            'min_resistance_dist_to_barycenter': 0.0,
            'min_dist_to_center': 0,
            'min_resistance_dist_to_center': 0.0,
            'max_dist_to_periphery': 0,
            'max_resistance_dist_to_periphery': 0.0
        }

        # Override distance metrics with SMARTS-specific values if available
        if smarts_metrics is not None:
            result.update({
                'min_dist_to_barycenter': smarts_metrics.get('min_dist_to_barycenter', 0),
                'min_resistance_dist_to_barycenter': smarts_metrics.get('min_resistance_dist_to_barycenter', 0.0),
                'min_dist_to_center': smarts_metrics.get('min_dist_to_center', 0),
                'min_resistance_dist_to_center': smarts_metrics.get('min_resistance_dist_to_center', 0.0),
                'max_dist_to_periphery': smarts_metrics.get('max_dist_to_periphery', 0),
                'max_resistance_dist_to_periphery': smarts_metrics.get('max_resistance_dist_to_periphery', 0.0)
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
