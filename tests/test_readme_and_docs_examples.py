"""Smoke tests verifying that README and documentation code examples still run.

Each test corresponds to one or more code blocks from:
  - README.md       (Quick Start, Embeddings, Multi-Halogen, Custom Configuration sections)
  - docs/getting_started.rst
  - docs/fingerprints.rst

The goal is to catch early when API changes break published examples.

Run with:  pytest tests/test_readme_and_docs_examples.py
"""

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# Shared SMILES constants — reused across tests to avoid repetition
# ---------------------------------------------------------------------------

# Four representative PFAS from the README embedding section
SMILES_FP = [
    "O=C(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",  # PFOA
    "O=C(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",          # PFHpA
    "OCC(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                              # 4:2 FTOH
    "OCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",               # 6:2 FTOH
]

# Mixed-halogen SMILES from the Multi-Halogen section
SMILES_MULTI_HAL = [
    "C(C(F)(F)F)F",            # fluorinated
    "ClC(Cl)(Cl)C(Cl)(Cl)Cl",  # chlorinated
    "BrC(Br)(Br)CBr",          # brominated
]

# Three molecules from the README benchmark block
SMILES_BENCHMARK = [
    "C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",
    "C1=CC(=CC=C1[N+](=O)[O-])OC2=C(C=C(C=C2)C(F)(F)F)[N+](=O)[O-]",
    "C(C(C(C(F)(F)S(=O)(=O)O)(F)F)(F)F)(C(C(F)(F)F)(F)F)(C(F)(F)F)C(F)(F)F",
]

# Three molecules used in the fingerprints.rst examples
SMILES_FP_DOCS = [
    "FC(F)(F)C(F)(F)C(=O)O",     # PFOA
    "FC(F)(F)C(F)(F)S(=O)(=O)O", # PFHxS
    "C(C(F)(F)F)F",               # HFC-134
]

# ---------------------------------------------------------------------------
# Module-scoped fixtures — parse once, reuse across many tests
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def results_fp():
    """Parsed PFASEmbeddingSet for SMILES_FP, default halogens."""
    from PFASGroups import parse_smiles
    return parse_smiles(SMILES_FP)


@pytest.fixture(scope="module")
def results_multi_hal():
    """Parsed PFASEmbeddingSet for SMILES_MULTI_HAL with all four halogens."""
    from PFASGroups import parse_smiles
    return parse_smiles(SMILES_MULTI_HAL, halogens=["F", "Cl", "Br", "I"])


# ===========================================================================
# README — Verify Installation / Quick Start
# ===========================================================================


def test_readme_quickstart_single_smiles():
    """README 'Verify installation': parse_smiles returns a PFASEmbeddingSet; str() works."""
    from PFASGroups import parse_smiles
    from PFASGroups.PFASEmbeddings import PFASEmbeddingSet

    results = parse_smiles("FC(F)(F)C(F)(F)C(=O)O")
    assert isinstance(results, PFASEmbeddingSet)
    assert len(results) == 1
    str(results)    # must not raise
    str(results[0]) # must not raise


def test_readme_quickstart_list():
    """README Quick Start: parse_smiles with a list returns correct length."""
    from PFASGroups import parse_smiles

    results = parse_smiles(SMILES_FP)
    assert len(results) == len(SMILES_FP)
    str(results)
    str(results[0])


# ===========================================================================
# README — Metric control switches
# ===========================================================================


def test_readme_metric_control_no_metrics():
    """README benchmark block: compute_component_metrics=False — fastest path."""
    from PFASGroups import parse_smiles

    r = parse_smiles(SMILES_BENCHMARK, compute_component_metrics=False)
    assert len(r) == len(SMILES_BENCHMARK)


def test_readme_metric_control_limit_egr_zero():
    """README benchmark block: limit_effective_graph_resistance=0 skips EGR."""
    from PFASGroups import parse_smiles

    r = parse_smiles(SMILES_BENCHMARK, limit_effective_graph_resistance=0)
    assert len(r) == len(SMILES_BENCHMARK)


def test_readme_metric_control_limit_egr_threshold():
    """README benchmark block: limit_effective_graph_resistance=20 uses size threshold."""
    from PFASGroups import parse_smiles

    r = parse_smiles(SMILES_BENCHMARK, limit_effective_graph_resistance=20)
    assert len(r) == len(SMILES_BENCHMARK)


# ===========================================================================
# README — to_array shapes and column_names consistency
# ===========================================================================


def test_readme_to_array_default_shape(results_fp):
    """README: arr_bin = results.to_array() — shape matches column_names()."""
    arr = results_fp.to_array()
    cols = results_fp.column_names()
    assert arr.ndim == 2
    assert arr.shape[0] == len(SMILES_FP)
    assert arr.shape[1] == len(cols)


def test_readme_to_array_oecd_selection(results_fp):
    """README: to_array(group_selection='oecd', component_metrics=['binary'])."""
    arr = results_fp.to_array(group_selection="oecd", component_metrics=["binary"])
    cols = results_fp.column_names(group_selection="oecd", component_metrics=["binary"])
    assert arr.shape[0] == len(SMILES_FP)
    assert arr.shape[1] == len(cols)


def test_readme_to_array_two_metrics_double_width(results_fp):
    """README: two component_metrics → exactly 2× default column count."""
    n_default = results_fp.to_array().shape[1]
    arr = results_fp.to_array(component_metrics=["binary", "effective_graph_resistance"])
    cols = results_fp.column_names(component_metrics=["binary", "effective_graph_resistance"])
    assert arr.shape[1] == len(cols)
    assert arr.shape[1] == 2 * n_default


def test_readme_to_array_n_spacer(results_fp):
    """README: component_metrics=['n_spacer'] — single block, same width as default."""
    n_default = results_fp.to_array().shape[1]
    arr = results_fp.to_array(component_metrics=["n_spacer"])
    assert arr.shape[1] == n_default


def test_readme_to_array_ring_size(results_fp):
    """README: component_metrics=['ring_size'] — single block, same width as default."""
    n_default = results_fp.to_array().shape[1]
    arr = results_fp.to_array(component_metrics=["ring_size"])
    assert arr.shape[1] == n_default


def test_readme_to_array_combined(results_fp):
    """README combined: 4 component metrics + 6 molecule metrics."""
    n_default = results_fp.to_array().shape[1]
    mol_metrics = [
        "n_components", "max_size", "mean_branching", "max_branching",
        "mean_component_fraction", "max_component_fraction",
    ]
    arr = results_fp.to_array(
        component_metrics=["binary", "effective_graph_resistance", "n_spacer", "ring_size"],
        molecule_metrics=mol_metrics,
    )
    cols = results_fp.column_names(
        component_metrics=["binary", "effective_graph_resistance", "n_spacer", "ring_size"],
        molecule_metrics=mol_metrics,
    )
    assert arr.shape == (len(SMILES_FP), n_default * 4 + len(mol_metrics))
    assert arr.shape[1] == len(cols)


def test_readme_column_names_are_strings(results_fp):
    """README: cols[:4] should be non-empty strings."""
    cols = results_fp.column_names(component_metrics=["binary", "effective_graph_resistance"])
    assert all(isinstance(c, str) and c for c in cols)


# ===========================================================================
# README — FINGERPRINT_PRESETS
# ===========================================================================


def test_readme_presets_access(results_fp):
    """README: FINGERPRINT_PRESETS['best'] is accessible and has a 'description'."""
    from PFASGroups import FINGERPRINT_PRESETS

    assert "best" in FINGERPRINT_PRESETS
    p = FINGERPRINT_PRESETS["best"]
    assert "description" in p
    assert "component_metrics" in p


def test_readme_preset_best_double_width(results_fp):
    """README: to_array(preset='best') → 2 blocks (binary + EGR)."""
    from PFASGroups import FINGERPRINT_PRESETS

    n_default = results_fp.to_array().shape[1]
    arr = results_fp.to_array(preset="best")
    cols = results_fp.column_names(preset="best")
    assert arr.shape[0] == len(SMILES_FP)
    assert arr.shape[1] == len(cols)
    # 'best' = ['binary', 'effective_graph_resistance'] → 2 blocks
    assert len(FINGERPRINT_PRESETS["best"]["component_metrics"]) == 2
    assert arr.shape[1] == 2 * n_default


def test_readme_preset_best3_expand(results_fp):
    """fingerprints.rst: expand preset manually, pass component_metrics directly."""
    from PFASGroups import FINGERPRINT_PRESETS

    p = FINGERPRINT_PRESETS["best_3"]
    arr = results_fp.to_array(component_metrics=p["component_metrics"])
    cols = results_fp.column_names(component_metrics=p["component_metrics"])
    assert arr.shape[0] == len(SMILES_FP)
    assert arr.shape[1] == len(cols)


def test_readme_all_named_presets_exist():
    """fingerprints.rst preset table: all documented presets are present."""
    from PFASGroups import FINGERPRINT_PRESETS

    for name in ["best", "best_2", "best_3", "best_4", "best_5", "binary", "count", "max_component"]:
        assert name in FINGERPRINT_PRESETS, f"Preset {name!r} missing"
        assert "component_metrics" in FINGERPRINT_PRESETS[name]


# ===========================================================================
# README — Filtering components by halogen, form, saturation
# ===========================================================================


def test_readme_filter_halogens_f_only():
    """README filtering section: parse_smiles(smiles, halogens='F')."""
    from PFASGroups import parse_smiles

    results = parse_smiles(SMILES_FP, halogens="F")
    assert len(results) == len(SMILES_FP)


def test_readme_filter_per_alkyl():
    """README filtering section: halogens='F', saturation='per', form='alkyl'."""
    from PFASGroups import parse_smiles

    results = parse_smiles(SMILES_FP, halogens="F", saturation="per", form="alkyl")
    assert len(results) == len(SMILES_FP)


def test_readme_filter_poly_cyclic():
    """README filtering section: halogens='F', saturation='poly', form='cyclic'."""
    from HalogenGroups import parse_smiles

    results = parse_smiles(SMILES_MULTI_HAL, halogens="F", saturation="poly", form="cyclic")
    assert len(results) == len(SMILES_MULTI_HAL)


def test_readme_filter_cyclic_form():
    """README Quick Start section: parse_smiles(smiles, form='cyclic')."""
    from PFASGroups import parse_smiles

    results = parse_smiles(SMILES_FP, form="cyclic")
    assert len(results) == len(SMILES_FP)


def test_readme_filter_multi_halogen():
    """README filtering section: halogens=['F', 'Cl']."""
    from HalogenGroups import parse_smiles

    results = parse_smiles(SMILES_MULTI_HAL, halogens=["F", "Cl"])
    assert len(results) == len(SMILES_MULTI_HAL)


# ===========================================================================
# README — Multi-halogen to_array with per-halogen blocks
# ===========================================================================


def test_readme_multi_halogen_to_array_shape(results_multi_hal):
    """README: r_hal.to_array(halogens=[...]) shape equals len(column_names(...))."""
    halogens = ["F", "Cl", "Br", "I"]
    arr = results_multi_hal.to_array(
        component_metrics=["effective_graph_resistance"],
        molecule_metrics=["n_components", "max_size", "mean_branching"],
        halogens=halogens,
    )
    cols = results_multi_hal.column_names(
        component_metrics=["effective_graph_resistance"],
        molecule_metrics=["n_components", "max_size", "mean_branching"],
        halogens=halogens,
    )
    assert arr.shape[0] == len(SMILES_MULTI_HAL)
    assert arr.shape[1] == len(cols)


def test_readme_multi_halogen_blocks_sum(results_multi_hal):
    """4-block width == sum of 4 individual single-halogen block widths."""
    halogens = ["F", "Cl", "Br", "I"]
    arr_4x = results_multi_hal.to_array(halogens=halogens)
    individual_widths = sum(
        results_multi_hal.to_array(halogens=h).shape[1] for h in halogens
    )
    assert arr_4x.shape[1] == individual_widths


def test_readme_multi_halogen_excludeHalogens_respected(results_multi_hal):
    """Each per-halogen block must exclude the matching halide-ion group."""
    cols_F = results_multi_hal.column_names(halogens="F")
    assert not any(c.startswith("fluoride") for c in cols_F), \
        "'fluoride' group should be excluded from the F block"
    assert any(c.startswith("bromide") for c in cols_F), \
        "'bromide' group should be included in the F block"

    cols_Cl = results_multi_hal.column_names(halogens="Cl")
    assert not any(c.startswith("chloride") for c in cols_Cl), \
        "'chloride' group should be excluded from the Cl block"

    cols_Br = results_multi_hal.column_names(halogens="Br")
    assert not any(c.startswith("bromide") for c in cols_Br), \
        "'bromide' group should be excluded from the Br block"

    cols_I = results_multi_hal.column_names(halogens="I")
    assert not any(c.startswith("iodide") for c in cols_I), \
        "'iodide' group should be excluded from the I block"


# ===========================================================================
# README — Multi-halogen Option A (HalogenGroups) and Option B (PFASgroups)
# ===========================================================================


def test_readme_option_a_halogengroups():
    """README Option A: `from HalogenGroups import parse_smiles` defaults to all halogens."""
    from HalogenGroups import parse_smiles

    results = parse_smiles(SMILES_MULTI_HAL)
    arr = results.to_array(
        group_selection="oecd",
        component_metrics=["binary"],
        halogens=["F", "Cl", "Br", "I"],
    )
    cols = results.column_names(
        group_selection="oecd",
        halogens=["F", "Cl", "Br", "I"],
    )
    assert arr.shape[0] == len(SMILES_MULTI_HAL)
    assert arr.shape[1] == len(cols)


def test_readme_option_b_pfasgroups():
    """README Option B: `from PFASGroups import parse_smiles` with explicit halogens list."""
    from PFASGroups import parse_smiles  # noqa: PLC0415  # only reached when installed

    results = parse_smiles(SMILES_MULTI_HAL, halogens=["F", "Cl", "Br", "I"])
    arr = results.to_array(
        group_selection="oecd",
        component_metrics=["binary"],
        halogens=["F", "Cl", "Br", "I"],
    )
    cols = results.column_names(
        group_selection="oecd",
        halogens=["F", "Cl", "Br", "I"],
    )
    assert arr.shape[1] == len(cols)


# ===========================================================================
# README — KL divergence comparison
# ===========================================================================


def test_readme_compare_kld():
    """README: results.compare_kld(other_results, method='minmax') returns a non-negative float."""
    from PFASGroups import parse_smiles

    other_smiles = SMILES_FP[:2] + [
        "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(O)=O",
        "OCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
    ]
    results = parse_smiles(SMILES_FP)
    other_results = parse_smiles(other_smiles)
    similarity = results.compare_kld(other_results, method="minmax")
    assert isinstance(similarity, float)
    assert similarity >= 0.0


# ===========================================================================
# README — prioritise_molecules
# ===========================================================================


def test_readme_prioritise_reference_based():
    """README: prioritise_molecules with reference list and return_scores=True."""
    from PFASGroups import prioritise_molecules

    reference = ["FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O"]
    prioritised, scores = prioritise_molecules(
        molecules=SMILES_FP,
        reference=reference,
        return_scores=True,
    )
    assert len(prioritised) == len(SMILES_FP)
    assert len(scores) == len(SMILES_FP)
    assert all("smiles" in p for p in prioritised)


def test_readme_prioritise_intrinsic():
    """README: prioritise_molecules with intrinsic fluorination properties."""
    from PFASGroups import prioritise_molecules

    prioritised = prioritise_molecules(
        molecules=SMILES_FP,
        a=1.0,
        b=5.0,
        percentile=90,
        return_scores=False,
    )
    assert len(prioritised) == len(SMILES_FP)
    assert all("smiles" in p for p in prioritised)


# ===========================================================================
# README — Custom configuration: preprocess_componentsSmarts
# ===========================================================================


def test_readme_custom_paths_preprocess():
    """README: get_componentSMARTSs + preprocess_componentsSmarts + parse."""
    from PFASGroups import get_componentSMARTSs, get_HalogenGroups, parse_smiles
    from PFASGroups.core import preprocess_componentsSmarts
    from PFASGroups.parser import preprocess_HalogenGroups

    custom_paths = get_componentSMARTSs()
    custom_paths = preprocess_componentsSmarts(custom_paths)
    custom_groups = get_HalogenGroups()
    pfas_groups, agg_groups = preprocess_HalogenGroups(custom_groups)
    results = parse_smiles(
        ["C(C(F)(F)F)F"],
        componentSmartss=custom_paths,
        pfas_groups=pfas_groups,
        agg_pfas_groups= agg_groups
    )
    assert len(results) == 1


# ===========================================================================
# README — Custom configuration: extend with HalogenGroup
# ===========================================================================


def test_readme_extend_with_custom_group():
    """README: append a custom HalogenGroup and parse with it."""
    from PFASGroups import get_compiled_HalogenGroups, HalogenGroup, parse_smiles

    groups = get_compiled_HalogenGroups()
    groups.append(HalogenGroup(
        id=999,
        name="My Custom Group",
        smarts={"[C](F)(F)F": 1},
        componentSmarts=None,
        componentSaturation="per",
        componentHalogens="F",
        componentForm="alkyl",
        constraints={"gte": {"F": 1}},
    ))
    results = parse_smiles(["FC(F)(F)C(F)(F)[N+](=O)[O-]"], pfas_groups=groups)
    assert len(results) == 1


def test_readme_custom_max_dist_from_comp():
    """README: HalogenGroup with custom max_dist_from_comp=3."""
    from PFASGroups import get_HalogenGroups, HalogenGroup, parse_smiles

    groups = get_HalogenGroups()
    customised_group = groups[1]
    customised_group['smarts'] = {x:2 for x,_ in customised_group['smarts'].items()}
    groups[1] = customised_group
    groups = [HalogenGroup(**x) for x in groups]
    groups.append(HalogenGroup(
        id=998,
        name="Extended Distance Group",
        smarts={"[#6$([#6][OH1])]": 1},
        componentSmarts=None,
        constraints={},
        max_dist_from_comp=3,
    ))
    results = parse_smiles(["OCC(F)(F)C(F)(F)C(F)(F)F"], pfas_groups=groups)
    assert len(results) == 1


# ===========================================================================
# README — Custom path types
# ===========================================================================


def test_readme_custom_path_types():
    """README: adds a Perchlorinated path type."""
    from PFASGroups import get_componentSMARTSs, parse_smiles
    from PFASGroups.core import preprocess_componentsSmarts

    paths = get_componentSMARTSs()
    paths['Cl']['alkyl'] = {"per":{
                    "component":"[C;X4](Cl)(Cl)!@!=!#[C;X4](Cl)(Cl)",
                    "name":"Custom Chloro"}}
    paths = preprocess_componentsSmarts(paths)

    results = parse_smiles(["ClC(Cl)(Cl)C(Cl)(Cl)C(=O)O"], halogens = 'Cl', componentSmartss=paths)
    assert len(results) == 1


# ===========================================================================
# README — to_sql persistence
# ===========================================================================


def test_readme_to_sql(tmp_path, results_fp):
    """README: results.to_sql(filename=...) creates a non-empty SQLite file."""
    db_path = tmp_path / "results.db"
    results_fp.to_sql(filename=str(db_path))
    assert db_path.exists()
    assert db_path.stat().st_size > 0


# ===========================================================================
# docs/fingerprints.rst — generate_fingerprint
# ===========================================================================


def test_docs_generate_fingerprint_preset_best():
    """fingerprints.rst: generate_fingerprint(smiles, preset='best') → (array, col_names)."""
    from PFASGroups import generate_fingerprint

    fp, col_names = generate_fingerprint(SMILES_FP_DOCS, preset="best")
    assert fp.ndim == 2
    assert fp.shape[0] == len(SMILES_FP_DOCS)
    assert fp.shape[1] == len(col_names)
    assert fp.shape[1] > 0


def test_docs_generate_fingerprint_single_returns_1d():
    """fingerprints.rst: single SMILES string is squeezed to 1-D."""
    from PFASGroups import generate_fingerprint

    fp, cols = generate_fingerprint("FC(F)(F)C(F)(F)C(=O)O", preset="best")
    assert fp.ndim == 1
    assert fp.shape[0] == len(cols)


def test_docs_tanimoto_between_fingerprints():
    """fingerprints.rst: Tanimoto similarity lies in [0, 1]."""
    from PFASGroups import generate_fingerprint

    fp, _ = generate_fingerprint(SMILES_FP_DOCS, preset="best")

    def tanimoto(a, b):
        a, b = a.astype(float), b.astype(float)
        denom = np.dot(a, a) + np.dot(b, b) - np.dot(a, b)
        return float(np.dot(a, b) / denom) if denom > 0 else 0.0

    t = tanimoto(fp[0], fp[1])
    assert 0.0 <= t <= 1.0 + 1e-9


def test_docs_generate_fingerprint_expand_preset():
    """fingerprints.rst: expand preset manually and pass component_metrics directly."""
    from PFASGroups import generate_fingerprint, FINGERPRINT_PRESETS

    p = FINGERPRINT_PRESETS["best_3"]
    fp, cols = generate_fingerprint(SMILES_FP_DOCS, component_metrics=p["component_metrics"])
    assert fp.shape[0] == len(SMILES_FP_DOCS)
    assert fp.shape[1] == len(cols)


# ===========================================================================
# docs/fingerprints.rst — to_fingerprint (deprecated wrapper)
# ===========================================================================


def test_docs_to_fingerprint_deprecated_warning(results_fp):
    """fingerprints.rst: results.to_fingerprint() emits DeprecationWarning."""
    import warnings

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        arr = results_fp.to_fingerprint(preset="best")
        assert any(issubclass(x.category, DeprecationWarning) for x in w)
    assert arr.shape[0] == len(SMILES_FP)


# ===========================================================================
# docs/getting_started.rst — output_format variants
# ===========================================================================


def test_docs_output_format_dataframe():
    """getting_started.rst: output_format='dataframe' returns a pandas DataFrame."""
    pytest.importorskip("pandas")
    import pandas as pd
    from PFASGroups import parse_smiles

    df = parse_smiles(["FC(F)(F)C(F)(F)C(=O)O", "C(C(F)(F)F)F"], output_format="dataframe")
    assert isinstance(df, pd.DataFrame)


def test_docs_output_format_csv():
    """getting_started.rst: output_format='csv' returns a non-empty CSV string."""
    from PFASGroups import parse_smiles

    csv_str = parse_smiles(["FC(F)(F)C(F)(F)C(=O)O", "C(C(F)(F)F)F"], output_format="csv")
    assert isinstance(csv_str, str)
    assert "\n" in csv_str


# ===========================================================================
# README — Dimensionality reduction (plot=False to suppress display)
# ===========================================================================


def test_readme_perform_pca(results_fp):
    """README: perform_pca(n_components=2, plot=False) → dict with 'transformed'."""
    pytest.importorskip("sklearn")
    pytest.importorskip("matplotlib")

    result = results_fp.perform_pca(n_components=2, plot=False)
    assert "transformed" in result
    assert result["transformed"].shape == (len(SMILES_FP), 2)


def test_readme_perform_pca_color_by_top_group(results_fp):
    """README: perform_pca(color_by='top_group') does not raise."""
    pytest.importorskip("sklearn")
    pytest.importorskip("matplotlib")

    result = results_fp.perform_pca(n_components=2, plot=False, color_by="top_group")
    assert "transformed" in result


def test_readme_perform_pca_color_by_labels(results_fp):
    """README: perform_pca(color_by=custom_labels) does not raise."""
    pytest.importorskip("sklearn")
    pytest.importorskip("matplotlib")

    labels = ["group_A", "group_A", "group_B", "group_B"]
    result = results_fp.perform_pca(n_components=2, plot=False, color_by=labels)
    assert "transformed" in result


def test_readme_perform_tsne(results_fp):
    """README: perform_tsne(perplexity=3, plot=False) → dict with 'transformed'."""
    pytest.importorskip("sklearn")
    pytest.importorskip("matplotlib")

    # perplexity must be < n_samples; use 3 for 4-molecule fixture
    result = results_fp.perform_tsne(perplexity=3, plot=False)
    assert "transformed" in result
    assert result["transformed"].shape[0] == len(SMILES_FP)


def test_readme_perform_umap(results_fp):
    """README: perform_umap(n_neighbors=2, plot=False) → dict with 'transformed'."""
    pytest.importorskip("umap")
    pytest.importorskip("matplotlib")

    result = results_fp.perform_umap(n_neighbors=2, plot=False)
    assert "transformed" in result
    assert result["transformed"].shape[0] == len(SMILES_FP)
