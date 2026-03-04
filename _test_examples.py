"""Test all article examples to verify they work with the current API."""
import sys
errors = []

# ── Example 1: Custom group ────────────────────────────────────────────────────
print("=== Example 1: Custom group ===")
try:
    from PFASGroups import get_compiled_HalogenGroups, HalogenGroup, parse_smiles
    groups = get_compiled_HalogenGroups()
    print(f"  get_compiled_HalogenGroups() → {len(groups)} groups")
    groups.append(HalogenGroup(
        id=200,
        name="Perfluoroalkyl nitrates",
        smarts={"[C$(C[ON+](=O)[O-])]": 1},
        componentSaturation="per",
        componentHalogens="F",
        componentForm="alkyl",
        constraints={"eq": {"N": 1}, "gte": {"F": 1}},
    ))
    results = parse_smiles(["FC(F)(F)C(F)(F)ON(=O)=O"], pfas_groups=groups)
    print(f"  parse_smiles → {type(results).__name__}, len={len(results)}")
    print("  ✓ Example 1 OK")
except Exception as e:
    print(f"  ✗ Example 1 FAILED: {e}")
    errors.append(("Example 1", str(e)))

# ── Example 2: Fingerprint generation ─────────────────────────────────────────
print("\n=== Example 2: Fingerprint generation ===")
try:
    from PFASGroups import parse_smiles as ps, generate_fingerprint

    smiles = [
        "FC(F)(F)C(F)(F)C(F)(F)C(=O)O",        # PFBA
        "FC(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O",    # PFBS
        "C(CF)(CF)C(F)(F)C(F)(F)OCC(F)(F)F",   # FTOH-like
    ]

    fp_matrix, info = generate_fingerprint(smiles, halogens='F', saturation='per')
    print(f"  generate_fingerprint shape: {fp_matrix.shape}")
    print(f"  info keys: {list(info.keys())}")
    print(f"  group_names[:3]: {info['group_names'][:3]}")

    results = ps(smiles)
    fp = results.to_fingerprint(group_selection='all', halogens='F', count_mode='binary')
    print(f"  to_fingerprint repr: {repr(fp)}")

    pca = fp.perform_pca(n_components=2, plot=False)
    print(f"  perform_pca → {type(pca).__name__}")

    try:
        tsne = fp.perform_tsne(perplexity=2, max_iter=300, plot=False)
        print(f"  perform_tsne → {type(tsne).__name__}")
    except Exception as e:
        print(f"  perform_tsne: {e}")
        errors.append(("Example 2 (tsne)", str(e)))

    try:
        umap_res = fp.perform_umap(n_neighbors=2, plot=False)
        print(f"  perform_umap → {type(umap_res).__name__}")
    except Exception as e:
        print(f"  perform_umap (may need umap-learn): {e}")
        # not fatal

    other_fp = ps(["FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O"]).to_fingerprint()
    kld = fp.compare_kld(other_fp, method='minmax')
    print(f"  compare_kld → {type(kld).__name__}: {kld}")
    print("  ✓ Example 2 OK")
except Exception as e:
    import traceback; traceback.print_exc()
    print(f"  ✗ Example 2 FAILED: {e}")
    errors.append(("Example 2", str(e)))

# ── Example 3: Multi-halogen ───────────────────────────────────────────────────
print("\n=== Example 3: Multi-halogen ===")
try:
    from HalogenGroups import parse_smiles as hal_parse
    from PFASGroups import parse_smiles as pfas_parse

    multi = hal_parse(["ClC(Cl)(Cl)CCl", "BrC(Br)(Br)CBr"])
    pfas  = pfas_parse(["FC(F)(F)C(F)(F)C(=O)O"])
    print(f"  HalogenGroups parse → {len(multi)} results")
    print(f"  PFASGroups parse    → {len(pfas)} results")
    print("  ✓ Example 3 OK")
except Exception as e:
    print(f"  ✗ Example 3 FAILED: {e}")
    errors.append(("Example 3", str(e)))

# ── Example 4: Prioritization ─────────────────────────────────────────────────
print("\n=== Example 4: Prioritization ===")
try:
    from PFASGroups import parse_smiles as ps2, prioritise_molecules
    known_persistent_pfas = [
        "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # PFOA
        "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O",  # PFBS
    ]
    unknown_library = [
        "FC(F)(F)C(F)(F)C(=O)O",          # PFBA
        "FC(F)(F)C(F)(F)C(F)(F)C(=O)O",   # PFPA
        "CCO",                              # ethanol (not PFAS)
    ]
    reference  = ps2(known_persistent_pfas)
    candidates = ps2(unknown_library)
    ranked, scores = prioritise_molecules(candidates, reference=reference,
                                          group_selection='oecd',
                                          return_scores=True)
    print(f"  prioritise_molecules → {len(ranked)} ranked, {len(scores)} scores")
    print("  ✓ Example 4 OK")
except Exception as e:
    import traceback; traceback.print_exc()
    print(f"  ✗ Example 4 FAILED: {e}")
    errors.append(("Example 4", str(e)))

# ── Summary ───────────────────────────────────────────────────────────────────
print("\n=== Summary ===")
if errors:
    print(f"FAILED: {len(errors)} error(s)")
    for name, msg in errors:
        print(f"  {name}: {msg}")
    sys.exit(1)
else:
    print("All examples passed ✓")
