#!/usr/bin/env python3
"""
2-D chemical space visualisation of OECD PFAS compounds.

Supports three reduction methods (--method):
    tsne   (default)  t-SNE via scikit-learn
    mds               Metric MDS via scikit-learn  (slow for > 2000 pts)
    umap              UMAP (requires umap-learn)

Fingerprint presets (--fp):
    full   (default)  EGR + n_spacer + ring_size + mol metrics  (~15 min first run)
    total    total + n_spacer + ring_size + mol metrics  (~15 min first run)
    binary            binary presence/absence only (< 2 min, good for exploration)

Categories are collapsed to the N most frequent groups (--top N, default 9).
The remainder are merged into "Other PFAS".
Convex hulls are drawn for every category with >= 4 points.
Fingerprints are cached to avoid recomputation on subsequent runs.

Output
------
    benchmark/imgs/<method>_oecd_pfas[_binary].png
    benchmark/data/<method>_oecd_pfas[_binary].csv

Usage
-----
    conda activate chem
    cd benchmark
    python scripts/umap_oecd_pfas.py --method tsne --fp binary   # quick
    python scripts/umap_oecd_pfas.py --method tsne               # full (cached)
    python scripts/umap_oecd_pfas.py --method mds --fp binary
    python scripts/umap_oecd_pfas.py --method umap --top 12
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR   = SCRIPT_DIR.parents[2]
DATA_DIR   = SCRIPT_DIR.parents[1] / "data"
TEST_DATA  = SCRIPT_DIR.parents[1] / "test_data"
IMGS_DIR   = SCRIPT_DIR.parents[1] / "imgs"

OECD_CSV     = TEST_DATA / "S25_OECDPFAS_list_22012019.csv"
OTHER_LABEL  = "Other PFAS"
HULL_MIN_PTS = 4

sys.path.insert(0, str(ROOT_DIR))
IMGS_DIR.mkdir(exist_ok=True)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _short_category(cat):
    """Truncate very long category labels for legend readability."""
    return cat if len(cat) <= 40 else cat[:38] + "..."


def _collapse_categories(categories, top_n):
    """Keep the top_n most frequent categories; relabel the rest as OTHER_LABEL."""
    top = {cat for cat, _ in Counter(categories).most_common(top_n)}
    return [c if c in top else OTHER_LABEL for c in categories]


def _draw_convex_hulls(ax, embedding_2d, categories, cat_to_color, min_pts=HULL_MIN_PTS):
    """Draw a filled convex hull outline for each category that has enough points."""
    from scipy.spatial import ConvexHull, QhullError
    for cat in sorted(set(categories)):
        mask = np.array([c == cat for c in categories])
        pts  = embedding_2d[mask]
        if len(pts) < min_pts:
            continue
        try:
            hull = ConvexHull(pts)
        except QhullError:
            continue
        verts = np.vstack([pts[hull.vertices], pts[hull.vertices[0]]])
        color = cat_to_color[cat]
        ax.fill(verts[:, 0], verts[:, 1], color=color, alpha=0.12, zorder=1)
        ax.plot(verts[:, 0], verts[:, 1], color=color, alpha=0.50,
                linewidth=0.9, zorder=2)


def _reduce_2d(X, method):
    """Reduce feature matrix X to 2-D using the requested method."""
    if method == "tsne":
        from sklearn.manifold import TSNE
        n_pca = min(50, X.shape[1], X.shape[0] - 1)
        if X.shape[1] > n_pca:
            from sklearn.decomposition import PCA
            print("  PCA pre-reduction to %d components ..." % n_pca)
            X = PCA(n_components=n_pca, random_state=42).fit_transform(X)
        reducer = TSNE(
            n_components=2,
            perplexity=min(30, len(X) - 1),
            learning_rate="auto",
            init="pca",
            random_state=42,
            max_iter=1000,
            verbose=1,
        )
        return reducer.fit_transform(X)

    elif method == "mds":
        from sklearn.manifold import MDS
        reducer = MDS(
            n_components=2,
            metric=True,
            n_init=4,
            max_iter=300,
            random_state=42,
            dissimilarity="euclidean",
            n_jobs=-1,
            verbose=1,
        )
        return reducer.fit_transform(X)

    elif method == "umap":
        try:
            from umap import UMAP
        except ImportError:
            raise ImportError(
                "umap-learn is required: conda install -c conda-forge umap-learn"
            )
        reducer = UMAP(
            n_components=2,
            n_neighbors=15,
            min_dist=0.1,
            metric="euclidean",
            random_state=42,
            verbose=True,
        )
        return reducer.fit_transform(X)

    else:
        raise ValueError("Unknown method %r. Choose tsne, mds, or umap." % method)


def _load_or_compute_fingerprint(fp_preset):
    """Return (X_nz, smiles_nz, names_nz, raw_categories) from cache or computed."""
    cache_X    = DATA_DIR / ("oecd_%s_fp.npy" % fp_preset)
    cache_meta = DATA_DIR / ("oecd_%s_meta.csv" % fp_preset)

    if cache_X.exists() and cache_meta.exists():
        print("  Loading cached fingerprint from %s ..." % cache_X.name)
        X_nz = np.load(cache_X)
        meta = pd.read_csv(cache_meta)
        return (
            X_nz,
            meta["smiles"].tolist(),
            meta["name"].tolist(),
            meta["raw_category"].tolist(),
        )

    # --- Compute from scratch ---
    print("Loading OECD PFAS list from %s ..." % OECD_CSV.name)
    oecd_df = pd.read_csv(OECD_CSV)
    smiles_col, name_col = "SMILES", "PREFERRED_NAME"
    valid = oecd_df[smiles_col].notna() & (oecd_df[smiles_col].str.strip() != "-")
    oecd_valid  = oecd_df[valid].reset_index(drop=True)
    smiles_list = oecd_valid[smiles_col].tolist()
    names_list  = oecd_valid[name_col].fillna("").tolist()
    print("  Valid SMILES: %d" % len(smiles_list))

    print("\nParsing with PFASGroups ...")
    from PFASGroups import parse_smiles
    from PFASGroups.getter import get_compiled_HalogenGroups
    results     = parse_smiles(smiles_list, halogens="F", progress=True)
    pfas_groups = get_compiled_HalogenGroups(halogens="F")
    print("  Parsed %d molecules." % len(results))

    parsed_smiles = [emb.get("smiles", smiles_list[i]) for i, emb in enumerate(results)]
    parsed_names  = names_list[: len(results)]

    print("\nComputing '%s' fingerprint ..." % fp_preset)
    if fp_preset == "binary":
        X = np.asarray(
            results.to_array(
                component_metrics=["binary"],
                pfas_groups=pfas_groups,
                progress=True,
            )
        )
    elif  fp_preset == 'total':  # total_component + n_spacer + ring_size + mol metrics
        mol_metrics = [
            "n_components", "total_size", "mean_size", "max_size",
            "mean_branching", "max_branching",
            "mean_eccentricity", "max_diameter",
            "mean_component_fraction", "max_component_fraction",
        ]
        X = np.asarray(
            results.to_array(
                component_metrics=["total_component", "n_spacer", "ring_size"],
                molecule_metrics=mol_metrics,
                pfas_groups=pfas_groups,
                progress=True,
            ))
    else:  # full: EGR + n_spacer + ring_size + mol metrics
        mol_metrics = [
            "n_components", "total_size", "mean_size", "max_size",
            "mean_branching", "max_branching",
            "mean_eccentricity", "max_diameter",
            "mean_component_fraction", "max_component_fraction",
        ]
        X = np.asarray(
            results.to_array(
                component_metrics=["effective_graph_resistance", "n_spacer", "ring_size"],
                molecule_metrics=mol_metrics,
                pfas_groups=pfas_groups,
                progress=True,
            )
        )
    print("  Fingerprint shape: %s" % str(X.shape))

    nonzero_mask = np.abs(X).sum(axis=1) > 0
    X_nz      = X[nonzero_mask]
    smiles_nz = [parsed_smiles[i] for i in range(len(X)) if nonzero_mask[i]]
    names_nz  = [parsed_names[i]  for i in range(len(X)) if nonzero_mask[i]]
    print("  Non-zero (matched): %d / %d" % (len(X_nz), len(X)))

    print("\nClassifying molecules ...")
    raw_categories = []
    for i, emb in enumerate(results):
        if not nonzero_mask[i]:
            continue
        cat, _ = emb.classify()
        raw_categories.append(_short_category(cat))

    # Save cache
    np.save(cache_X, X_nz)
    meta_df = pd.DataFrame({
        "smiles":       smiles_nz,
        "name":         names_nz,
        "raw_category": raw_categories,
    })
    meta_df.to_csv(cache_meta, index=False)
    print("  [cached] %s  (%s)" % (cache_X.name, cache_meta.name))

    return X_nz, smiles_nz, names_nz, raw_categories


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="2-D chemical space of OECD PFAS")
    parser.add_argument(
        "--method", choices=["tsne", "mds", "umap"], default="tsne",
        help="Dimensionality-reduction method (default: tsne)",
    )
    parser.add_argument(
        "--fp", choices=["full", "binary", "total"], default="full",
        help="Fingerprint: 'full' (EGR+spacer+ring+mol), 'total' (total_component+spacer+ring+mol), or 'binary' (fast)",
    )
    parser.add_argument(
        "--top", type=int, default=9,
        help="Top N categories to show; rest -> 'Other PFAS' (default: 9)",
    )
    args   = parser.parse_args()
    method = args.method
    fp     = args.fp
    top_n  = args.top

    # --- Load / compute fingerprint ---
    X_nz, smiles_nz, names_nz, raw_categories = _load_or_compute_fingerprint(fp)

    # --- Collapse to top_n categories ---
    categories  = _collapse_categories(raw_categories, top_n)
    unique_cats = sorted(set(categories))
    print("\nCategories (top %d):" % top_n)
    for cat, cnt in Counter(categories).most_common():
        print("  %5d  %s" % (cnt, cat))

    # --- Standardise features ---
    from sklearn.preprocessing import StandardScaler
    X_scaled = StandardScaler().fit_transform(X_nz)

    # --- Dimensionality reduction ---
    print("\nRunning %s (n=%d) ..." % (method.upper(), len(X_scaled)))
    embedding_2d = _reduce_2d(X_scaled, method)
    print("  Output shape: %s" % str(embedding_2d.shape))

    # --- Build colour map (tab10 qualitative palette) ---
    freq_order  = [c for c, _ in Counter(categories).most_common() if c != OTHER_LABEL]
    sorted_cats = freq_order + ([OTHER_LABEL] if OTHER_LABEL in unique_cats else [])
    tab10 = plt.cm.tab10.colors
    cat_to_color = {}
    idx = 0
    for cat in sorted_cats:
        if cat == OTHER_LABEL:
            cat_to_color[cat] = "#aaaaaa"
        else:
            cat_to_color[cat] = tab10[idx % len(tab10)]
            idx += 1
    point_colors = [cat_to_color[c] for c in categories]

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(11, 8))
    _draw_convex_hulls(ax, embedding_2d, categories, cat_to_color)
    ax.scatter(
        embedding_2d[:, 0], embedding_2d[:, 1],
        c=point_colors, s=18, alpha=0.82,
        edgecolors="none", linewidths=0, zorder=3,
    )
    method_label = {"tsne": "t-SNE", "mds": "MDS", "umap": "UMAP"}[method]
    fp_label     = {"full": "EGR+spacer+ring+mol", "binary": "binary", "total": "total_component+spacer+ring+mol"}[fp]
    ax.set_xlabel("%s 1" % method_label, fontsize=12)
    ax.set_ylabel("%s 2" % method_label, fontsize=12)
    ax.set_title(
        "OECD PFAS -- chemical space  (%s, %s)" % (method_label, fp_label),
        fontsize=13, fontweight="bold",
    )
    handles = [
        mpatches.Patch(facecolor=cat_to_color[c], label=c,
                       edgecolor="#555", linewidth=0.8)
        for c in sorted_cats
    ]
    ax.legend(
        handles=handles, title="PFASGroups category",
        loc="upper left", bbox_to_anchor=(1.01, 1.0),
        borderaxespad=0, frameon=True, fontsize=9, title_fontsize=10,
    )
    ax.grid(True, linewidth=0.3, alpha=0.4, zorder=0)
    plt.tight_layout()

    suffix   = "" if fp == "full" else ("_%s" % fp)
    out_stem = "%s_oecd_pfas%s" % (method, suffix)
    out_img  = IMGS_DIR / ("%s.png" % out_stem)
    fig.savefig(out_img, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print("\n[saved] %s" % out_img.relative_to(IMGS_DIR.parent))

    out_csv = DATA_DIR / ("%s.csv" % out_stem)
    pd.DataFrame({
        "smiles":           smiles_nz,
        "name":             names_nz,
        "category":         categories,
        "raw_category":     raw_categories,
        ("%s_1" % method):  embedding_2d[:, 0],
        ("%s_2" % method):  embedding_2d[:, 1],
    }).to_csv(out_csv, index=False)
    print("[saved] %s" % out_csv.relative_to(DATA_DIR.parent))
    print("\nDone.")


if __name__ == "__main__":
    main()
