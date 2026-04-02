"""
gui/utils/chemspace.py
─────────────────────
Builds a self-contained Plotly HTML string for the chemical space visualisation.

Supports UMAP, PCA, and t-SNE projections of PFASGroups embedding arrays.
Hover tooltip includes molecule name, SMILES, and (optionally) a user label.
"""
from __future__ import annotations

import base64
from typing import Any, Callable, Optional

import numpy as np


def build_chemspace_html(
    embedding_set,
    method: str = "UMAP",
    preset: str = "best",
    label_col: Optional[str] = None,
    df_labels=None,          # pd.DataFrame with a column matching label_col
    n_neighbours: int = 15,
    min_dist: float = 0.1,
    perplexity: float = 30.0,
    progress_cb: Optional[Callable[[int], None]] = None,
) -> str:
    """Compute 2-D projection and return a full Plotly HTML page.

    Parameters
    ----------
    embedding_set : PFASEmbeddingSet
    method : "UMAP" | "PCA" | "t-SNE"
    preset : key from FINGERPRINT_PRESETS
    label_col : name of label column in df_labels (used for colour)
    df_labels : DataFrame aligned with embedding_set (same length)
    n_neighbours / min_dist : UMAP hyper-parameters
    perplexity : t-SNE hyper-parameter
    progress_cb : optional int callback (0-100)
    """
    import plotly.graph_objects as go
    from sklearn.preprocessing import StandardScaler

    if progress_cb:
        progress_cb(5)

    # ── Feature matrix ────────────────────────────────────────────────────
    arr = embedding_set.to_array(preset=preset, progress=False)
    X = np.asarray(arr, dtype=float)

    # Replace NaN/inf
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)

    X = StandardScaler().fit_transform(X)

    if progress_cb:
        progress_cb(25)

    # ── Dimensionality reduction ──────────────────────────────────────────
    method = method.upper()
    if method == "UMAP":
        try:
            from umap import UMAP
        except ImportError as exc:
            raise ImportError("umap-learn is not installed.") from exc
        reducer = UMAP(
            n_components=2,
            n_neighbors=n_neighbours,
            min_dist=min_dist,
            random_state=42,
        )
        X2 = reducer.fit_transform(X)
    elif method == "PCA":
        from sklearn.decomposition import PCA
        X2 = PCA(n_components=2, random_state=42).fit_transform(X)
    elif method in ("T-SNE", "TSNE", "T_SNE"):
        from sklearn.manifold import TSNE
        X2 = TSNE(
            n_components=2,
            perplexity=min(perplexity, max(5, len(X) - 1)),
            random_state=42,
        ).fit_transform(X)
    else:
        raise ValueError(f"Unknown method: {method!r}")

    if progress_cb:
        progress_cb(70)

    # ── Labels and hover ─────────────────────────────────────────────────
    names = [str(m.get("name") or m.get("smiles") or i)
             for i, m in enumerate(embedding_set)]
    smiles_list = [str(m.get("smiles") or "") for m in embedding_set]

    if label_col and df_labels is not None and label_col in df_labels.columns:
        colours = df_labels[label_col].astype(str).tolist()
    else:
        colours = [_top_group(m) for m in embedding_set]

    hover_text = [
        f"<b>{n}</b><br>SMILES: {s}<br>Label: {c}"
        for n, s, c in zip(names, smiles_list, colours)
    ]

    # ── Build figure ─────────────────────────────────────────────────────
    unique_colours = sorted(set(colours))
    colour_map = _make_colour_map(unique_colours)

    traces = []
    for colour in unique_colours:
        mask = [c == colour for c in colours]
        idx = [i for i, m in enumerate(mask) if m]
        traces.append(go.Scatter(
            x=X2[idx, 0],
            y=X2[idx, 1],
            mode="markers",
            name=colour,
            marker=dict(color=colour_map[colour], size=8, opacity=0.8),
            text=[hover_text[i] for i in idx],
            hovertemplate="%{text}<extra></extra>",
            customdata=[smiles_list[i] for i in idx],
        ))

    fig = go.Figure(
        data=traces,
        layout=go.Layout(
            title=f"{method} chemical space (preset: {preset})",
            xaxis_title=f"{method}-1",
            yaxis_title=f"{method}-2",
            legend_title="Label",
            hovermode="closest",
            paper_bgcolor="white",
            plot_bgcolor="#fafafa",
            font=dict(family="Inter, Arial, sans-serif"),
            margin=dict(l=40, r=20, t=50, b=40),
        ),
    )

    if progress_cb:
        progress_cb(90)

    html = fig.to_html(full_html=True, include_plotlyjs="cdn")

    if progress_cb:
        progress_cb(100)

    return html


# ── Helpers ───────────────────────────────────────────────────────────────

def _top_group(embedding) -> str:
    """Return the name of the most frequently matched PFAS group, or 'None'."""
    try:
        groups = embedding.get("groups") or []
        if not groups:
            return "None"
        best = max(groups, key=lambda g: g.get("n_components", 0))
        return str(best.get("name") or "None")
    except Exception:
        return "None"


def _make_colour_map(labels: list[str]) -> dict[str, str]:
    """Map label strings to hex colours using a qualitative palette."""
    palette = [
        "#306DBA", "#E15D0B", "#9D206C", "#51127C",
        "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
        "#66A61E", "#E6AB02", "#A6761D", "#666666",
    ]
    return {lbl: palette[i % len(palette)] for i, lbl in enumerate(labels)}
