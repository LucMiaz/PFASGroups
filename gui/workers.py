"""
QThread worker base classes for all long-running operations.

Every worker emits:
  progress(int)        – 0-100 completion percentage
  status(str)          – short human-readable status message
  result(object)       – arbitrary payload delivered on success
  error(str)           – error message on failure

Usage
-----
    worker = ClassifyWorker(smiles_list, options)
    worker.progress.connect(progress_bar.setValue)
    worker.status.connect(status_label.setText)
    worker.result.connect(handle_result)
    worker.error.connect(show_error)
    worker.finished.connect(worker.deleteLater)
    worker.start()
"""
from __future__ import annotations

from typing import Any

from PySide6.QtCore import QThread, Signal


# ---------------------------------------------------------------------------
# Base worker
# ---------------------------------------------------------------------------

class Worker(QThread):
    """Abstract base class for all background workers."""

    progress = Signal(int)      # 0-100
    status   = Signal(str)      # human-readable message
    result   = Signal(object)   # arbitrary payload
    error    = Signal(str)      # error message (empty string = cancelled)

    def run(self):
        raise NotImplementedError


# ---------------------------------------------------------------------------
# Classification worker
# ---------------------------------------------------------------------------

class ClassifyWorker(Worker):
    """Parse a list of SMILES strings through PFASGroups."""

    def __init__(
        self,
        smiles: list[str],
        names: list[str] | None = None,
        halogens: str | list[str] = "F",
        compute_metrics: bool = True,
        check_definitions: bool = True,
    ):
        super().__init__()
        self._smiles = smiles
        self._names = names or []
        self._halogens = halogens
        self._compute_metrics = compute_metrics
        self._check_definitions = check_definitions

    def run(self):
        try:
            from PFASGroups import parse_smiles
            self.status.emit("Loading PFASGroups…")
            self.progress.emit(2)

            total = len(self._smiles)
            if total == 0:
                self.error.emit("No SMILES strings provided.")
                return

            results = []
            for i, smi in enumerate(self._smiles):
                pct = int(5 + 90 * i / total)
                self.progress.emit(pct)
                self.status.emit(f"Processing molecule {i + 1}/{total}…")
                try:
                    res = parse_smiles(
                        [smi],
                        halogens=self._halogens,
                        compute_component_metrics=self._compute_metrics,
                    )
                    results.append(res)
                except Exception as exc:
                    results.append(None)

            self.status.emit("Combining results…")
            self.progress.emit(97)

            # Merge individual PFASEmbeddingSet objects
            from PFASGroups import parse_smiles as _ps
            merged = _ps(
                self._smiles,
                halogens=self._halogens,
                compute_component_metrics=self._compute_metrics,
            )
            if self._names:
                for i, emb in enumerate(merged):
                    if i < len(self._names):
                        emb["name"] = self._names[i]

            self.progress.emit(100)
            self.status.emit("Done.")
            self.result.emit(merged)

        except Exception as exc:
            self.error.emit(f"Classification failed: {exc}")


# ---------------------------------------------------------------------------
# Definition-test worker
# ---------------------------------------------------------------------------

class DefinitionTestWorker(Worker):
    """Test a single molecule against a set of PFASDefinitions."""

    def __init__(self, smiles: str, custom_definitions: list | None = None):
        super().__init__()
        self._smiles = smiles
        self._custom_defs = custom_definitions  # list of PFASDefinition objects

    def run(self):
        try:
            from rdkit import Chem
            from PFASGroups.parser import parse_definitions_in_mol
            from PFASGroups import get_PFASDefinitions

            self.status.emit("Parsing molecule…")
            self.progress.emit(10)

            mol = Chem.MolFromSmiles(self._smiles)
            if mol is None:
                self.error.emit(f"Invalid SMILES: {self._smiles!r}")
                return

            self.status.emit("Loading definitions…")
            self.progress.emit(30)

            definitions = self._custom_defs if self._custom_defs else get_PFASDefinitions()

            self.status.emit("Testing definitions…")
            self.progress.emit(50)

            matched = parse_definitions_in_mol(mol)
            matched_ids = {d.id for d in matched}

            results = []
            for d in definitions:
                results.append({
                    "id": d.id,
                    "name": d.name,
                    "description": d.description,
                    "passed": d.id in matched_ids,
                })

            self.progress.emit(100)
            self.status.emit("Done.")
            self.result.emit(results)

        except Exception as exc:
            self.error.emit(f"Definition test failed: {exc}")


# ---------------------------------------------------------------------------
# Custom SMARTS/ComponentSmarts tester worker
# ---------------------------------------------------------------------------

class SmartsTestWorker(Worker):
    """Test a custom HalogenGroup definition against a SMILES string."""

    def __init__(self, smiles: str, group_dict: dict):
        super().__init__()
        self._smiles = smiles
        self._group_dict = group_dict

    def run(self):
        try:
            from rdkit import Chem
            from PFASGroups import HalogenGroup

            self.status.emit("Parsing molecule…")
            self.progress.emit(20)

            mol = Chem.MolFromSmiles(self._smiles)
            if mol is None:
                self.error.emit(f"Invalid SMILES: {self._smiles!r}")
                return

            self.status.emit("Building group…")
            self.progress.emit(40)

            try:
                group = HalogenGroup(**self._group_dict)
            except Exception as exc:
                self.error.emit(f"Invalid group definition: {exc}")
                return

            self.status.emit("Running test…")
            self.progress.emit(60)

            test_result = group.test(mol)

            self.progress.emit(100)
            self.status.emit("Done.")
            self.result.emit(test_result)

        except Exception as exc:
            self.error.emit(f"SMARTS test failed: {exc}")


# ---------------------------------------------------------------------------
# Prioritisation worker
# ---------------------------------------------------------------------------

class PrioritiseWorker(Worker):
    """Run prioritise_molecules on a PFASEmbeddingSet."""

    def __init__(
        self,
        embedding_set,
        reference=None,
        mode: str = "intrinsic",
        a: float = 1.0,
        b: float = 1.0,
        percentile: float = 90.0,
    ):
        super().__init__()
        self._embs = embedding_set
        self._reference = reference
        self._mode = mode
        self._a = a
        self._b = b
        self._percentile = percentile

    def run(self):
        try:
            from PFASGroups import prioritise_molecules, get_priority_statistics

            self.status.emit("Computing priority scores…")
            self.progress.emit(20)

            kwargs = dict(
                a=self._a,
                b=self._b,
                percentile=self._percentile,
                return_scores=True,
            )
            if self._reference is not None:
                kwargs["reference"] = self._reference

            ranked, scores = prioritise_molecules(self._embs, **kwargs)

            self.status.emit("Computing statistics…")
            self.progress.emit(80)

            stats = get_priority_statistics(ranked)

            self.progress.emit(100)
            self.status.emit("Done.")
            self.result.emit({"ranked": ranked, "scores": scores, "stats": stats})

        except Exception as exc:
            self.error.emit(f"Prioritisation failed: {exc}")


# ---------------------------------------------------------------------------
# Chemical space worker
# ---------------------------------------------------------------------------

class ChemSpaceWorker(Worker):
    """Compute dimensionality reduction and return a Plotly HTML string."""

    def __init__(
        self,
        embedding_set,
        method: str = "UMAP",
        preset: str = "best",
        label_col: str | None = None,
        df_labels=None,       # optional extra DataFrame with label_col
        n_neighbours: int = 15,
        min_dist: float = 0.1,
    ):
        super().__init__()
        self._embs = embedding_set
        self._method = method
        self._preset = preset
        self._label_col = label_col
        self._df_labels = df_labels
        self._n_neighbours = n_neighbours
        self._min_dist = min_dist

    def run(self):
        try:
            from gui.utils.chemspace import build_chemspace_html

            self.status.emit("Building fingerprint matrix…")
            self.progress.emit(15)

            html = build_chemspace_html(
                self._embs,
                method=self._method,
                preset=self._preset,
                label_col=self._label_col,
                df_labels=self._df_labels,
                n_neighbours=self._n_neighbours,
                min_dist=self._min_dist,
                progress_cb=lambda p, s: (self.progress.emit(p), self.status.emit(s)),
            )

            self.progress.emit(100)
            self.status.emit("Done.")
            self.result.emit(html)

        except Exception as exc:
            self.error.emit(f"Chemical space computation failed: {exc}")


# ---------------------------------------------------------------------------
# Modelling worker
# ---------------------------------------------------------------------------

class ModelWorker(Worker):
    """Run ML benchmark (HistGBM + Bayesian t-test) in background."""

    def __init__(
        self,
        feature_sets: dict,   # name → np.ndarray
        y,                    # 1-D array of target labels
        cv_splits: int = 5,
        cv_repeats: int = 3,
    ):
        super().__init__()
        self._feature_sets = feature_sets
        self._y = y
        self._cv_splits = cv_splits
        self._cv_repeats = cv_repeats

    def run(self):
        try:
            from gui.utils.modelling import run_benchmark

            self.status.emit("Running cross-validation…")
            self.progress.emit(10)

            results = run_benchmark(
                self._feature_sets,
                self._y,
                cv_splits=self._cv_splits,
                cv_repeats=self._cv_repeats,
                progress_cb=lambda p, s: (self.progress.emit(p), self.status.emit(s)),
            )

            self.progress.emit(100)
            self.status.emit("Done.")
            self.result.emit(results)

        except Exception as exc:
            self.error.emit(f"Modelling failed: {exc}")
