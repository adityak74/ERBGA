"""JSON-in, JSON-out functions for agents, the CLI and the MCP server.

Every public function takes plain values (str, int, float, bool, list, dict)
and returns a dict that ``json.dumps`` can serialise. Failures raise
``ERBGAError``; ``error.to_dict()`` gives
``{"ok": false, "error": {"code": ..., "message": ..., "hint": ...}}``.

A *graph spec* can be any of:

- a bundled dataset name: ``"karate"``, ``"dolphins"``, ``"polbooks"``, ``"football"``
- a file path: ``.gml``, ``.json`` (see below), or an edge list (``.txt``,
  ``.csv``, ``.tsv``, ``.edgelist``; one ``u v`` pair per line, ``#`` comments)
- an inline dict, or the same as a JSON string: ``{"edges": [[u, v], ...],
  "nodes": [...]}`` (``nodes`` is optional), or networkx node-link JSON
  ``{"nodes": [{"id": ...}], "links": [{"source": ..., "target": ...}]}``
- a ``Graph`` object (Python callers only)
"""

from __future__ import annotations

import json
import re
from dataclasses import asdict, fields
from importlib import resources
from pathlib import Path
from typing import Any

import numpy as np

from .ga import ERBGA, GAConfig
from .graph import Graph, read_gml
from .objectives import OBJECTIVES, get_objective
from .partition import canonical_cut, decode

BUNDLED_DATASETS = ("karate", "dolphins", "polbooks", "football")
CONFIG_FIELDS = tuple(f.name for f in fields(GAConfig))


class ERBGAError(Exception):
    """An error with a stable machine-readable ``code`` and a fix-it ``hint``.

    Codes: GRAPH_NOT_FOUND, INVALID_GRAPH, EMPTY_GRAPH, UNKNOWN_OBJECTIVE,
    INVALID_CONFIG, INVALID_PARTITION, and for ``erbga evolve`` LLM_UNAVAILABLE
    and LLM_ERROR.
    """

    def __init__(self, code: str, message: str, hint: str = ""):
        super().__init__(message)
        self.code, self.message, self.hint = code, message, hint

    def to_dict(self) -> dict:
        return {"ok": False, "error": {"code": self.code, "message": self.message,
                                       "hint": self.hint}}


# ---- graph input ------------------------------------------------------

def load_graph(spec: Any) -> Graph:
    """Turn a graph spec (see module docstring) into a ``Graph``."""
    if isinstance(spec, Graph):
        graph = spec
    elif isinstance(spec, dict):
        graph = _graph_from_dict(spec, name="inline")
    elif isinstance(spec, str) and spec.lstrip().startswith("{"):
        try:
            graph = _graph_from_dict(json.loads(spec), name="inline")
        except json.JSONDecodeError as exc:
            raise ERBGAError("INVALID_GRAPH",
                             f"graph looks like JSON but is invalid: {exc}") from None
    elif isinstance(spec, (str, Path)):
        graph = _graph_from_path_or_name(str(spec))
    else:
        raise ERBGAError("INVALID_GRAPH", f"unsupported graph spec of type {type(spec).__name__}",
                         "pass a dataset name, a file path, or {'edges': [[u, v], ...]}")
    if graph.m == 0:
        raise ERBGAError("EMPTY_GRAPH", f"graph {graph.name!r} has no edges",
                         "community detection needs at least one edge")
    return graph


def _graph_from_path_or_name(spec: str) -> Graph:
    if spec in BUNDLED_DATASETS:
        ref = resources.files("erbga") / "data" / f"{spec}.gml"
        with resources.as_file(ref) as path:
            return read_gml(path)
    path = Path(spec).expanduser()
    if not path.is_file():
        raise ERBGAError("GRAPH_NOT_FOUND", f"no such file or bundled dataset: {spec!r}",
                         f"bundled datasets: {', '.join(BUNDLED_DATASETS)}")
    suffix = path.suffix.lower()
    try:
        if suffix == ".gml":
            return read_gml(path)
        if suffix == ".json":
            return _graph_from_dict(json.loads(path.read_text()), name=path.stem)
        return read_edge_list(path)
    except ERBGAError:
        raise
    except Exception as exc:
        raise ERBGAError("INVALID_GRAPH", f"could not read {spec!r}: {exc}",
                         "supported: .gml, .json, or an edge list with one 'u v' per line") from exc


def _graph_from_dict(data: dict, name: str) -> Graph:
    edges = data.get("edges", data.get("links"))
    if edges is None:
        raise ERBGAError("INVALID_GRAPH", "graph dict needs an 'edges' (or 'links') list",
                         "example: {'edges': [[1, 2], [2, 3]]}")
    pairs = []
    for i, e in enumerate(edges):
        if isinstance(e, dict) and "source" in e and "target" in e:
            pairs.append((e["source"], e["target"]))
        elif isinstance(e, (list, tuple)) and len(e) >= 2:
            pairs.append((e[0], e[1]))
        else:
            raise ERBGAError("INVALID_GRAPH",
                             f"edge #{i} is not [u, v] or {{source, target}}: {e!r}")
    nodes = data.get("nodes")
    if nodes is not None:
        nodes = [n["id"] if isinstance(n, dict) else n for n in nodes]
    return Graph.from_edges(pairs, nodes=nodes, name=str(data.get("name", name)))


def read_edge_list(path: str | Path) -> Graph:
    """Read 'u v' (or 'u,v' / tab-separated) pairs; extra columns are ignored."""
    pairs = []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line[0] in "#%":
            continue
        parts = [p for p in re.split(r"[,\s;]+", line) if p]
        if len(parts) < 2:
            raise ERBGAError("INVALID_GRAPH", f"edge list line has fewer than two fields: {line!r}")
        pairs.append((_parse_id(parts[0]), _parse_id(parts[1])))
    return Graph.from_edges(pairs, name=Path(path).stem)


def _parse_id(token: str):
    try:
        return int(token)
    except ValueError:
        return token


def _json_id(x):
    x = x.item() if hasattr(x, "item") else x
    return x if isinstance(x, (int, str)) else str(x)


# ---- tools ------------------------------------------------------------

def list_objectives() -> dict:
    """Names and one-line descriptions of the built-in objectives."""
    return {"ok": True, "objectives": [
        {"name": name, "description": _first_paragraph(fn.__doc__)}
        for name, fn in OBJECTIVES.items()
    ]}


def graph_summary(graph: Any) -> dict:
    """Size, density, degree statistics and connectivity of a graph."""
    g = load_graph(graph)
    comps = np.bincount(decode(g, np.zeros(g.m, dtype=bool)))
    deg = g.degree
    return {
        "ok": True,
        "graph": g.name,
        "nodes": g.n,
        "edges": g.m,
        "density": round(2 * g.m / (g.n * (g.n - 1)), 6) if g.n > 1 else 0.0,
        "degree": {"min": int(deg.min()), "max": int(deg.max()),
                   "mean": round(float(deg.mean()), 3)},
        "connected_components": int(len(comps)),
        "largest_component": int(comps.max()),
        "isolated_nodes": int(np.sum(deg == 0)),
    }


def detect_communities(graph: Any, objective: str = "modularity", *,
                       max_nodes_per_community: int | None = None,
                       include_history: bool = False, progress=None, operators=(),
                       **config) -> dict:
    """Run the GA and return communities plus run metadata.

    ``config`` keys are the ``GAConfig`` fields: pop_size, generations,
    elite_frac, tournament_size, crossover_rate, crossover, mutation_rate,
    init, local_search, patience, seed. ``max_nodes_per_community`` truncates
    the node lists (the size is always exact) to keep output small.
    ``progress(generation, (best, mean))`` is called after each generation.
    ``operators`` adds mutation operators, see ``ERBGA``.
    """
    g = load_graph(graph)
    fn = _objective(objective)
    unknown = sorted(set(config) - set(CONFIG_FIELDS))
    if unknown:
        raise ERBGAError("INVALID_CONFIG", f"unknown config option(s): {', '.join(unknown)}",
                         f"valid options: {', '.join(CONFIG_FIELDS)}")
    try:
        cfg = GAConfig(**config)
    except (ValueError, TypeError) as exc:
        raise ERBGAError("INVALID_CONFIG", str(exc),
                         "see GAConfig in erbga/ga.py for ranges") from exc

    name = objective if isinstance(objective, str) else getattr(fn, "__name__", "custom")
    result = ERBGA(g, fn, cfg, operators).run(progress)
    best = [b for b, _ in result.history]
    stopped_early = result.generations < cfg.generations
    communities = _community_list(g, result.labels, max_nodes_per_community)
    out = {
        "ok": True,
        "summary": (f"{len(communities)} communities in {g.name} ({g.n} nodes, {g.m} edges); "
                    f"{name} = {result.fitness:.6f}; largest community has "
                    f"{communities[0]['size']} nodes"),
        "graph": {"name": g.name, "nodes": g.n, "edges": g.m},
        "objective": name,
        "fitness": result.fitness,
        "num_communities": len(communities),
        "communities": communities,
        "run": {
            "generations": result.generations,
            "stopped_early": stopped_early,
            "evaluations": result.evaluations,
            "seconds": round(result.seconds, 3),
            "initial_best": best[0],
            "generation_of_best": int(np.argmax(np.array(best) >= result.fitness - 1e-12)),
        },
        "config": asdict(cfg),
    }
    if include_history:
        out["history"] = {"best": best, "mean": [m for _, m in result.history]}
    return out


def score_partition(graph: Any, partition: Any, objectives: str | list[str] | None = None) -> dict:
    """Score a given partition without running the GA.

    ``partition`` is ``{node_id: community}`` or a list of node-id lists.
    Node ids are matched as strings, so JSON object keys work. Every node must
    be assigned exactly once. ``objectives`` defaults to all built-ins.
    """
    g = load_graph(graph)
    labels = _labels_from_partition(g, partition)
    names = list(OBJECTIVES) if objectives is None else (
        [objectives] if isinstance(objectives, str) else list(objectives))
    scores = {name: float(_objective(name)(g, labels)) for name in names}

    pieces = decode(g, canonical_cut(g, labels))
    disconnected = sorted({int(c) for c in np.unique(labels)
                           if len(np.unique(pieces[labels == c])) > 1})
    return {
        "ok": True,
        "graph": {"name": g.name, "nodes": g.n, "edges": g.m},
        "num_communities": int(len(np.unique(labels))),
        "scores": scores,
        "connected": not disconnected,
        "disconnected_communities": disconnected,
    }


# ---- helpers ----------------------------------------------------------

def _objective(name):
    try:
        return get_objective(name)
    except ValueError as exc:
        raise ERBGAError("UNKNOWN_OBJECTIVE", str(exc), "call list_objectives() for names") from exc


def _community_list(g: Graph, labels: np.ndarray, limit: int | None) -> list[dict]:
    out = []
    for c in range(int(labels.max()) + 1):
        members = [_json_id(x) for x in g.node_ids[labels == c]]
        entry = {"size": len(members), "nodes": members}
        if limit is not None and len(members) > limit:
            entry["nodes"], entry["truncated"] = members[:limit], True
        out.append(entry)
    out.sort(key=lambda e: -e["size"])
    for i, entry in enumerate(out):
        entry["id"] = i
    return [{"id": e.pop("id"), **e} for e in out]


def _labels_from_partition(g: Graph, partition: Any) -> np.ndarray:
    index = {str(_json_id(nid)): i for i, nid in enumerate(g.node_ids)}
    labels = np.full(g.n, -1, dtype=np.int64)
    if isinstance(partition, dict):
        items = [(node, comm) for node, comm in partition.items()]
    elif isinstance(partition, (list, tuple)):
        items = [(node, c) for c, group in enumerate(partition) for node in group]
    else:
        raise ERBGAError("INVALID_PARTITION", "partition must be a dict or a list of lists",
                         "examples: {'1': 0, '2': 0, '3': 1} or [[1, 2], [3]]")
    comm_ids: dict[str, int] = {}
    unknown = []
    for node, comm in items:
        i = index.get(str(node))
        if i is None:
            unknown.append(node)
            continue
        if labels[i] != -1:
            raise ERBGAError("INVALID_PARTITION", f"node {node!r} is assigned more than once")
        labels[i] = comm_ids.setdefault(str(comm), len(comm_ids))
    if unknown:
        raise ERBGAError("INVALID_PARTITION",
                         f"{len(unknown)} node id(s) not in the graph, e.g. {unknown[:5]}",
                         "node ids must match the graph's ids (compared as strings)")
    missing = [_json_id(g.node_ids[i]) for i in np.flatnonzero(labels == -1)]
    if missing:
        raise ERBGAError("INVALID_PARTITION",
                         f"{len(missing)} node(s) have no community, e.g. {missing[:5]}",
                         "assign every node; use its own community for a singleton")
    return labels


def _first_paragraph(doc: str | None) -> str:
    if not doc:
        return ""
    return " ".join(doc.strip().split("\n\n")[0].split())
