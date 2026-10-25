"""Undirected graph storage and a small GML reader."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

_TOKEN = re.compile(r'"[^"]*"|\[|\]|[^\s\[\]]+')


@dataclass
class Graph:
    """Simple undirected, unweighted graph with contiguous node indices 0..n-1.

    ``node_ids[i]`` is the id the node had in the input file. ``src``/``dst``
    hold one entry per undirected edge with ``src < dst``.
    """

    node_ids: np.ndarray
    src: np.ndarray
    dst: np.ndarray
    name: str = ""
    degree: np.ndarray = field(init=False)
    incident_ptr: np.ndarray = field(init=False)
    incident_edge: np.ndarray = field(init=False)
    incident_nbr: np.ndarray = field(init=False)

    def __post_init__(self) -> None:
        n = self.n
        self.degree = np.bincount(self.src, minlength=n) + np.bincount(self.dst, minlength=n)

        # CSR adjacency: for node v, incident_edge[incident_ptr[v]:incident_ptr[v+1]]
        # are the edge indices touching v and incident_nbr the matching neighbours.
        ends = np.concatenate([self.src, self.dst])
        others = np.concatenate([self.dst, self.src])
        edge_idx = np.concatenate([np.arange(self.m), np.arange(self.m)])
        order = np.argsort(ends, kind="stable")
        self.incident_edge = edge_idx[order]
        self.incident_nbr = others[order]
        self.incident_ptr = np.concatenate([[0], np.cumsum(self.degree)])

    @property
    def n(self) -> int:
        return len(self.node_ids)

    @property
    def m(self) -> int:
        return len(self.src)

    @classmethod
    def from_edges(cls, edges, nodes=None, name: str = "") -> Graph:
        """Build a graph from (u, v) pairs of arbitrary hashable node ids.

        Self-loops and duplicate edges are dropped.
        """
        node_list = list(nodes) if nodes is not None else []
        index = {v: i for i, v in enumerate(node_list)}
        pairs = set()
        for u, v in edges:
            for x in (u, v):
                if x not in index:
                    index[x] = len(node_list)
                    node_list.append(x)
            a, b = index[u], index[v]
            if a != b:
                pairs.add((min(a, b), max(a, b)))
        arr = np.array(sorted(pairs), dtype=np.int64).reshape(-1, 2)
        return cls(np.array(node_list), arr[:, 0].copy(), arr[:, 1].copy(), name)

    def neighbours(self, v: int) -> np.ndarray:
        return self.incident_nbr[self.incident_ptr[v]:self.incident_ptr[v + 1]]

    def incident_edges(self, v: int) -> np.ndarray:
        return self.incident_edge[self.incident_ptr[v]:self.incident_ptr[v + 1]]


def read_gml(path: str | Path) -> Graph:
    """Read the node ids and edges of a GML file (attributes are ignored)."""
    text = Path(path).read_text()
    tokens = _TOKEN.findall(text)
    nodes, edges = [], []
    stack: list[dict] = []
    key = None
    for tok in tokens:
        if tok == "[":
            stack.append({"_kind": key})
            key = None
        elif tok == "]":
            block = stack.pop()
            if block["_kind"] == "node":
                nodes.append(block["id"])
            elif block["_kind"] == "edge":
                edges.append((block["source"], block["target"]))
        elif key is None:
            key = tok
        else:
            if stack:
                stack[-1][key] = _parse_value(tok)
            key = None
    if not stack and not nodes and not edges:
        raise ValueError(f"no graph found in {path}")
    return Graph.from_edges(edges, nodes=nodes, name=Path(path).stem)


def _parse_value(tok: str):
    if tok.startswith('"'):
        return tok[1:-1]
    try:
        return int(tok)
    except ValueError:
        try:
            return float(tok)
        except ValueError:
            return tok
