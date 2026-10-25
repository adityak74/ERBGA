"""Community-quality objectives. Higher is better for all of them.

An objective is any callable ``f(graph, labels) -> float`` where ``labels[i]``
is the community (0..k-1) of node i. Register new ones in ``OBJECTIVES`` or
pass the callable straight to the GA.
"""

from __future__ import annotations

from typing import Callable

import numpy as np

from .graph import Graph

Objective = Callable[[Graph, np.ndarray], float]


def _intra_edges(graph: Graph, labels: np.ndarray, k: int) -> np.ndarray:
    lu = labels[graph.src]
    inside = lu == labels[graph.dst]
    return np.bincount(lu[inside], minlength=k)


def modularity(graph: Graph, labels: np.ndarray) -> float:
    """Newman-Girvan modularity: sum_c [ L_c/m - (d_c / 2m)^2 ]."""
    m = graph.m
    if m == 0:
        return 0.0
    k = int(labels.max()) + 1
    intra = _intra_edges(graph, labels, k)
    deg = np.bincount(labels, weights=graph.degree, minlength=k)
    return float(np.sum(intra / m - (deg / (2.0 * m)) ** 2))


def qs(graph: Graph, labels: np.ndarray) -> float:
    """Density-based score: sum_c [ L_c/m - n_c(n_c-1) / (N(N-1)) ].

    The expected share of edges in a community is its share of all node pairs,
    not its share of degree, as in the 2018 C++ ``Network::q_calc``. The C++
    version skipped two-node communities; this one counts every community.
    """
    m, n = graph.m, graph.n
    if m == 0 or n < 2:
        return 0.0
    k = int(labels.max()) + 1
    intra = _intra_edges(graph, labels, k)
    size = np.bincount(labels, minlength=k).astype(float)
    return float(np.sum(intra / m - size * (size - 1) / (n * (n - 1.0))))


OBJECTIVES: dict[str, Objective] = {
    "modularity": modularity,
    "qs": qs,
}


def get_objective(objective: str | Objective) -> Objective:
    if callable(objective):
        return objective
    try:
        return OBJECTIVES[objective]
    except KeyError:
        raise ValueError(
            f"unknown objective {objective!r}; choose from {sorted(OBJECTIVES)}") from None
