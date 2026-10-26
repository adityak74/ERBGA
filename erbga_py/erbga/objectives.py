"""Community-quality objectives. Higher is better for all of them.

An objective is any callable ``f(graph, labels) -> float`` where ``labels[i]``
is the community (0..k-1) of node i. Register new ones in ``OBJECTIVES`` or
pass the callable straight to the GA.
"""

from __future__ import annotations

from collections.abc import Callable

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


# ---- incremental scoring of single-node moves ---------------------------

class CommunityStats:
    """Labels plus per-community degree sums and sizes, kept current under
    single-node moves, so a move can be scored in O(degree) instead of O(m)."""

    def __init__(self, graph: Graph, labels: np.ndarray):
        self.labels = labels.copy()
        k = int(labels.max()) + 1
        self.degree = np.bincount(labels, weights=graph.degree, minlength=k)
        self.size = np.bincount(labels, minlength=k)

    def links(self, graph: Graph, v: int) -> tuple[np.ndarray, np.ndarray]:
        """Communities adjacent to v and how many of v's edges go to each."""
        return np.unique(self.labels[graph.neighbours(v)], return_counts=True)

    def move(self, graph: Graph, v: int, target: int) -> None:
        source, kv = self.labels[v], graph.degree[v]
        self.degree[source] -= kv
        self.degree[target] += kv
        self.size[source] -= 1
        self.size[target] += 1
        self.labels[v] = target


# A move gain gives, for node v and each community in ``targets``, the
# objective change from moving v there: ``links[i]`` is the number of v's edges
# into ``targets[i]`` and ``own`` the number into v's current community.
# Signature: gain(graph, stats, v, targets, links, own) -> array.
# Both gains below treat the partition as given. If the move disconnects v's
# old community, decoding splits it, and for these objectives splitting a
# community with no edges between the pieces only raises the score, so the
# gain is a lower bound on the decoded result.

def modularity_gain(graph: Graph, stats: CommunityStats, v: int, targets: np.ndarray,
                    links: np.ndarray, own: int) -> np.ndarray:
    m, kv = graph.m, graph.degree[v]
    d_own = stats.degree[stats.labels[v]] - kv
    return (links - own) / m - kv * (stats.degree[targets] - d_own) / (2.0 * m * m)


def qs_gain(graph: Graph, stats: CommunityStats, v: int, targets: np.ndarray,
            links: np.ndarray, own: int) -> np.ndarray:
    m, n = graph.m, graph.n
    n_own = stats.size[stats.labels[v]] - 1
    return (links - own) / m - 2.0 * (stats.size[targets] - n_own) / (n * (n - 1.0))


MOVE_GAINS = {modularity: modularity_gain, qs: qs_gain}


def get_objective(objective: str | Objective) -> Objective:
    if callable(objective):
        return objective
    try:
        return OBJECTIVES[objective]
    except KeyError:
        raise ValueError(
            f"unknown objective {objective!r}; choose from {sorted(OBJECTIVES)}") from None
