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


# ---- incremental scoring of moves ----------------------------------------

# A move gain scores moving one unit (a node, or at coarser levels of
# multilevel local search a whole group of nodes) from its community into each
# of several neighbouring communities, without recomputing the objective.
#
#     gain(m, n, unit_degree, unit_size, links, own_links,
#          target_degree, target_size, own_degree, own_size) -> array
#
# m, n: edges and nodes of the original graph. unit_degree / unit_size: the
# unit's summed degree and node count. links[i]: edges between the unit and
# target community i; own_links: edges between the unit and the rest of its own
# community. target_*: degree sum and node count of each target; own_*: the
# same for the unit's community *without* the unit.
#
# Gains are exact for the labelling as given. If a move disconnects a
# community, decoding splits it, and for these objectives splitting a
# community with no edges between the pieces only raises the score, so the
# gain is a lower bound on the decoded result.

def modularity_gain(m, n, unit_degree, unit_size, links, own_links,
                    target_degree, target_size, own_degree, own_size):
    return ((links - own_links) / m
            - unit_degree * (target_degree - own_degree) / (2.0 * m * m))


def qs_gain(m, n, unit_degree, unit_size, links, own_links,
            target_degree, target_size, own_degree, own_size):
    def pairs(x):
        return x * (x - 1.0)

    size_change = (pairs(target_size + unit_size) - pairs(target_size)
                   + pairs(own_size) - pairs(own_size + unit_size))
    return (links - own_links) / m - size_change / (n * (n - 1.0))


MOVE_GAINS = {modularity: modularity_gain, qs: qs_gain}


def get_objective(objective: str | Objective) -> Objective:
    if callable(objective):
        return objective
    try:
        return OBJECTIVES[objective]
    except KeyError:
        raise ValueError(
            f"unknown objective {objective!r}; choose from {sorted(OBJECTIVES)}") from None
