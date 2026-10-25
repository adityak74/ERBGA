"""Decoding an edge-cut chromosome into communities."""

from __future__ import annotations

import numpy as np

from .graph import Graph


def connected_components(n: int, src: np.ndarray, dst: np.ndarray) -> np.ndarray:
    """Label the connected components of the graph on n nodes with the given edges.

    Returns labels in 0..k-1, numbered in order of each component's lowest node.
    Uses min-label hooking plus pointer jumping, so every step is vectorised.
    """
    labels = np.arange(n)
    if len(src) == 0:
        return labels
    while True:
        low = np.minimum(labels[src], labels[dst])
        np.minimum.at(labels, src, low)
        np.minimum.at(labels, dst, low)
        while True:  # pointer jumping: point every node at its root
            jumped = labels[labels]
            if np.array_equal(jumped, labels):
                break
            labels = jumped
        if np.array_equal(labels[src], labels[dst]):
            break
    return np.unique(labels, return_inverse=True)[1]


def decode(graph: Graph, cut: np.ndarray) -> np.ndarray:
    """Communities are the connected components left after removing cut edges."""
    keep = ~cut
    return connected_components(graph.n, graph.src[keep], graph.dst[keep])


def canonical_cut(graph: Graph, labels: np.ndarray) -> np.ndarray:
    """The minimal cut that yields ``labels``: exactly the inter-community edges.

    Cut edges that lie inside a community change nothing in the decoded
    partition. Restoring them makes genotype and phenotype one-to-one.
    """
    return labels[graph.src] != labels[graph.dst]


def is_connected_partition(graph: Graph, labels: np.ndarray) -> bool:
    """True if every community induces a connected subgraph."""
    return np.array_equal(decode(graph, canonical_cut(graph, labels)),
                          np.unique(labels, return_inverse=True)[1])
