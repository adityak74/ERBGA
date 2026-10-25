import numpy as np
import pytest

from erbga import Graph, decode, modularity, qs

nx = pytest.importorskip("networkx")


def to_nx(g):
    G = nx.Graph()
    G.add_nodes_from(range(g.n))
    G.add_edges_from(zip(g.src.tolist(), g.dst.tolist()))
    return G


def as_sets(labels):
    return [set(np.flatnonzero(labels == c).tolist()) for c in range(labels.max() + 1)]


def test_modularity_matches_networkx(karate, football):
    rng = np.random.default_rng(1)
    for g in (karate, football):
        G = to_nx(g)
        for p in (0.1, 0.3, 0.6):
            labels = decode(g, rng.random(g.m) < p)
            expected = nx.community.modularity(G, as_sets(labels))
            assert modularity(g, labels) == pytest.approx(expected, abs=1e-12)


def test_modularity_known_values(karate):
    assert modularity(karate, np.zeros(karate.n, dtype=int)) == pytest.approx(0.0)
    # Zachary's club split (Mr. Hi vs Officer); unweighted Q is 0.3582
    G = nx.karate_club_graph()
    club = np.array([G.nodes[i]["club"] == "Officer" for i in range(34)], dtype=int)
    order = np.argsort(karate.node_ids)  # GML ids 1..34 -> networkx 0..33
    labels = np.empty(karate.n, dtype=int)
    labels[order] = club
    assert modularity(karate, labels) == pytest.approx(0.3582, abs=1e-4)


def test_qs_by_hand():
    # two triangles joined by one edge: N=6, m=7
    g = Graph.from_edges([(0, 1), (1, 2), (0, 2), (3, 4), (4, 5), (3, 5), (2, 3)])
    labels = np.array([0, 0, 0, 1, 1, 1])
    expected = 2 * (3 / 7 - 3 * 2 / (6 * 5))
    assert qs(g, labels) == pytest.approx(expected)
