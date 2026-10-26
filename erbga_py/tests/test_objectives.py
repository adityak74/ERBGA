import numpy as np
import pytest

from erbga import Graph, decode, modularity, qs

nx = pytest.importorskip("networkx")


def to_nx(g):
    G = nx.Graph()
    G.add_nodes_from(range(g.n))
    G.add_edges_from(zip(g.src.tolist(), g.dst.tolist(), strict=True))
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


@pytest.mark.parametrize("name", ["modularity", "qs"])
@pytest.mark.parametrize("coarse", [False, True])
def test_move_gain_matches_full_recompute(football, name, coarse):
    """Every gain equals the full objective change, for single nodes (level 0)
    and for whole groups of nodes on an aggregated level."""
    from erbga.multilevel import Level
    from erbga.objectives import MOVE_GAINS, OBJECTIVES

    objective = OBJECTIVES[name]
    gain = MOVE_GAINS[objective]
    rng = np.random.default_rng(5)
    level, membership = Level.from_graph(football), np.arange(football.n)
    if coarse:  # units are 40 random groups of nodes
        level, membership = level.aggregate(rng.integers(0, 40, football.n))
    part = rng.integers(0, 8, level.n)  # communities of units, not necessarily connected
    for _ in range(150):
        v = int(rng.integers(level.n))
        lo, hi = level.ptr[v], level.ptr[v + 1]
        targets, inverse = np.unique(part[level.nbr[lo:hi]], return_inverse=True)
        links = np.bincount(inverse, weights=level.weight[lo:hi])
        own = targets == part[v]
        if own.all():
            continue
        k = max(part.max(), targets.max()) + 1
        degree = np.bincount(part, weights=level.degree, minlength=k)
        size = np.bincount(part, weights=level.size, minlength=k)
        predicted = gain(football.m, football.n, level.degree[v], level.size[v], links[~own],
                         links[own].sum(), degree[targets[~own]], size[targets[~own]],
                         degree[part[v]] - level.degree[v], size[part[v]] - level.size[v])
        before = objective(football, part[membership])
        for t, p in zip(targets[~own], predicted, strict=True):
            trial = part.copy()
            trial[v] = t
            assert p == pytest.approx(objective(football, trial[membership]) - before, abs=1e-12)
        part[v] = targets[~own][rng.integers((~own).sum())]


def test_aggregate_preserves_totals(football):
    from erbga.multilevel import Level

    base = Level.from_graph(football)
    part = np.random.default_rng(0).integers(0, 12, football.n)
    level, groups = base.aggregate(part)
    assert level.n == len(np.unique(part))
    assert level.degree.sum() == pytest.approx(2 * football.m)
    assert level.size.sum() == football.n
    inside = (part[football.src] == part[football.dst]).sum()
    assert level.weight.sum() == pytest.approx(2 * (football.m - inside))  # both directions
