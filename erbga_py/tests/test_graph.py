import numpy as np

from erbga import Graph


def test_reads_datasets(karate, football):
    assert (karate.n, karate.m) == (34, 78)
    assert (football.n, football.m) == (115, 613)
    assert karate.degree.sum() == 2 * karate.m
    assert np.all(karate.src < karate.dst)


def test_node_ids_preserved(karate):
    assert sorted(karate.node_ids.tolist()) == list(range(1, 35))


def test_from_edges_drops_loops_and_duplicates():
    g = Graph.from_edges([("a", "b"), ("b", "a"), ("c", "c"), ("b", "c")])
    assert g.n == 3
    assert g.m == 2


def test_adjacency_matches_edges(karate):
    for v in range(karate.n):
        for e, w in zip(karate.incident_edges(v), karate.neighbours(v)):
            assert {karate.src[e], karate.dst[e]} == {v, w}
