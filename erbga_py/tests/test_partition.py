import numpy as np

from erbga import Graph, canonical_cut, connected_components, decode
from erbga.partition import is_connected_partition


def path_and_triangle():
    # 0-1-2 path, 3-4-5 triangle, 6 isolated
    return Graph.from_edges([(0, 1), (1, 2), (3, 4), (4, 5), (3, 5)], nodes=range(7))


def test_connected_components():
    g = path_and_triangle()
    labels = connected_components(g.n, g.src, g.dst)
    assert labels.tolist() == [0, 0, 0, 1, 1, 1, 2]


def test_long_path_converges():
    n = 1000
    labels = connected_components(n, np.arange(n - 1)[::-1].copy(), np.arange(1, n)[::-1].copy())
    assert np.all(labels == 0)


def test_cut_everything_gives_singletons(karate):
    labels = decode(karate, np.ones(karate.m, dtype=bool))
    assert labels.max() + 1 == karate.n


def test_cut_nothing_gives_one_community(karate):
    labels = decode(karate, np.zeros(karate.m, dtype=bool))
    assert np.all(labels == 0)


def test_canonical_cut_is_fixed_point(karate):
    rng = np.random.default_rng(0)
    for _ in range(20):
        labels = decode(karate, rng.random(karate.m) < 0.4)
        cut = canonical_cut(karate, labels)
        assert np.array_equal(decode(karate, cut), labels)
        assert is_connected_partition(karate, labels)


def test_detects_disconnected_partition():
    g = path_and_triangle()
    labels = np.array([0, 1, 0, 2, 2, 2, 3])  # nodes 0 and 2 only meet through 1
    assert not is_connected_partition(g, labels)
