import json

import numpy as np
import pytest

from erbga import ERBGA, GAConfig, Graph, detect_communities, modularity
from erbga.cli import main
from erbga.partition import is_connected_partition

KARATE_BEST_MODULARITY = 0.4198


def test_reaches_known_optimum_on_karate(karate):
    result = detect_communities(karate, "modularity", seed=0)
    assert result.fitness == pytest.approx(KARATE_BEST_MODULARITY, abs=1e-4)
    assert result.num_communities == 4
    assert modularity(karate, result.labels) == pytest.approx(result.fitness)


def test_communities_are_connected(football):
    result = detect_communities(football, "modularity", seed=0, generations=30)
    assert is_connected_partition(football, result.labels)


def test_same_seed_same_result(karate):
    a = detect_communities(karate, seed=7, generations=20)
    b = detect_communities(karate, seed=7, generations=20)
    assert a.fitness == b.fitness
    assert np.array_equal(a.labels, b.labels)
    assert a.history == b.history


@pytest.mark.parametrize("objective", ["modularity", "qs"])
@pytest.mark.parametrize("local_search", [True, False])
def test_best_never_gets_worse(karate, objective, local_search):
    cfg = GAConfig(seed=3, generations=40, local_search=local_search, patience=None)
    result = ERBGA(karate, objective, cfg).run()
    best = [b for b, _ in result.history]
    assert all(later >= earlier - 1e-12 for earlier, later in zip(best, best[1:]))


@pytest.mark.parametrize("init", ["random", "degree"])
@pytest.mark.parametrize("crossover", ["uniform", "one_point"])
def test_config_options_run(karate, init, crossover):
    result = detect_communities(karate, seed=0, generations=5, init=init, crossover=crossover)
    assert result.fitness > 0


@pytest.mark.parametrize("bad", [
    {"pop_size": 0}, {"tournament_size": 0}, {"generations": -1}, {"elite_frac": 1.5},
    {"crossover": "two_point"}, {"init": "smart"}, {"patience": 0},
])
def test_invalid_config_rejected(bad):
    with pytest.raises(ValueError):
        GAConfig(**bad)


@pytest.mark.parametrize("crossover", ["uniform", "one_point"])
def test_single_edge_graph(crossover):
    g = Graph.from_edges([(0, 1)])
    result = detect_communities(g, seed=0, generations=3, pop_size=4, crossover=crossover)
    assert result.num_communities == 1


def test_custom_objective(karate):
    calls = []

    def fewest_communities(graph, labels):
        calls.append(1)
        return -float(labels.max())

    result = detect_communities(karate, fewest_communities, seed=0, generations=10)
    assert calls
    assert result.num_communities == 1


def test_cli(tmp_path, capsys):
    from conftest import DATA
    out = tmp_path / "r.json"
    assert main([str(DATA / "karate.gml"), "--seed", "0", "--json", str(out)]) == 0
    assert "best modularity = 0.41979" in capsys.readouterr().out
    assert json.loads(out.read_text())["num_communities"] == 4
