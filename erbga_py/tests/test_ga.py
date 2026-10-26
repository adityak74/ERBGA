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
    assert all(later >= earlier - 1e-12 for earlier, later in zip(best, best[1:], strict=False))


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


@pytest.mark.parametrize("objective", ["modularity", "qs"])
def test_fast_refine_never_worse_and_connected(football, objective):
    ga = ERBGA(football, objective, GAConfig(seed=0))
    rng = np.random.default_rng(1)
    for _ in range(5):
        start, start_fit, _ = ga.evaluate(rng.random(football.m) < 0.5)
        cut, fit = ga.refine(start)
        assert fit >= start_fit - 1e-12
        _, exact, labels = ga.evaluate(cut)
        assert fit == pytest.approx(exact)
        assert is_connected_partition(football, labels)


def test_custom_objective_uses_exact_refine(karate):
    from erbga.objectives import modularity as mod

    def wrapped(graph, labels):  # not in MOVE_GAINS, so the slow path runs
        return mod(graph, labels)

    ga = ERBGA(karate, wrapped, GAConfig(seed=0))
    cut, fit = ga.refine(np.ones(karate.m, dtype=bool))  # start from singletons
    assert fit > 0.3


@pytest.mark.parametrize("bad", [
    {"islands": 0}, {"migration_interval": 0}, {"migrants": 0},
    {"pop_size": 10, "islands": 5, "migrants": 2},  # islands of 2 can't send 2 migrants
])
def test_invalid_island_config_rejected(bad):
    with pytest.raises(ValueError):
        GAConfig(**bad)


def test_islands_split_population_and_stay_monotone(karate):
    cfg = GAConfig(seed=1, pop_size=23, islands=4, migration_interval=3,
                   generations=12, patience=None)
    result = ERBGA(karate, config=cfg).run()
    best = [b for b, _ in result.history]
    assert all(later >= earlier - 1e-12 for earlier, later in zip(best, best[1:], strict=False))
    assert is_connected_partition(karate, result.labels)


def test_islands_reproducible(karate):
    a = detect_communities(karate, seed=4, islands=3, generations=15)
    b = detect_communities(karate, seed=4, islands=3, generations=15)
    assert a.history == b.history and np.array_equal(a.labels, b.labels)


def test_migration_replaces_worst_with_neighbours_best(karate):
    ga = ERBGA(karate, config=GAConfig(seed=0, islands=3, migrants=2))
    rng = np.random.default_rng(0)
    islands = []
    for i in range(3):
        pop = rng.random((5, karate.m)) < 0.5
        islands.append((pop, np.array([0.1, 0.5, 0.3, 0.9, 0.2]) + i))
    before = [(p.copy(), f.copy()) for p, f in islands]
    ga._migrate(islands)
    for i in range(3):
        src_pop, src_fit = before[i - 1]
        top2 = np.argsort(-src_fit)[:2]
        worst2 = np.argsort(before[i][1])[:2]  # indices 0 and 4
        assert np.array_equal(islands[i][1][worst2], src_fit[top2])
        assert np.array_equal(islands[i][0][worst2], src_pop[top2])
        untouched = [j for j in range(5) if j not in worst2]
        assert np.array_equal(islands[i][1][untouched], before[i][1][untouched])


def test_aggregation_escapes_singleton_trap():
    """From all-singletons, node moves alone stall on a fine partition; moving
    whole communities on the aggregated graph reaches the planted structure."""
    rng = np.random.default_rng(0)
    n, blocks = 600, 10
    block = np.repeat(np.arange(blocks), n // blocks)
    iu, ju = np.triu_indices(n, 1)
    keep = rng.random(len(iu)) < np.where(block[iu] == block[ju], 50 / n, 2 / n)
    g = Graph.from_edges(zip(iu[keep].tolist(), ju[keep].tolist(), strict=True), nodes=range(n))
    planted = modularity(g, block)

    start = np.ones(g.m, dtype=bool)
    _, node_only = ERBGA(g, config=GAConfig(seed=0, aggregate=False)).refine(start)
    cut, multilevel = ERBGA(g, config=GAConfig(seed=0)).refine(start)
    assert node_only < planted - 0.2  # the trap
    assert multilevel > planted - 0.01  # escaped it
    assert is_connected_partition(g, ERBGA(g).evaluate(cut)[2])


def test_cli_no_aggregate_flag(capsys):
    from conftest import DATA
    assert main([str(DATA / "karate.gml"), "--seed", "0", "--no-aggregate",
                 "--format", "json"]) == 0
    assert json.loads(capsys.readouterr().out)["config"]["aggregate"] is False


@pytest.mark.parametrize("theta", [0.0, 0.01, 1.0])
@pytest.mark.parametrize("objective", ["modularity", "qs"])
def test_refined_subcommunities_are_nested_and_connected(football, objective, theta):
    """Leiden refinement only merges a unit into an adjacent sub-community of
    its own community, so every refined piece is connected and nested."""
    from erbga.multilevel import Level, local_moving, refine_partition
    from erbga.objectives import MOVE_GAINS, OBJECTIVES

    gain = MOVE_GAINS[OBJECTIVES[objective]]
    rng = np.random.default_rng(3)
    level = Level.from_graph(football)
    part, _ = local_moving(level, np.arange(football.n), gain, football.m, football.n, rng)
    refined = refine_partition(level, part, gain, football.m, football.n, rng, theta)
    for r in np.unique(refined):
        members = np.flatnonzero(refined == r)
        assert len(np.unique(part[members])) == 1  # nested in one community
    assert is_connected_partition(football, np.unique(refined, return_inverse=True)[1])
    assert len(np.unique(refined)) < football.n  # it did merge something


def test_refinement_flag(karate):
    on = ERBGA(karate, config=GAConfig(seed=0, refinement=True)).refine(
        np.ones(karate.m, dtype=bool))
    off = ERBGA(karate, config=GAConfig(seed=0, refinement=False)).refine(
        np.ones(karate.m, dtype=bool))
    assert on[1] > 0.4 and off[1] > 0.4


def test_cli_refinement_flag(capsys):
    from conftest import DATA
    assert main([str(DATA / "karate.gml"), "--seed", "0", "--refinement", "--format", "json"]) == 0
    assert json.loads(capsys.readouterr().out)["config"]["refinement"] is True


def test_randomized_refinement_reproducible_and_validated(karate):
    cfg = dict(seed=5, refinement=True, refinement_theta=0.5, generations=10)
    a, b = detect_communities(karate, **cfg), detect_communities(karate, **cfg)
    assert a.history == b.history and np.array_equal(a.labels, b.labels)
    with pytest.raises(ValueError):
        GAConfig(refinement_theta=-0.1)
