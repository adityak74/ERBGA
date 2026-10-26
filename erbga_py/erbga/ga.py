"""The genetic algorithm.

Chromosome: one bit per edge, True = edge is cut. Communities are the connected
components of the graph that remains, so every community is connected.

After each child is built it is repaired to its canonical cut (only the
inter-community edges stay cut). That makes each chromosome a unique partition,
so fitness can be cached by chromosome and duplicates are easy to spot.
"""

from __future__ import annotations

import time
from collections import deque
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field

import numpy as np

from .graph import Graph
from .multilevel import Level, multilevel
from .objectives import MOVE_GAINS, Objective, get_objective
from .partition import canonical_cut, decode

# fn(graph, labels, rng) -> new labels
Operator = Callable[[Graph, np.ndarray, np.random.Generator], np.ndarray]


@dataclass
class GAConfig:
    pop_size: int = 100
    generations: int = 200
    elite_frac: float = 0.1
    tournament_size: int = 3
    crossover_rate: float = 0.8
    crossover: str = "uniform"  # "uniform" or "one_point"
    mutation_rate: float = 0.9  # chance a child gets mutated at all
    init: str = "random"  # "random" or "degree" (the 2018 degree-biased start)
    local_search: bool = True  # greedily refine each new best individual (memetic step)
    aggregate: bool = True  # local search also moves whole communities (Louvain-style)
    patience: int | None = 100  # stop after this many generations without improvement
    seed: int | None = None
    # island model: pop_size is split across islands that evolve separately; every
    # migration_interval generations each island's best `migrants` replace the
    # next island's worst (ring topology)
    islands: int = 1
    migration_interval: int = 25
    migrants: int = 1

    def __post_init__(self) -> None:
        checks = [
            (self.pop_size >= 1, "pop_size must be >= 1"),
            (self.islands >= 1, "islands must be >= 1"),
            (self.migration_interval >= 1, "migration_interval must be >= 1"),
            (self.migrants >= 1, "migrants must be >= 1"),
            (self.islands <= 1 or self.pop_size // self.islands > self.migrants,
             "each island needs more than `migrants` members: raise pop_size or lower islands"),
            (self.generations >= 0, "generations must be >= 0"),
            (0 <= self.elite_frac <= 1, "elite_frac must be in [0, 1]"),
            (1 <= self.tournament_size, "tournament_size must be >= 1"),
            (0 <= self.crossover_rate <= 1, "crossover_rate must be in [0, 1]"),
            (0 <= self.mutation_rate <= 1, "mutation_rate must be in [0, 1]"),
            (self.crossover in ("uniform", "one_point"),
             "crossover must be 'uniform' or 'one_point'"),
            (self.init in ("random", "degree"), "init must be 'random' or 'degree'"),
            (self.patience is None or self.patience >= 1, "patience must be >= 1 or None"),
        ]
        for ok, message in checks:
            if not ok:
                raise ValueError(message)


@dataclass
class GAResult:
    fitness: float
    labels: np.ndarray  # community of each node, by node index
    cut: np.ndarray
    communities: dict  # original node id -> community
    history: list = field(default_factory=list)  # (best, mean) per generation
    generations: int = 0
    evaluations: int = 0
    seconds: float = 0.0

    @property
    def num_communities(self) -> int:
        return int(self.labels.max()) + 1


class ERBGA:
    def __init__(self, graph: Graph, objective: str | Objective = "modularity",
                 config: GAConfig | None = None, operators: Sequence[tuple[float, Operator]] = ()):
        """``operators`` are extra mutation operators as (probability, fn) pairs.

        ``fn(graph, labels, rng)`` gets a copy of a parent's community labels and
        returns new labels; communities that end up disconnected are split.
        Each mutation picks an extra operator with its probability, otherwise a
        built-in one.
        """
        if graph.m == 0:
            raise ValueError("graph has no edges")
        if sum(w for w, _ in operators) > 1 or any(w < 0 for w, _ in operators):
            raise ValueError("operator probabilities must be >= 0 and sum to at most 1")
        self.graph = graph
        self.objective = get_objective(objective)
        self.config = config or GAConfig()
        self.operators = list(operators)
        self.rng = np.random.default_rng(self.config.seed)
        self._cache: dict[bytes, float] = {}  # canonical cut -> fitness
        self._base_level: Level | None = None  # graph as a multilevel Level, built lazily
        self.evaluations = 0

    # ---- evaluation --------------------------------------------------

    def evaluate(self, cut: np.ndarray) -> tuple[np.ndarray, float, np.ndarray]:
        """Repair ``cut`` to canonical form; return (cut, fitness, labels)."""
        labels = decode(self.graph, cut)
        cut = canonical_cut(self.graph, labels)
        key = np.packbits(cut).tobytes()
        fitness = self._cache.get(key)
        if fitness is None:
            self.evaluations += 1
            fitness = float(self.objective(self.graph, labels))
            self._cache[key] = fitness
        return cut, fitness, labels

    # ---- initialisation ----------------------------------------------

    def _initial_cut(self) -> np.ndarray:
        g, rng = self.graph, self.rng
        if self.config.init == "degree":
            # Edges between unusual-degree nodes are more likely to start cut,
            # following the 2018 C++ initialiser.
            d = g.degree.astype(float)
            rate = (d.max() - d) / d.max() + np.abs(d.mean() - d) / d.max()
            ru, rv = rate[g.src], rate[g.dst]
            factor = np.where((ru == 0) | (rv == 0), ru + rv, (ru + rv) / 2)
            p = np.clip(rng.uniform(0.0, 0.75) + factor, 0, 1)
        elif self.config.init == "random":
            # Vary the cut density per individual to cover coarse and fine partitions.
            p = rng.uniform(0.05, 0.6)
        else:
            raise ValueError(f"unknown init {self.config.init!r}")
        return rng.random(g.m) < p

    # ---- variation operators -----------------------------------------

    def _tournament(self, fitness: np.ndarray) -> int:
        picks = self.rng.integers(0, len(fitness), self.config.tournament_size)
        return int(picks[np.argmax(fitness[picks])])

    def _crossover(self, a: np.ndarray, b: np.ndarray) -> np.ndarray:
        if self.config.crossover == "uniform":
            return np.where(self.rng.random(len(a)) < 0.5, a, b)
        if self.config.crossover == "one_point":
            if len(a) < 2:
                return a.copy()
            point = self.rng.integers(1, len(a))
            return np.concatenate([a[:point], b[point:]])
        raise ValueError(f"unknown crossover {self.config.crossover!r}")

    def _move_node(self, cut: np.ndarray, labels: np.ndarray) -> None:
        """Move a random node into a neighbour's community (or isolate it)."""
        g, rng = self.graph, self.rng
        v = int(rng.integers(g.n))
        edges, nbrs = g.incident_edges(v), g.neighbours(v)
        if len(edges) == 0:
            return
        cut[edges] = True
        if rng.random() < 0.9:
            target = labels[nbrs[rng.integers(len(nbrs))]]
            cut[edges[labels[nbrs] == target]] = False

    def _merge(self, cut: np.ndarray, labels: np.ndarray) -> None:
        """Merge two adjacent communities by restoring every edge between them."""
        boundary = np.flatnonzero(cut)
        if len(boundary) == 0:
            return
        e = boundary[self.rng.integers(len(boundary))]
        a, b = labels[self.graph.src[e]], labels[self.graph.dst[e]]
        lu, lv = labels[self.graph.src], labels[self.graph.dst]
        cut[((lu == a) & (lv == b)) | ((lu == b) & (lv == a))] = False

    def _split(self, cut: np.ndarray, labels: np.ndarray) -> None:
        """Split a random community in two by growing regions from two seed nodes."""
        g, rng = self.graph, self.rng
        sizes = np.bincount(labels)
        candidates = np.flatnonzero(sizes > 1)
        if len(candidates) == 0:
            return
        c = candidates[rng.integers(len(candidates))]
        members = np.flatnonzero(labels == c)
        seeds = rng.choice(members, 2, replace=False)
        side = {int(seeds[0]): 0, int(seeds[1]): 1}
        queue = deque(int(s) for s in seeds)
        while queue:  # breadth-first, so the two regions grow at the same rate
            v = queue.popleft()
            for w in g.neighbours(v):
                w = int(w)
                if labels[w] == c and w not in side:
                    side[w] = side[v]
                    queue.append(w)
        lu, lv = labels[g.src], labels[g.dst]
        inside = np.flatnonzero((lu == c) & (lv == c))
        for e in inside:
            if side[int(g.src[e])] != side[int(g.dst[e])]:
                cut[e] = True

    def _flip_bits(self, cut: np.ndarray) -> None:
        m = self.graph.m
        flips = self.rng.random(m) < 1.0 / m
        flips[self.rng.integers(m)] = True
        cut ^= flips

    def _mutate(self, cut: np.ndarray) -> np.ndarray:
        if self.rng.random() >= self.config.mutation_rate:
            return cut
        cut, _, labels = self.evaluate(cut)
        cut = cut.copy()
        op = self.rng.random()
        for weight, fn in self.operators:
            if op < weight:
                new_labels = np.asarray(fn(self.graph, labels.copy(), self.rng))
                if (new_labels.shape != labels.shape
                        or not np.issubdtype(new_labels.dtype, np.integer)):
                    raise ValueError(f"operator {getattr(fn, '__name__', fn)!r} must return "
                                     f"an integer array of shape {labels.shape}")
                return canonical_cut(self.graph, new_labels)
            op -= weight
        # rescale so the built-in thresholds below see a uniform draw again
        op /= 1 - sum(w for w, _ in self.operators) or 1
        if op < 0.55:
            for _ in range(self.rng.geometric(0.5)):
                self._move_node(cut, labels)
        elif op < 0.75:
            self._split(cut, labels)
        elif op < 0.9:
            self._merge(cut, labels)
        else:
            self._flip_bits(cut)
        return cut

    def refine(self, cut: np.ndarray) -> tuple[np.ndarray, float]:
        """Greedy local search from ``cut``; never returns anything worse.

        Objectives in ``MOVE_GAINS`` get Louvain-style multilevel search:
        single-node moves scored in O(degree), then (with ``aggregate``) moves
        of whole communities on the collapsed graph. The result is decoded once
        per round, splitting any community a move disconnected. Other
        objectives fall back to a full evaluation per single-node move.
        """
        gain = MOVE_GAINS.get(self.objective)
        if gain is None:
            return self._refine_exact(cut)
        g = self.graph
        if self._base_level is None:
            self._base_level = Level.from_graph(g)
        cut, best, labels = self.evaluate(cut)
        while True:
            new_labels = multilevel(self._base_level, labels, gain, g.m, g.n, self.rng,
                                    aggregate=self.config.aggregate)
            # decoding splits any community a move disconnected; that only adds gain
            t_cut, t_fit, t_labels = self.evaluate(canonical_cut(g, new_labels))
            if t_fit <= best + 1e-12:
                return cut, best
            cut, best, labels = t_cut, t_fit, t_labels

    def _refine_exact(self, cut: np.ndarray) -> tuple[np.ndarray, float]:
        """``refine`` for objectives without a move gain: evaluate every move."""
        g = self.graph
        cut, best, labels = self.evaluate(cut)
        improved = True
        while improved:
            improved = False
            for v in self.rng.permutation(g.n):
                for target in np.unique(labels[g.neighbours(v)]):
                    if target == labels[v]:
                        continue
                    trial = labels.copy()
                    trial[v] = target
                    t_cut, t_fit, t_labels = self.evaluate(canonical_cut(g, trial))
                    if t_fit > best + 1e-12:
                        cut, best, labels, improved = t_cut, t_fit, t_labels, True
        return cut, best

    # ---- main loop ---------------------------------------------------

    def _generation(self, pop: np.ndarray, fit: np.ndarray, refined: set[bytes]):
        """One generation of one island: refine its best, keep elites, breed the rest."""
        cfg = self.config
        size = len(pop)
        n_elite = max(1, int(round(cfg.elite_frac * size)))
        if cfg.local_search:
            top = int(np.argmax(fit))
            key = np.packbits(pop[top]).tobytes()
            if key not in refined:
                pop[top], fit[top] = self.refine(pop[top])
                refined.add(key)
                refined.add(np.packbits(pop[top]).tobytes())
        order = np.argsort(-fit, kind="stable")
        new_pop, new_fit = list(pop[order[:n_elite]]), list(fit[order[:n_elite]])
        while len(new_pop) < size:
            a = pop[self._tournament(fit)]
            if self.rng.random() < cfg.crossover_rate:
                child = self._crossover(a, pop[self._tournament(fit)])
            else:
                child = a.copy()
            cut, f, _ = self.evaluate(self._mutate(child))
            new_pop.append(cut)
            new_fit.append(f)
        return np.array(new_pop), np.array(new_fit)

    def _migrate(self, islands: list) -> None:
        """Ring migration: island i's best ``migrants`` replace island i+1's worst."""
        k = self.config.migrants
        emigrants = []
        for pop, fit in islands:
            best = np.argsort(-fit, kind="stable")[:k]
            emigrants.append((pop[best].copy(), fit[best].copy()))
        for i, (pop, fit) in enumerate(islands):
            src_pop, src_fit = emigrants[i - 1]  # i - 1 wraps to the last island
            worst = np.argsort(fit, kind="stable")[:k]
            pop[worst], fit[worst] = src_pop, src_fit

    def run(self, callback=None) -> GAResult:
        cfg, start = self.config, time.perf_counter()
        base, extra = divmod(cfg.pop_size, cfg.islands)
        sizes = [base + (i < extra) for i in range(cfg.islands)]

        islands = []
        for size in sizes:
            pop, fit = [], []
            for _ in range(size):
                cut, f, _ = self.evaluate(self._initial_cut())
                pop.append(cut)
                fit.append(f)
            islands.append((np.array(pop), np.array(fit)))

        def stats():
            all_fit = np.concatenate([fit for _, fit in islands])
            return float(all_fit.max()), float(all_fit.mean())

        refined: set[bytes] = set()
        history = [stats()]
        best, stale, gen = history[0][0], 0, 0
        for gen in range(1, cfg.generations + 1):
            islands = [self._generation(pop, fit, refined) for pop, fit in islands]
            if cfg.islands > 1 and gen % cfg.migration_interval == 0:
                self._migrate(islands)
            history.append(stats())
            if callback:
                callback(gen, history[-1])

            if history[-1][0] > best + 1e-12:
                best, stale = history[-1][0], 0
            else:
                stale += 1
                if cfg.patience is not None and stale >= cfg.patience:
                    break

        pop, fit = max(islands, key=lambda island: island[1].max())
        top = int(np.argmax(fit))
        cut, f, labels = self.evaluate(pop[top])
        return GAResult(
            fitness=float(f),
            labels=labels,
            cut=cut,
            communities={_plain(nid): int(c)
                         for nid, c in zip(self.graph.node_ids, labels, strict=True)},
            history=history,
            generations=gen,
            evaluations=self.evaluations,
            seconds=time.perf_counter() - start,
        )


def detect_communities(graph: Graph, objective: str | Objective = "modularity",
                       **config) -> GAResult:
    """Convenience wrapper: ``detect_communities(g, "modularity", seed=0)``."""
    return ERBGA(graph, objective, GAConfig(**config)).run()


def _plain(x):
    return x.item() if hasattr(x, "item") else x
