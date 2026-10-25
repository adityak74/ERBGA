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
from dataclasses import dataclass, field

import numpy as np

from .graph import Graph
from .objectives import Objective, get_objective
from .partition import canonical_cut, decode


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
    patience: int | None = 50  # stop after this many generations without improvement
    seed: int | None = None


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
                 config: GAConfig | None = None):
        if graph.m == 0:
            raise ValueError("graph has no edges")
        self.graph = graph
        self.objective = get_objective(objective)
        self.config = config or GAConfig()
        self.rng = np.random.default_rng(self.config.seed)
        self._cache: dict[bytes, tuple[float, np.ndarray]] = {}
        self.evaluations = 0

    # ---- evaluation --------------------------------------------------

    def evaluate(self, cut: np.ndarray) -> tuple[np.ndarray, float, np.ndarray]:
        """Repair ``cut`` to canonical form; return (cut, fitness, labels)."""
        labels = decode(self.graph, cut)
        cut = canonical_cut(self.graph, labels)
        key = np.packbits(cut).tobytes()
        hit = self._cache.get(key)
        if hit is None:
            self.evaluations += 1
            hit = (self.objective(self.graph, labels), labels)
            self._cache[key] = hit
        return cut, hit[0], hit[1]

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
        """Greedy local search: move single nodes to neighbouring communities
        while the objective improves. Works for any objective, at the price of
        a full objective evaluation per candidate move."""
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

    def run(self, callback=None) -> GAResult:
        cfg, start = self.config, time.perf_counter()
        n_elite = max(1, int(round(cfg.elite_frac * cfg.pop_size)))

        pop, fit = [], []
        for _ in range(cfg.pop_size):
            cut, f, _ = self.evaluate(self._initial_cut())
            pop.append(cut)
            fit.append(f)
        pop, fit = np.array(pop), np.array(fit)

        refined: set[bytes] = set()
        history = [(float(fit.max()), float(fit.mean()))]
        best, stale, gen = fit.max(), 0, 0
        for gen in range(1, cfg.generations + 1):
            if cfg.local_search:
                top = int(np.argmax(fit))
                key = np.packbits(pop[top]).tobytes()
                if key not in refined:
                    pop[top], fit[top] = self.refine(pop[top])
                    refined.add(key)
                    refined.add(np.packbits(pop[top]).tobytes())
            order = np.argsort(-fit, kind="stable")
            new_pop, new_fit = list(pop[order[:n_elite]]), list(fit[order[:n_elite]])
            while len(new_pop) < cfg.pop_size:
                a = pop[self._tournament(fit)]
                if self.rng.random() < cfg.crossover_rate:
                    child = self._crossover(a, pop[self._tournament(fit)])
                else:
                    child = a.copy()
                cut, f, _ = self.evaluate(self._mutate(child))
                new_pop.append(cut)
                new_fit.append(f)
            pop, fit = np.array(new_pop), np.array(new_fit)
            history.append((float(fit.max()), float(fit.mean())))
            if callback:
                callback(gen, history[-1])

            if fit.max() > best + 1e-12:
                best, stale = fit.max(), 0
            else:
                stale += 1
                if cfg.patience is not None and stale >= cfg.patience:
                    break

        top = int(np.argmax(fit))
        cut, f, labels = self.evaluate(pop[top])
        return GAResult(
            fitness=float(f),
            labels=labels,
            cut=cut,
            communities={_plain(nid): int(c) for nid, c in zip(self.graph.node_ids, labels)},
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
