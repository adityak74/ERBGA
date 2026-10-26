"""Multilevel (Louvain-style) local search over a partition.

Local moving alone gets stuck once every single-node move would lower the
score, even when merging two whole communities would raise it. Aggregation
fixes that: collapse each community into one weighted node and run local
moving again, so a single move now relocates a whole community.

Everything here works on labels and returns labels; ``ERBGA.refine`` turns
the result back into a canonical cut, and decoding splits any community that
ended up disconnected, so ERBGA's connectivity guarantee is unchanged.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .graph import Graph


@dataclass
class Level:
    """A (possibly aggregated) graph: unit i stands for ``size[i]`` original
    nodes with summed degree ``degree[i]``. Edges are stored once per direction
    in CSR form with weights; edges inside a unit are dropped, since no gain
    depends on them."""

    degree: np.ndarray
    size: np.ndarray
    ptr: np.ndarray
    nbr: np.ndarray
    weight: np.ndarray

    @property
    def n(self) -> int:
        return len(self.degree)

    @classmethod
    def from_graph(cls, graph: Graph) -> Level:
        return cls(graph.degree.astype(float), np.ones(graph.n), graph.incident_ptr,
                   graph.incident_nbr, np.ones(len(graph.incident_nbr)))

    def aggregate(self, part: np.ndarray) -> tuple[Level, np.ndarray]:
        """Collapse each group of ``part`` into one unit.

        Returns the new level and the compacted group id of every unit here.
        """
        groups = np.unique(part, return_inverse=True)[1]
        k = int(groups.max()) + 1
        src = np.repeat(groups, np.diff(self.ptr))
        dst = groups[self.nbr]
        between = src != dst
        key = src[between] * k + dst[between]
        keys, inverse = np.unique(key, return_inverse=True)
        weight = np.bincount(inverse, weights=self.weight[between])
        a = keys // k
        level = Level(
            degree=np.bincount(groups, weights=self.degree, minlength=k),
            size=np.bincount(groups, weights=self.size, minlength=k),
            ptr=np.concatenate([[0], np.cumsum(np.bincount(a, minlength=k))]),
            nbr=keys % k,
            weight=weight,
        )
        return level, groups


def local_moving(level: Level, part: np.ndarray, gain, m: int, n: int,
                 rng: np.random.Generator) -> tuple[np.ndarray, bool]:
    """Move units to the neighbouring community with the best positive gain,
    pass after pass, until a pass moves nothing. Returns (part, moved_any)."""
    part = part.copy()
    k = int(part.max()) + 1
    comm_degree = np.bincount(part, weights=level.degree, minlength=k)
    comm_size = np.bincount(part, weights=level.size, minlength=k)
    moved_any = False
    while True:
        moved = False
        for v in rng.permutation(level.n):
            lo, hi = level.ptr[v], level.ptr[v + 1]
            if lo == hi:
                continue
            targets, inverse = np.unique(part[level.nbr[lo:hi]], return_inverse=True)
            links = np.bincount(inverse, weights=level.weight[lo:hi])
            source = part[v]
            own = targets == source
            if own.all():
                continue
            gains = gain(m, n, level.degree[v], level.size[v], links[~own],
                         links[own].sum(), comm_degree[targets[~own]],
                         comm_size[targets[~own]], comm_degree[source] - level.degree[v],
                         comm_size[source] - level.size[v])
            i = int(np.argmax(gains))
            if gains[i] > 1e-12:
                target = targets[~own][i]
                comm_degree[source] -= level.degree[v]
                comm_size[source] -= level.size[v]
                comm_degree[target] += level.degree[v]
                comm_size[target] += level.size[v]
                part[v] = target
                moved = moved_any = True
        if not moved:
            return part, moved_any


def multilevel(base: Level, labels: np.ndarray, gain, m: int, n: int,
               rng: np.random.Generator, aggregate: bool = True) -> np.ndarray:
    """Louvain-style search starting from ``labels`` (one per original node).

    Level 0 moves single nodes. With ``aggregate``, each later level treats the
    previous level's communities as units, so moves merge whole communities;
    it stops at the first coarser level where nothing moves.
    """
    level, part = base, labels.copy()
    membership = np.arange(base.n)  # original node -> unit at the current level
    depth = 0
    while True:
        part, moved = local_moving(level, part, gain, m, n, rng)
        if (not moved and depth > 0) or not aggregate:
            return part[membership]
        level, groups = level.aggregate(part)
        membership = groups[membership]
        part = np.arange(level.n)
        depth += 1
