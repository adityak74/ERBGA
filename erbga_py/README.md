# ERBGA (Python)

Connected, objective-agnostic community detection with a genetic algorithm.
This is a clean rewrite of the 2017–2018 C++ ERBGA in `../bfs_orig`.

## The idea

A chromosome has **one bit per edge**: `1` means the edge is cut. Communities
are the connected components left after the cut edges are removed. Two things
follow:

- **Every community is connected by construction.** Louvain does not give that
  guarantee; it was Leiden's headline fix in 2019.
- **The fitness function is any callable** `f(graph, labels) -> float`. Nothing
  in the search assumes modularity. `modularity` and `qs` (the density-based
  score from the thesis) are built in, and you can pass your own.

## What changed from the C++ version

| | C++ (2018) | Python (this) |
|---|---|---|
| Decoding | BFS, plus a shell process and a file write on every evaluation | vectorised union-find (numpy), no I/O |
| Encoding | any bit string; many strings give the same partition | **canonical repair**: only edges between communities stay cut, so one chromosome = one partition, and fitness is cached by chromosome |
| Mutation | flips each bit with probability 0.5 (random search) | node move, community split, community merge, rare bit flips |
| Gene repair | puts back cut edges at high-degree nodes | replaced by canonical repair; optional greedy node-move refinement of the best individual each generation (`local_search`) |
| Selection | tournament that always includes an elite | standard tournament + elitism; best fitness never decreases (tested) |
| Reproducibility | `srand(time)` | `seed=` gives identical runs (tested) |
| Karate, best modularity | 0.342 after 100 gens, 28 s | **0.4198** (known optimum), < 1 s |

## Install and run

```bash
cd erbga_py
pip install -e ".[test]"

erbga data/karate.gml --seed 0
erbga data/football.gml --objective qs --pop 200 --gen 300 --json result.json -v
python -m erbga --help          # same thing without installing the script
```

```python
from erbga import read_gml, detect_communities

g = read_gml("data/karate.gml")
r = detect_communities(g, "modularity", seed=0)
r.fitness, r.num_communities, r.communities   # {node_id: community}

# any objective works
def my_score(graph, labels):
    ...
detect_communities(g, my_score, seed=0, pop_size=200)
```

## Benchmark

`python benchmarks/run_datasets.py --seeds 10`. Modularity, 10 seeds per row,
default settings (population 100, at most 200 generations, stop after 50
without improvement). "Hit optimum" counts runs within 5e-5 of the proven
maximum (Aloise et al., 2010).

| graph | n | m | variant | best | mean ± sd | hit optimum | mean time (s) |
|---|---|---|---|---|---|---|---|
| karate | 34 | 78 | GA + local search | 0.4198 | 0.4198 ± 0.0000 | 10/10 | 0.70 |
| karate | 34 | 78 | plain GA | 0.4198 | 0.4198 ± 0.0000 | 10/10 | 0.91 |
| dolphins | 62 | 159 | GA + local search | 0.5285 | 0.5265 ± 0.0026 | 2/10 | 1.03 |
| dolphins | 62 | 159 | plain GA | 0.5285 | 0.5276 ± 0.0010 | 5/10 | 1.78 |
| polbooks | 105 | 441 | GA + local search | 0.5272 | 0.5263 ± 0.0025 | 5/10 | 1.20 |
| polbooks | 105 | 441 | plain GA | 0.5272 | 0.5272 ± 0.0001 | 8/10 | 2.50 |
| football | 115 | 613 | GA + local search | 0.6046 | 0.6045 ± 0.0001 | 5/10 | 2.86 |
| football | 115 | 613 | plain GA | 0.6046 | 0.6004 ± 0.0033 | 1/10 | 3.85 |

Local search is faster, and better on football, but it converges early and
sometimes gets stuck (dolphins, polbooks). Balancing the two is open.

## Layout

```
erbga/graph.py       Graph (CSR adjacency) and GML reader
erbga/partition.py   union-find decoding, canonical cut, connectivity check
erbga/objectives.py  modularity, qs, objective registry
erbga/ga.py          GAConfig, ERBGA (operators + main loop), detect_communities
erbga/cli.py         command line
tests/               pytest suite (modularity checked against networkx)
benchmarks/          dataset benchmark
data/                karate, dolphins, polbooks, football
```

## Next steps

- **Speed.** Objectives are recomputed from scratch on every evaluation, and
  local search evaluates every candidate move that way. Incremental (delta)
  updates for modularity and Qs would make refinement about n times cheaper,
  which is needed past a few thousand nodes.
- **LFR benchmarks** with NMI/ARI against Louvain, Leiden and Infomap.
- **Multi-objective** search (NSGA-II over modularity and Qs).
- **Weighted graphs.** The GML reader ignores edge weights for now.
- **Qs definition.** The C++ `q_calc` skipped two-node communities; this
  version counts them. Check which matches the intended definition.
