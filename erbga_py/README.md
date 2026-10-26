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

erbga detect karate --seed 0                 # bundled datasets work by name
erbga detect football --objective qs --pop 200 --gen 300 --json result.json -v
erbga detect my_graph.txt --format json      # edge list, .gml or .json
erbga summary polbooks
python -m erbga --help                       # same thing without installing the script
```

```python
from erbga import read_gml, detect_communities

g = read_gml("erbga/data/karate.gml")
r = detect_communities(g, "modularity", seed=0)
r.fitness, r.num_communities, r.communities   # {node_id: community}

# any objective works
def my_score(graph, labels):
    ...
detect_communities(g, my_score, seed=0, pop_size=200)
```

## Use with LLMs

Full reference for models: [`llms.txt`](llms.txt). Guidance for coding agents
working on this repo: [`AGENTS.md`](AGENTS.md).

**MCP server.** `pip install -e ".[mcp]"`, then add it to your client:

```json
{"mcpServers": {"erbga": {"command": "erbga", "args": ["mcp"]}}}
```

Tools: `list_datasets`, `list_objectives`, `graph_summary`, `detect_communities`,
`score_partition`. A graph can be a dataset name, a file path, or inline
`{"edges": [[u, v], ...]}`. Errors come back as
`{"ok": false, "error": {"code", "message", "hint"}}` so the model can fix the
input and retry.

**Agent-readable CLI.** Every subcommand takes `--format json`. stdout is then
one JSON object (the result or the error), progress goes to stderr, and the exit
code is 0 or 1. `--max-nodes N` keeps output small on big graphs.

**LLM-evolved operators.** `erbga evolve` runs a FunSearch-style loop. Claude
writes `mutate(graph, labels, rng)` operators; each is checked against an
import allowlist, then scored in a subprocess with time and memory limits
against a fixed-seed baseline. The best is saved to `best_operator.py`:

```bash
pip install -e ".[evolve]"     # needs ANTHROPIC_API_KEY or `ant auth login`
erbga evolve --iterations 20 --datasets karate dolphins polbooks
erbga detect football --operator evolve_runs/<run>/best_operator.py
```

It uses `claude-opus-5-5` at effort `high` with server-side refusal fallback on
(`--model` to change). This executes LLM-written code. The allowlist and limits
are not a security sandbox, so run it in a container if that matters to you.

## Benchmark

`python benchmarks/run_datasets.py --seeds 10`. Modularity, 10 seeds per row,
default settings (population 100, at most 200 generations, stop after 100
without improvement). "Hit optimum" counts runs within 5e-5 of the proven
maximum (Aloise et al., 2010).

| graph | n | m | variant | best | mean ± sd | hit optimum | mean time (s) |
|---|---|---|---|---|---|---|---|
| karate | 34 | 78 | GA + local search | 0.4198 | 0.4198 ± 0.0000 | 10/10 | 1.44 |
| karate | 34 | 78 | plain GA | 0.4198 | 0.4198 ± 0.0000 | 10/10 | 1.80 |
| dolphins | 62 | 159 | GA + local search | 0.5285 | 0.5274 ± 0.0007 | 2/10 | 2.37 |
| dolphins | 62 | 159 | plain GA | 0.5285 | 0.5276 ± 0.0010 | 5/10 | 2.93 |
| polbooks | 105 | 441 | GA + local search | 0.5272 | 0.5272 ± 0.0002 | 8/10 | 2.20 |
| polbooks | 105 | 441 | plain GA | 0.5272 | 0.5272 ± 0.0000 | 10/10 | 3.61 |
| football | 115 | 613 | GA + local search | 0.6046 | 0.6045 ± 0.0001 | 6/10 | 3.99 |
| football | 115 | 613 | plain GA | 0.6046 | 0.6004 ± 0.0033 | 1/10 | 4.45 |

On these small graphs the plain GA, which runs without local search, hits the
optimum more often but takes longer. Local search reaches a good partition
quickly and then converges, which can stop it just short of the optimum
(dolphins). On large graphs local search is what makes ERBGA practical; see Scaling below.

### Scaling

`python benchmarks/scaling.py`: planted-partition graphs with 10 equal blocks
and average degree about 7, default settings, one seed. "Planted Q" is the
modularity of the blocks the graph was generated from.

| n | m | planted Q | best Q | communities | generations | seconds |
|---|---|---|---|---|---|---|
| 500 | 1668 | 0.6363 | 0.6400 | 10 | 107 | 4.6 |
| 1000 | 3407 | 0.6290 | 0.6323 | 10 | 111 | 8.7 |
| 2000 | 6857 | 0.6374 | 0.6411 | 14 | 116 | 17.2 |
| 5000 | 16961 | 0.6344 | 0.6370 | 17 | 101 | 32.8 |

Two things make this fast:

- **Incremental move gains.** Each candidate move is scored in O(degree) from
  per-community totals (`MOVE_GAINS` in `objectives.py`) instead of
  re-decoding the graph. On the 2000-node graph at 30 generations, that cut
  local search from 360 s to 9 s; it had been 95% of the runtime.
- **Multilevel (Louvain-style) local search** (`multilevel.py`). After
  single-node moves stall, each community is collapsed into one weighted node
  and moves continue on the smaller graph, so a single move merges two whole
  communities. Starting from all singletons, node moves alone stall far below
  the planted partition, and aggregation gets past it:

  | n | planted Q | node moves only | with aggregation |
  |---|---|---|---|
  | 1000 | 0.6290 | 0.2760 (324 communities) | 0.6323 (10), 0.5 s |
  | 2000 | 0.6374 | 0.2485 (690) | 0.6408 (14), 1.7 s |
  | 5000 | 0.6344 | 0.1788 (2078) | 0.6382 (17), 6.5 s |

  The result is mapped back to the original nodes and decoded, so a move that
  disconnects a community makes the decoder split it, and the connectivity
  guarantee holds. `--no-aggregate` (or `GAConfig(aggregate=False)`) turns
  aggregation off.

### Leiden-style refinement (opt-in)

`--refinement` (or `GAConfig(refinement=True)`) adds Leiden's refinement
step to multilevel local search. Before each aggregation, every community is
rebuilt from singletons by merging each unit into an adjacent sub-community
of the *same* community, and those connected sub-communities, not whole
communities, become the next level's units. Each unit starts the next level
in its parent community, so later moves can shift part of a community
instead of only merging whole ones.

As in Leiden, a unit picks among the merges that don't lower the score at
random, with probability proportional to exp(gain / θ), where gains are in
edge units and θ = `refinement_theta` (default 0.01, Leiden's default, set
with `--refinement-theta`). θ = 0 always takes the best merge.

Measured with default settings otherwise:

| | Louvain (default) | greedy refinement (θ = 0) | **randomized refinement (θ = 0.01)** |
|---|---|---|---|
| Optimum hits, 4 datasets × 10 seeds | 26 | 27 | **28** |
| Football: hits, worst seed | 6/10, 0.6044 | 5/10, 0.5972 | **7/10, 0.6044** |
| Dolphins: hits, worst seed | 2/10, **0.5265** | 4/10, 0.5201 | 3/10, 0.5201 |
| 1,000 / 2,000 / 5,000 nodes, mean of 3 seeds | 0.6322 / 0.6407 / 0.6376 | 0.6325 / 0.6413 / 0.6384 | **0.6325 / 0.6413 / 0.6384** |
| 1,000 / 2,000 / 5,000 nodes, mean time | 11.9 / 20.7 / 41.1 s | 9.2 / 18.2 / 43.0 s | 10.4 / 18.4 / 55.5 s |

The randomized merge fixes greedy refinement's football trap (2 of 10 seeds
had settled at 0.5972, with two conferences merged). It gives the most
optimum hits and the best large-graph scores. One dolphins seed in ten still
lands at 0.5201, which the Louvain default never does, so refinement stays
opt-in. Use it for graphs of 1,000+ nodes, or run a few seeds and keep the
best. (Choosing Louvain or Leiden at random on each local search call did
worse than either: 25 hits.)

### Island model

`GAConfig(islands=4)` (or `erbga detect GRAPH --islands 4`) splits the
population into islands that evolve separately. Every `migration_interval`
generations (default 25), each island's best `migrants` (default 1) replace
the next island's worst, in a ring. `pop_size` is the total across islands,
so turning islands on doesn't change the budget. `islands=1`, the default,
reproduces the single-population runs exactly.

Islands help when one population converges too early, and can cost a little
where it doesn't. Optimum hits out of 10 seeds, population 100, no early
stopping. These were measured before multilevel local search was added:

| | dolphins | football | mean time |
|---|---|---|---|
| 1 island | 2/10 | 8/10 | 3.2 s / 3.8 s |
| 4 islands | **4/10** | 5/10 | 3.4 s / 4.2 s |
| 4 islands, 2 migrants every 10 generations | **4/10** | 7/10 | 3.4 s / 4.1 s |

With default early stopping, summed over all four datasets (karate,
dolphins, polbooks, football), the score is 25 hits for 1 island and 24–25
for 2–5 islands: islands move hits from football to dolphins rather than
adding them. Try islands when repeated seeds stall at slightly different
scores. This ports the island model from the 2018 `multiple-subpops`
branch, without its bugs.

## Layout

```
erbga/graph.py       Graph (CSR adjacency) and GML reader
erbga/partition.py   union-find decoding, canonical cut, connectivity check
erbga/objectives.py  modularity, qs, objective registry
erbga/ga.py          GAConfig, ERBGA (operators + main loop), detect_communities
erbga/api.py         JSON-in/JSON-out functions, structured errors, graph formats
erbga/cli.py         command line (detect, score, summary, objectives, mcp, evolve)
erbga/mcp_server.py  MCP tools for agents
erbga/evolve.py      LLM-evolved mutation operators
erbga/data/          karate, dolphins, polbooks, football
tests/               pytest suite (modularity checked against networkx)
benchmarks/          dataset benchmark
llms.txt, AGENTS.md  docs for models using / agents editing the package
```

## Next steps

- **Adaptive islands.** Islands help on dolphins and hurt on football. Varying
  the operator mix or the migration rate per island might keep the gain
  without the cost.

- **Speed past ~10k nodes.** Local search is now incremental, but each GA
  evaluation still decodes the whole chromosome (O(m) numpy work), and that is
  what dominates at scale. Options: update components incrementally for
  single-node moves, or move the hot loops to compiled code.
- **Leiden refinement on small graphs.** With the randomized merge, `--refinement`
  gives the most optimum hits overall, but 1 dolphins seed in 10 settles at
  0.5201. Leiden's well-connectedness condition (merge only into
  sub-communities well connected to the rest of their community) is the
  untried part of Leiden that might close this gap.
- **Multilevel search for custom objectives.** Aggregated moves currently need a
  move gain; custom objectives use single-node moves with full evaluation.
- **LFR benchmarks** with NMI/ARI against Louvain, Leiden and Infomap.
- **Multi-objective** search (NSGA-II over modularity and Qs).
- **Weighted graphs.** The GML reader ignores edge weights for now.
- **Qs definition.** The C++ `q_calc` skipped two-node communities; this
  version counts them. Check which matches the intended definition.
