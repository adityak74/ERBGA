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
default settings (population 100, at most 200 generations, stop after 50
without improvement). "Hit optimum" counts runs within 5e-5 of the proven
maximum (Aloise et al., 2010).

| graph | n | m | variant | best | mean ± sd | hit optimum | mean time (s) |
|---|---|---|---|---|---|---|---|
| karate | 34 | 78 | GA + local search | 0.4198 | 0.4198 ± 0.0000 | 10/10 | 0.83 |
| karate | 34 | 78 | plain GA | 0.4198 | 0.4198 ± 0.0000 | 10/10 | 1.01 |
| dolphins | 62 | 159 | GA + local search | 0.5285 | 0.5272 ± 0.0007 | 1/10 | 1.26 |
| dolphins | 62 | 159 | plain GA | 0.5285 | 0.5276 ± 0.0010 | 5/10 | 2.06 |
| polbooks | 105 | 441 | GA + local search | 0.5272 | 0.5270 ± 0.0006 | 7/10 | 1.29 |
| polbooks | 105 | 441 | plain GA | 0.5272 | 0.5272 ± 0.0001 | 8/10 | 2.69 |
| football | 115 | 613 | GA + local search | 0.6046 | 0.6045 ± 0.0001 | 7/10 | 3.11 |
| football | 115 | 613 | plain GA | 0.6046 | 0.6004 ± 0.0033 | 1/10 | 4.36 |

Local search is faster, and better on football, but it converges early and
can stop short of the optimum (dolphins). Balancing the two is open.

### Scaling

`python benchmarks/scaling.py`: planted-partition graphs with 10 equal blocks
and average degree about 7, default settings, one seed. "Planted Q" is the
modularity of the blocks the graph was generated from.

| n | m | planted Q | best Q | communities | generations | seconds |
|---|---|---|---|---|---|---|
| 500 | 1668 | 0.6363 | 0.6400 | 10 | 112 | 5.2 |
| 1000 | 3407 | 0.6290 | 0.6313 | 10 | 200 | 15.6 |
| 2000 | 6857 | 0.6374 | 0.6404 | 16 | 200 | 36.2 |

Local search scores each single-node move in O(degree) from per-community
totals (`MOVE_GAINS` in `objectives.py`) instead of re-decoding the graph.
At 30 generations on the 2000-node graph that took local search from 360 s
to 9 s; it was 95% of the runtime before.

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

- **Speed past ~10k nodes.** Local search is now incremental, but each GA
  evaluation still decodes the whole chromosome (O(m) numpy work), and that is
  what dominates at scale. Options: update components incrementally for
  single-node moves, or move the hot loops to compiled code.
- **Coarsening.** Single-node moves cannot merge whole communities, so local
  search alone stalls on fine partitions; the GA's merge operator does that
  work today. A Louvain-style aggregation step would do it directly.
- **LFR benchmarks** with NMI/ARI against Louvain, Leiden and Infomap.
- **Multi-objective** search (NSGA-II over modularity and Qs).
- **Weighted graphs.** The GML reader ignores edge weights for now.
- **Qs definition.** The C++ `q_calc` skipped two-node communities; this
  version counts them. Check which matches the intended definition.
