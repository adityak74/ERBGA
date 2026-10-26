# AGENTS.md — working on erbga_py

Guidance for coding agents (Claude Code, Codex, Cursor, …) changing this package.
For *using* ERBGA as a tool, see `llms.txt` and the MCP server instead.

## Commands

Run from `erbga_py/`:

```bash
pip install -e ".[test]"          # numpy + pytest, networkx, mcp
pytest -q                          # full suite, ~15 s
ruff check .                       # lint (config in pyproject.toml)
python benchmarks/run_datasets.py --seeds 10   # quality benchmark, ~1 min
python benchmarks/scaling.py                   # 500-2000 node graphs, ~1 min
erbga detect karate --seed 0       # smoke test; must print modularity 0.419790
```

CI (`.github/workflows/erbga-py.yml`) runs ruff and pytest on Python 3.10–3.13.
A change is done when both pass.

## Layout

| File | Role |
|---|---|
| `erbga/graph.py` | `Graph` (edges `src < dst`, CSR adjacency), `read_gml` |
| `erbga/partition.py` | decode a cut into communities (vectorised union-find), `canonical_cut` |
| `erbga/objectives.py` | `modularity`, `qs`, the `OBJECTIVES` registry |
| `erbga/ga.py` | `GAConfig`, `ERBGA` (operators, main loop), `detect_communities` |
| `erbga/api.py` | JSON-in/JSON-out layer, `ERBGAError` codes, graph input formats |
| `erbga/cli.py` | `erbga` command (subcommands, `--format json`) |
| `erbga/mcp_server.py` | MCP tools; thin wrappers over `api` |
| `erbga/evolve.py` | LLM-evolved mutation operators (FunSearch-style loop) |
| `erbga/data/` | bundled GML datasets (karate, dolphins, polbooks, football) |

Dependencies go one way: `graph` ← `partition` ← `objectives` ← `ga` ← `api` ←
`cli` / `mcp_server` / `evolve`. Only numpy is required at runtime; `mcp` and
`anthropic` are optional extras and must stay lazily imported outside their modules.

## Invariants (tests enforce these; don't break them)

- A chromosome is one bool per edge, `True` = cut. Communities are the connected
  components after removing cut edges, so **every community is connected**.
- `ERBGA.refine` never returns a worse fitness than its input, and the fitness
  it returns equals `evaluate` on the returned cut.
- `ERBGA.evaluate` always returns the **canonical cut** (only edges between
  communities are cut). The fitness cache is keyed on it.
- Elitism: the best fitness never decreases from one generation to the next.
- `islands=1` reproduces the single-population run exactly (same random draws);
  island code must not consume randomness when there is only one island.
- Same `seed` ⇒ identical run. Draw all randomness from `self.rng`. Adding
  operators must not change the random stream when they are not used: see the
  rescaling in `ERBGA._mutate`.
- Objectives are `f(graph, labels) -> float`, higher is better, with labels in
  `0..k-1`.
- `api` functions return JSON-serialisable dicts with `"ok": true`, and raise
  `ERBGAError(code, message, hint)` on bad input. Never let a raw traceback reach
  CLI `--format json` or MCP callers. Add new error codes to the `ERBGAError` docstring.

## Conventions

- Match the surrounding style: short docstrings that say *why*, numpy vectorised
  where it is on the hot path (`decode`, objectives), plain Python elsewhere.
- New objective: add a function to `objectives.py`, register it in `OBJECTIVES`,
  write a hand-computed test in `tests/test_objectives.py`, and document it in `llms.txt`.
  For fast local search, also add a move gain to `MOVE_GAINS` (it must equal the
  exact objective change for a single-node move on the undecoded labels, and the
  objective must not decrease when a disconnected community is split, so the gain
  is a lower bound after decoding). `test_move_gain_matches_full_recompute` checks
  the first property. Without a gain, `refine` falls back to full evaluation per move.
- New GA option: add it to `GAConfig` with validation in `__post_init__`, a CLI flag
  in `cli.py`, and (if agents should set it) a parameter of the MCP `detect_communities`.
- Keep MCP tool docstrings accurate. The model reads them to decide what to call.
- `evolve.py` runs LLM-written code. Keep the allowlist (`check_candidate`), the
  subprocess, and the resource limits intact. Never `exec` candidate code in the
  parent process (`load_operator` is only for files the user chose to trust).

## Known gaps

- Each GA evaluation decodes the whole chromosome (O(m)); that dominates past a
  few thousand nodes. Local search is already incremental (see `MOVE_GAINS` under Conventions).
- Edge weights in GML input are ignored.
- The C++ `q_calc` skipped two-node communities in Qs; `qs` here counts them. Which
  definition is intended is still open.
