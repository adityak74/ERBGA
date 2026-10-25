"""LLM-evolved mutation operators (FunSearch-style).

Loop: show Claude the best operators found so far with their scores, ask for a
better one, check it statically, score it in a separate process, keep the best.

An operator is ``mutate(graph, labels, rng) -> labels``: it gets a copy of a
parent's community labels and returns new ones. It runs inside the GA next to
the built-in operators (with probability ``operator_weight`` per mutation).

A candidate is scored by running the GA with a small budget on fixed datasets
and seeds and averaging the best objective value reached. The baseline is the
same runs without the candidate, so the improvement is a paired comparison.

SAFETY: this executes code written by an LLM. Candidates must pass an import
and builtin allowlist, and run in a subprocess with a timeout and memory/CPU
limits, but that is not a security sandbox. Run it in a container or VM when
in doubt.
"""

from __future__ import annotations

import ast
import json
import os
import re
import subprocess
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import TextIO
from collections.abc import Callable

DEFAULT_MODEL = "claude-opus-5-5"

ALLOWED_IMPORTS = {"numpy", "math", "random", "collections", "itertools", "heapq",
                   "functools", "operator"}
FORBIDDEN_NAMES = {"open", "exec", "eval", "compile", "__import__", "globals", "locals",
                   "vars", "input", "breakpoint", "exit", "quit", "help", "memoryview",
                   "getattr", "setattr", "delattr"}
FORBIDDEN_ATTRS = {"load", "save", "savez", "savez_compressed", "savetxt", "loadtxt",
                   "genfromtxt", "fromfile", "tofile", "memmap", "ctypeslib", "DataSource"}

SEED_OPERATOR = '''\
import numpy as np


def mutate(graph, labels, rng):
    """Move one random node into the community of a random neighbour."""
    v = int(rng.integers(graph.n))
    nbrs = graph.neighbours(v)
    if len(nbrs) > 0:
        labels[v] = labels[nbrs[rng.integers(len(nbrs))]]
    return labels
'''

SYSTEM_PROMPT = """\
You design mutation operators for ERBGA, a genetic algorithm for community
detection in graphs. Your operator runs alongside built-in operators (node move,
community split, community merge, bit flips) and is scored by the best
objective value the GA reaches on benchmark graphs with a small budget.

Write exactly one Python function:

    def mutate(graph, labels, rng):
        ...
        return labels

- labels: numpy int array, labels[i] = community of node i (a copy; you may edit it)
- rng: numpy.random.Generator; use it for all randomness
- graph: read-only, with
    graph.n, graph.m           node and edge counts
    graph.src, graph.dst       int arrays, edge e joins src[e] and dst[e]
    graph.degree               int array of node degrees
    graph.neighbours(v)        int array of v's neighbours
    graph.incident_edges(v)    int array of v's edge indices
- Return an int array of length graph.n. Community ids need not be contiguous.
  A community that ends up disconnected is split into its connected pieces.
- Only import: numpy, math, random, collections, itertools, heapq, functools,
  operator. No file, network or OS access.
- Must be fast: it runs thousands of times per evaluation. Prefer O(degree) or
  O(community size) work per call; avoid recomputing the objective.

The objective is {objective}. Good operators usually make targeted, local changes
that the objective is likely to reward, while still allowing big jumps sometimes.

Reply with a short rationale (2-3 sentences), then the code in one ```python block.
"""


@dataclass
class EvolveConfig:
    iterations: int = 10
    datasets: tuple[str, ...] = ("karate", "dolphins", "polbooks")
    seeds: int = 3
    objective: str = "modularity"
    pop_size: int = 30
    generations: int = 30
    operator_weight: float = 0.5
    timeout: float = 120.0
    memory_mb: int = 2048
    model: str = DEFAULT_MODEL
    effort: str = "high"
    top_k: int = 3
    out_dir: str = "evolve_runs"


@dataclass
class Candidate:
    id: int
    code: str
    score: float | None = None
    per_dataset: dict = field(default_factory=dict)
    error: str | None = None
    rationale: str = ""
    seconds: float = 0.0

    @property
    def valid(self) -> bool:
        return self.score is not None


# ---- static checks ----------------------------------------------------

def check_candidate(code: str) -> str | None:
    """Return a reason the code is rejected, or None if it may be run."""
    try:
        tree = ast.parse(code)
    except SyntaxError as exc:
        return f"syntax error: {exc}"
    if not any(isinstance(n, ast.FunctionDef) and n.name == "mutate" for n in tree.body):
        return "no top-level function named mutate"
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            mods = [a.name.split(".")[0] for a in node.names]
        elif isinstance(node, ast.ImportFrom):
            mods = [(node.module or "").split(".")[0]]
        else:
            mods = []
        bad = [m for m in mods if m not in ALLOWED_IMPORTS]
        if bad:
            return f"import of {bad[0]!r} is not allowed"
        if isinstance(node, ast.Name) and node.id in FORBIDDEN_NAMES:
            return f"use of {node.id!r} is not allowed"
        if isinstance(node, ast.Attribute) and (node.attr.startswith("__")
                                                or node.attr in FORBIDDEN_ATTRS):
            return f"attribute {node.attr!r} is not allowed"
    return None


def extract_code(text: str) -> tuple[str, str]:
    """Split a model reply into (rationale, code)."""
    match = re.search(r"```(?:python)?\s*\n(.*?)```", text, re.DOTALL)
    if not match:
        return "", text.strip()
    return text[:match.start()].strip(), match.group(1).strip() + "\n"


# ---- evaluation in a subprocess ---------------------------------------

_BOOT = ("import sys; sys.path.insert(0, sys.argv[1]); "
         "from erbga.evolve import _worker_main; _worker_main()")


def score_candidate(code: str | None, cfg: EvolveConfig) -> dict:
    """Score ``code`` (or the baseline, if None) in a child process.

    Returns {"score": float, "per_dataset": {...}} or {"error": str}.
    """
    request = json.dumps({"code": code, "config": asdict(cfg)})
    package_parent = str(Path(__file__).resolve().parent.parent)
    try:
        proc = subprocess.run(
            [sys.executable, "-I", "-c", _BOOT, package_parent],
            input=request, capture_output=True, text=True, timeout=cfg.timeout,
            preexec_fn=_limits(cfg) if os.name == "posix" else None,
        )
    except subprocess.TimeoutExpired:
        return {"error": f"timed out after {cfg.timeout:.0f}s"}
    lines = proc.stdout.strip().splitlines()
    if proc.returncode != 0 or not lines:
        tail = (proc.stderr.strip().splitlines() or ["(no output)"])[-1]
        return {"error": f"crashed: {tail[:500]}"}
    return json.loads(lines[-1])


def _limits(cfg: EvolveConfig):
    def apply():
        import resource
        mem = cfg.memory_mb * 1024 * 1024
        resource.setrlimit(resource.RLIMIT_AS, (mem, mem))
        cpu = int(cfg.timeout) + 5
        resource.setrlimit(resource.RLIMIT_CPU, (cpu, cpu))
    return apply


def _worker_main() -> None:
    """Child process: run the GA benchmark for one candidate, print JSON."""
    import numpy as np

    from .api import load_graph
    from .ga import ERBGA, GAConfig

    request = json.loads(sys.stdin.read())
    cfg = EvolveConfig(**{**request["config"], "datasets": tuple(request["config"]["datasets"])})
    operators = []
    if request["code"] is not None:
        namespace: dict = {}
        try:
            exec(compile(request["code"], "<candidate>", "exec"), namespace)
        except Exception as exc:
            print(json.dumps({"error": f"failed to load: {type(exc).__name__}: {exc}"}))
            return
        operators = [(cfg.operator_weight, namespace["mutate"])]

    per_dataset = {}
    try:
        for name in cfg.datasets:
            graph = load_graph(name)
            runs = []
            for seed in range(cfg.seeds):
                ga_cfg = GAConfig(pop_size=cfg.pop_size, generations=cfg.generations,
                                  local_search=False, patience=None, seed=seed)
                runs.append(ERBGA(graph, cfg.objective, ga_cfg, operators).run().fitness)
            per_dataset[name] = float(np.mean(runs))
    except Exception as exc:
        print(json.dumps({"error": f"failed while running: {type(exc).__name__}: {exc}"}))
        return
    print(json.dumps({"score": float(np.mean(list(per_dataset.values()))),
                      "per_dataset": per_dataset}))


# ---- the LLM ------------------------------------------------------------

Proposer = Callable[[str, str], str]  # (system, prompt) -> reply text


class ClaudeProposer:
    """Asks Claude for a new operator. Credentials come from the environment
    (ANTHROPIC_API_KEY or an ``ant auth login`` profile)."""

    def __init__(self, model: str = DEFAULT_MODEL, effort: str = "high"):
        import anthropic

        self.client = anthropic.Anthropic()
        self.model, self.effort = model, effort

    def __call__(self, system: str, prompt: str) -> str:
        import anthropic

        from .api import ERBGAError

        try:
            return self._ask(system, prompt)
        except (anthropic.AuthenticationError, TypeError) as exc:
            # TypeError: the SDK found no credentials at all
            raise ERBGAError("LLM_UNAVAILABLE", f"cannot call Claude: {exc}",
                             "set ANTHROPIC_API_KEY or run `ant auth login`") from exc
        except anthropic.APIError as exc:
            raise ERBGAError("LLM_ERROR", f"Claude API error: {exc}",
                             "check the model id and your account limits, then retry") from exc

    def _ask(self, system: str, prompt: str) -> str:
        response = self.client.beta.messages.create(
            model=self.model,
            max_tokens=16000,
            output_config={"effort": self.effort},
            # if a safety classifier declines, the API retries on a fallback model
            betas=["server-side-fallback-2026-07-01"],
            fallbacks="default",
            system=system,
            messages=[{"role": "user", "content": prompt}],
        )
        if response.stop_reason == "refusal":
            return ""
        return "".join(block.text for block in response.content if block.type == "text")


def build_prompt(baseline: float, best: list[Candidate], failures: list[Candidate]) -> str:
    parts = [f"Baseline score with built-in operators only: {baseline:.6f}. Higher is better.",
             "", "Best operators so far (worst first):"]
    for c in sorted(best, key=lambda c: c.score):
        parts += ["", f"### score {c.score:.6f}  "
                      f"({', '.join(f'{k} {v:.4f}' for k, v in c.per_dataset.items())})",
                  "```python", c.code.rstrip(), "```"]
    if failures:
        parts += ["", "Recent candidates that failed (avoid these mistakes):"]
        parts += [f"- {c.error}" for c in failures]
    parts += ["", "Write a new mutate operator that scores higher than all of the above."]
    return "\n".join(parts)


# ---- main loop ---------------------------------------------------------

def run_evolution(cfg: EvolveConfig, proposer: Proposer | None = None,
                  log: TextIO | None = sys.stdout) -> dict:
    """Run the evolution loop and return a JSON-safe summary.

    Writes ``log.jsonl`` (one line per candidate) and ``best_operator.py`` to a
    new timestamped directory under ``cfg.out_dir``.
    """
    proposer = proposer or ClaudeProposer(cfg.model, cfg.effort)
    run_dir = Path(cfg.out_dir) / time.strftime("%Y%m%d-%H%M%S")
    run_dir.mkdir(parents=True, exist_ok=True)
    log_path = run_dir / "log.jsonl"

    def say(msg: str) -> None:
        if log:
            print(msg, file=log, flush=True)

    baseline = score_candidate(None, cfg)
    if "error" in baseline:
        raise RuntimeError(f"baseline evaluation failed: {baseline['error']}")
    say(f"baseline {baseline['score']:.6f}  {baseline['per_dataset']}")

    candidates: list[Candidate] = []

    def evaluate(cand: Candidate) -> None:
        start = time.perf_counter()
        problem = check_candidate(cand.code)
        result = {"error": f"rejected: {problem}"} if problem else score_candidate(cand.code, cfg)
        cand.seconds = round(time.perf_counter() - start, 2)
        cand.score = result.get("score")
        cand.per_dataset = result.get("per_dataset", {})
        cand.error = result.get("error")
        candidates.append(cand)
        with log_path.open("a") as f:
            f.write(json.dumps(asdict(cand)) + "\n")
        status = f"score {cand.score:.6f}" if cand.valid else f"invalid: {cand.error}"
        say(f"candidate {cand.id}: {status} ({cand.seconds}s)")

    system = SYSTEM_PROMPT.format(objective=cfg.objective)
    evaluate(Candidate(id=0, code=SEED_OPERATOR, rationale="seed operator"))
    for i in range(1, cfg.iterations + 1):
        valid = sorted((c for c in candidates if c.valid), key=lambda c: -c.score)
        failures = [c for c in candidates[-3:] if not c.valid]
        reply = proposer(system, build_prompt(baseline["score"], valid[:cfg.top_k], failures))
        rationale, code = extract_code(reply)
        if not reply:
            evaluate(Candidate(id=i, code="", error="empty reply or refusal"))
            continue
        evaluate(Candidate(id=i, code=code, rationale=rationale))

    valid = sorted((c for c in candidates if c.valid), key=lambda c: -c.score)
    best = valid[0]
    best_path = run_dir / "best_operator.py"
    best_path.write_text(
        f"# ERBGA operator, score {best.score:.6f} vs baseline {baseline['score']:.6f}\n"
        f"# datasets {list(cfg.datasets)}, {cfg.seeds} seeds, objective {cfg.objective}\n"
        f"# use: erbga detect GRAPH --operator {best_path.name}\n\n{best.code}")
    return {
        "ok": True,
        "baseline": baseline,
        "best": {"id": best.id, "score": best.score, "per_dataset": best.per_dataset,
                 "rationale": best.rationale},
        "improvement": best.score - baseline["score"],
        "candidates": len(candidates),
        "valid_candidates": len(valid),
        "best_path": str(best_path),
        "log_path": str(log_path),
    }


def load_operator(path: str | Path):
    """Load ``mutate`` from an operator file (e.g. best_operator.py).

    The file is checked against the same allowlist as candidates, then run in
    this process: only load files you trust.
    """
    code = Path(path).read_text()
    problem = check_candidate(code)
    if problem:
        raise ValueError(f"{path}: {problem}")
    namespace: dict = {}
    exec(compile(code, str(path), "exec"), namespace)
    return namespace["mutate"]
