import json

import numpy as np
import pytest

from erbga import ERBGA, GAConfig, Graph
from erbga.evolve import (SEED_OPERATOR, EvolveConfig, build_prompt, check_candidate,
                          extract_code, load_operator, run_evolution, score_candidate)

FAST = dict(datasets=("karate",), seeds=1, pop_size=10, generations=5, timeout=60)


@pytest.mark.parametrize("code, reason", [
    ("def mutate(g, l, r) return l", "syntax"),
    ("def other(g, l, r): return l", "no top-level function"),
    ("import os\ndef mutate(g, l, r): return l", "'os'"),
    ("from subprocess import run\ndef mutate(g, l, r): return l", "'subprocess'"),
    ("def mutate(g, l, r):\n    open('x')\n    return l", "'open'"),
    ("def mutate(g, l, r):\n    return l.__class__", "'__class__'"),
    ("import numpy as np\ndef mutate(g, l, r):\n    np.save('x', l)\n    return l", "'save'"),
])
def test_check_candidate_rejects(code, reason):
    assert reason in check_candidate(code)


def test_seed_operator_passes_check():
    assert check_candidate(SEED_OPERATOR) is None


def test_extract_code():
    rationale, code = extract_code("Idea.\n```python\ndef mutate(g, l, r):\n    return l\n```\n")
    assert rationale == "Idea." and code.startswith("def mutate")


def test_score_candidate_and_failures():
    cfg = EvolveConfig(**FAST)
    assert "score" in score_candidate(SEED_OPERATOR, cfg)
    assert "shape" in score_candidate("def mutate(g, l, r):\n    return l[:2]\n", cfg)["error"]
    slow = EvolveConfig(**{**FAST, "timeout": 3})
    assert "timed out" in score_candidate("def mutate(g, l, r):\n    while True: pass\n",
                                          slow)["error"]


def test_operator_hook_is_used(karate):
    calls = []

    def op(graph, labels, rng):
        calls.append(1)
        return labels

    ERBGA(karate, config=GAConfig(seed=0, generations=3, pop_size=10),
          operators=[(1.0, op)]).run()
    assert calls


def test_operator_weights_validated(karate):
    with pytest.raises(ValueError):
        ERBGA(karate, operators=[(0.7, len), (0.6, len)])


def test_run_evolution_with_fake_llm(tmp_path):
    replies = iter([
        "Swap a node.\n```python\n" + SEED_OPERATOR + "```",
        "```python\nimport os\ndef mutate(g, l, r):\n    return l\n```",
        "",
    ])
    prompts = []

    def fake(system, prompt):
        prompts.append(prompt)
        return next(replies)

    cfg = EvolveConfig(**FAST, iterations=3, out_dir=str(tmp_path))
    out = run_evolution(cfg, proposer=fake, log=None)
    json.dumps(out)
    assert out["candidates"] == 4 and out["valid_candidates"] == 2
    assert "Baseline score" in prompts[0] and "def mutate" in prompts[0]
    assert "not allowed" in prompts[2]  # the failure is fed back
    log = [json.loads(line) for line in open(out["log_path"])]
    assert len(log) == 4
    mutate = load_operator(out["best_path"])
    path = Graph.from_edges([(0, 1), (1, 2)])
    assert len(mutate(path, np.array([0, 1, 2]), np.random.default_rng(0))) == 3


def test_build_prompt_orders_worst_first():
    from erbga.evolve import Candidate
    cands = [Candidate(1, "a", 0.5, {"k": 0.5}), Candidate(2, "b", 0.4, {"k": 0.4})]
    prompt = build_prompt(0.3, cands, [])
    assert prompt.index("0.400000") < prompt.index("0.500000")
