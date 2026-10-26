"""Time ERBGA on planted-partition graphs of growing size.

    python benchmarks/scaling.py --sizes 500 1000 2000

Each graph has 10 equal blocks; edge probabilities are scaled so the average
degree stays near 7. Prints the planted partition's modularity for reference.
"""

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from erbga import Graph, detect_communities, modularity  # noqa: E402


def planted_partition(n: int, blocks: int, p_in: float, p_out: float, seed: int = 0):
    rng = np.random.default_rng(seed)
    block = np.repeat(np.arange(blocks), n // blocks)
    iu, ju = np.triu_indices(len(block), 1)
    same = block[iu] == block[ju]
    keep = rng.random(len(same)) < np.where(same, p_in, p_out)
    edges = zip(iu[keep].tolist(), ju[keep].tolist(), strict=True)
    return Graph.from_edges(edges, nodes=range(len(block)), name=f"planted{n}"), block


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sizes", type=int, nargs="+", default=[500, 1000, 2000])
    p.add_argument("--seeds", type=int, default=1)
    args = p.parse_args()

    print("| n | m | planted Q | best Q | communities | generations | seconds |")
    print("|---|---|---|---|---|---|---|")
    for n in args.sizes:
        g, block = planted_partition(n, 10, 50 / n, 2 / n)
        for seed in range(args.seeds):
            r = detect_communities(g, seed=seed)
            print(f"| {g.n} | {g.m} | {modularity(g, block):.4f} | {r.fitness:.4f} | "
                  f"{r.num_communities} | {r.generations} | {r.seconds:.1f} |", flush=True)


if __name__ == "__main__":
    main()
