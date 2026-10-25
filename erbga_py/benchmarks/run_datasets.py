"""Run ERBGA over the bundled datasets with several seeds and print a Markdown table.

    python benchmarks/run_datasets.py --seeds 10
"""

import argparse
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from erbga import detect_communities, read_gml  # noqa: E402

# Proven maximum modularity (Aloise et al., 2010, column generation).
BEST_KNOWN = {"karate": 0.4198, "dolphins": 0.5285, "polbooks": 0.5272, "football": 0.6046}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--seeds", type=int, default=10)
    p.add_argument("--objective", default="modularity")
    args = p.parse_args()

    print("| graph | n | m | variant | best | mean ± sd | hit optimum | mean time (s) |")
    print("|---|---|---|---|---|---|---|---|")
    for name, best_known in BEST_KNOWN.items():
        g = read_gml(ROOT / "data" / f"{name}.gml")
        for variant, local_search in (("GA + local search", True), ("plain GA", False)):
            runs = [detect_communities(g, args.objective, seed=s, local_search=local_search)
                    for s in range(args.seeds)]
            fit = [r.fitness for r in runs]
            hits = sum(f >= best_known - 5e-5 for f in fit)
            sd = statistics.stdev(fit) if len(fit) > 1 else 0.0
            secs = statistics.mean(r.seconds for r in runs)
            hit = f"{hits}/{len(fit)}" if args.objective == "modularity" else "–"
            print(f"| {name} | {g.n} | {g.m} | {variant} | {max(fit):.4f} | "
                  f"{statistics.mean(fit):.4f} ± {sd:.4f} | {hit} | {secs:.2f} |", flush=True)


if __name__ == "__main__":
    main()
