"""Command line entry point: ``python -m erbga data/karate.gml``."""

from __future__ import annotations

import argparse
import json
import sys

from . import __version__
from .ga import ERBGA, GAConfig
from .graph import read_gml
from .objectives import OBJECTIVES


def build_parser() -> argparse.ArgumentParser:
    d = GAConfig()
    p = argparse.ArgumentParser(prog="erbga", description=__doc__)
    p.add_argument("input", help="graph in GML format")
    p.add_argument("-o", "--objective", default="modularity", choices=sorted(OBJECTIVES))
    p.add_argument("--pop", type=int, default=d.pop_size, help="population size")
    p.add_argument("--gen", type=int, default=d.generations, help="max generations")
    p.add_argument("--elite", type=float, default=d.elite_frac, help="elite fraction")
    p.add_argument("--tournament", type=int, default=d.tournament_size)
    p.add_argument("--crossover", default=d.crossover, choices=["uniform", "one_point"])
    p.add_argument("--crossover-rate", type=float, default=d.crossover_rate)
    p.add_argument("--mutation-rate", type=float, default=d.mutation_rate)
    p.add_argument("--init", default=d.init, choices=["random", "degree"])
    p.add_argument("--no-local-search", action="store_true", help="plain GA, no refinement")
    p.add_argument("--patience", type=int, default=d.patience,
                   help="stop after N generations without improvement (0 = never)")
    p.add_argument("--seed", type=int, default=None)
    p.add_argument("--json", metavar="PATH", help="write the result as JSON")
    p.add_argument("-v", "--verbose", action="store_true", help="print every generation")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    graph = read_gml(args.input)
    config = GAConfig(
        pop_size=args.pop, generations=args.gen, elite_frac=args.elite,
        tournament_size=args.tournament, crossover=args.crossover,
        crossover_rate=args.crossover_rate, mutation_rate=args.mutation_rate,
        init=args.init, local_search=not args.no_local_search,
        patience=args.patience or None, seed=args.seed,
    )
    print(f"{graph.name}: {graph.n} nodes, {graph.m} edges, objective={args.objective}")

    def progress(gen, stats):
        print(f"  gen {gen:4d}  best {stats[0]:.6f}  mean {stats[1]:.6f}")

    result = ERBGA(graph, args.objective, config).run(progress if args.verbose else None)
    print(f"best {args.objective} = {result.fitness:.6f} with {result.num_communities} communities "
          f"({result.generations} generations, {result.evaluations} evaluations, "
          f"{result.seconds:.2f}s)")

    if args.json:
        with open(args.json, "w") as f:
            json.dump({
                "graph": graph.name, "objective": args.objective, "seed": args.seed,
                "fitness": result.fitness, "num_communities": result.num_communities,
                "communities": {str(k): v for k, v in result.communities.items()},
                "history": result.history, "generations": result.generations,
                "evaluations": result.evaluations, "seconds": result.seconds,
            }, f, indent=2)
        print(f"wrote {args.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
