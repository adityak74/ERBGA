"""ERBGA command line.

    erbga detect karate --seed 0                 # find communities
    erbga detect graph.txt --format json         # machine-readable output
    erbga score karate partition.json            # score a given partition
    erbga summary football                       # graph statistics
    erbga objectives                             # list objectives
    erbga mcp                                    # run the MCP server (stdio)
    erbga evolve --iterations 10                 # LLM-evolved mutation operators

``erbga GRAPH ...`` without a subcommand means ``erbga detect GRAPH ...``.
With ``--format json`` the result, or an error object, is the only thing
written to stdout; progress goes to stderr. Exit codes: 0 ok, 1 error.
"""

from __future__ import annotations

import argparse
import json
import sys

from . import __version__, api
from .ga import GAConfig
from .objectives import OBJECTIVES

SUBCOMMANDS = ("detect", "score", "summary", "objectives", "mcp", "evolve")


def build_parser() -> argparse.ArgumentParser:
    d = GAConfig()
    fmt = argparse.ArgumentParser(add_help=False)
    fmt.add_argument("--format", choices=["text", "json"], default="text",
                     help="output format (json: one JSON object on stdout)")

    p = argparse.ArgumentParser(prog="erbga", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    det = sub.add_parser("detect", parents=[fmt], help="find communities")
    det.add_argument("graph", help="dataset name, .gml/.json/edge-list path")
    det.add_argument("-o", "--objective", default="modularity", choices=sorted(OBJECTIVES))
    det.add_argument("--pop", type=int, default=d.pop_size, help="population size")
    det.add_argument("--gen", type=int, default=d.generations, help="max generations")
    det.add_argument("--elite", type=float, default=d.elite_frac, help="elite fraction")
    det.add_argument("--tournament", type=int, default=d.tournament_size)
    det.add_argument("--crossover", default=d.crossover, choices=["uniform", "one_point"])
    det.add_argument("--crossover-rate", type=float, default=d.crossover_rate)
    det.add_argument("--mutation-rate", type=float, default=d.mutation_rate)
    det.add_argument("--init", default=d.init, choices=["random", "degree"])
    det.add_argument("--no-local-search", action="store_true", help="plain GA, no refinement")
    det.add_argument("--patience", type=int, default=d.patience,
                     help="stop after N generations without improvement (0 = never)")
    det.add_argument("--seed", type=int, default=None)
    det.add_argument("--max-nodes", type=int, default=None, metavar="N",
                     help="list at most N nodes per community in the output")
    det.add_argument("--history", action="store_true", help="include per-generation history")
    det.add_argument("--operator", metavar="FILE",
                     help="extra mutation operator file from 'erbga evolve' (runs its code)")
    det.add_argument("--operator-weight", type=float, default=0.3,
                     help="probability of using --operator per mutation")
    det.add_argument("--json", metavar="PATH", help="also write the full result as JSON")
    det.add_argument("-v", "--verbose", action="store_true",
                     help="print every generation (to stderr with --format json)")

    sc = sub.add_parser("score", parents=[fmt], help="score a given partition")
    sc.add_argument("graph")
    sc.add_argument("partition",
                    help="JSON file or inline JSON: {node: community} or [[nodes], ...]")
    sc.add_argument("-o", "--objective", action="append", choices=sorted(OBJECTIVES),
                    help="objective(s) to compute (default: all)")

    sm = sub.add_parser("summary", parents=[fmt], help="graph statistics")
    sm.add_argument("graph")

    sub.add_parser("objectives", parents=[fmt], help="list objectives")

    mc = sub.add_parser("mcp", help="run the MCP server")
    mc.add_argument("--transport", default="stdio", choices=["stdio", "sse", "streamable-http"])

    ev = sub.add_parser("evolve", parents=[fmt],
                        help="evolve mutation operators with an LLM (runs LLM-written code)")
    ev.add_argument("--iterations", type=int, default=10)
    ev.add_argument("-o", "--objective", default="modularity", choices=sorted(OBJECTIVES))
    ev.add_argument("--model", default=None, help="Claude model id (default: see evolve.py)")
    ev.add_argument("--datasets", nargs="+", default=["karate", "dolphins", "polbooks"])
    ev.add_argument("--seeds", type=int, default=3, help="GA runs per dataset per candidate")
    ev.add_argument("--timeout", type=float, default=120, help="seconds per candidate")
    ev.add_argument("--out", default="evolve_runs", help="directory for the log and best operator")
    return p


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] not in SUBCOMMANDS and argv[0] not in ("-h", "--help", "--version"):
        argv.insert(0, "detect")  # backwards compatible: erbga GRAPH ...
    args = build_parser().parse_args(argv)
    as_json = getattr(args, "format", "text") == "json"
    try:
        return COMMANDS[args.command](args, as_json)
    except api.ERBGAError as err:
        if as_json:
            print(json.dumps(err.to_dict()))
        else:
            print(f"error [{err.code}]: {err.message}", file=sys.stderr)
            if err.hint:
                print(f"hint: {err.hint}", file=sys.stderr)
        return 1


def _emit(data: dict, as_json: bool, text: str) -> int:
    print(json.dumps(data, indent=2) if as_json else text)
    return 0


def cmd_detect(args, as_json: bool) -> int:
    graph = api.load_graph(args.graph)
    if not as_json:
        print(f"{graph.name}: {graph.n} nodes, {graph.m} edges, objective={args.objective}")
    config = dict(
        pop_size=args.pop, generations=args.gen, elite_frac=args.elite,
        tournament_size=args.tournament, crossover=args.crossover,
        crossover_rate=args.crossover_rate, mutation_rate=args.mutation_rate,
        init=args.init, local_search=not args.no_local_search,
        patience=args.patience or None, seed=args.seed,
    )
    stream = sys.stderr if as_json else sys.stdout

    def progress(gen, stats):
        print(f"  gen {gen:4d}  best {stats[0]:.6f}  mean {stats[1]:.6f}", file=stream)

    operators = ()
    if args.operator:
        from .evolve import load_operator
        try:
            operators = [(args.operator_weight, load_operator(args.operator))]
        except (OSError, ValueError) as exc:
            raise api.ERBGAError("INVALID_CONFIG", f"cannot load operator: {exc}") from exc
    result = api.detect_communities(graph, args.objective,
                                    max_nodes_per_community=args.max_nodes,
                                    include_history=args.history or bool(args.json),
                                    progress=progress if args.verbose else None,
                                    operators=operators, **config)
    if args.json:
        with open(args.json, "w") as f:
            json.dump(result, f, indent=2)
    run = result["run"]
    text = (f"best {args.objective} = {result['fitness']:.6f} with "
            f"{result['num_communities']} communities ({run['generations']} generations, "
            f"{run['evaluations']} evaluations, {run['seconds']:.2f}s)")
    if args.json:
        text += f"\nwrote {args.json}"
    return _emit(result, as_json, text)


def cmd_score(args, as_json: bool) -> int:
    partition = _read_json_arg(args.partition)
    result = api.score_partition(args.graph, partition, args.objective)
    lines = [f"{result['num_communities']} communities, "
             f"{'all connected' if result['connected'] else 'some disconnected'}"]
    lines += [f"{name} = {value:.6f}" for name, value in result["scores"].items()]
    return _emit(result, as_json, "\n".join(lines))


def cmd_summary(args, as_json: bool) -> int:
    s = api.graph_summary(args.graph)
    text = (f"{s['graph']}: {s['nodes']} nodes, {s['edges']} edges, density {s['density']}\n"
            f"degree min/mean/max {s['degree']['min']}/{s['degree']['mean']}/"
            f"{s['degree']['max']}, {s['connected_components']} connected component(s), "
            f"{s['isolated_nodes']} isolated")
    return _emit(s, as_json, text)


def cmd_objectives(args, as_json: bool) -> int:
    result = api.list_objectives()
    text = "\n".join(f"{o['name']}: {o['description']}" for o in result["objectives"])
    return _emit(result, as_json, text)


def cmd_mcp(args, as_json: bool) -> int:
    from .mcp_server import main as mcp_main
    mcp_main(args.transport)
    return 0


def cmd_evolve(args, as_json: bool) -> int:
    from .evolve import EvolveConfig, run_evolution
    cfg = EvolveConfig(iterations=args.iterations, objective=args.objective,
                       datasets=tuple(args.datasets),
                       seeds=args.seeds, timeout=args.timeout, out_dir=args.out)
    if args.model:
        cfg.model = args.model
    log = sys.stderr if as_json else sys.stdout
    result = run_evolution(cfg, log=log)
    text = (f"baseline score {result['baseline']['score']:.6f}, "
            f"best score {result['best']['score']:.6f} "
            f"({result['improvement']:+.6f}); best operator saved to {result['best_path']}")
    return _emit(result, as_json, text)


def _read_json_arg(value: str):
    stripped = value.strip()
    try:
        if stripped.startswith(("{", "[")):
            return json.loads(stripped)
        with open(value) as f:
            return json.load(f)
    except FileNotFoundError:
        raise api.ERBGAError("INVALID_PARTITION", f"no such partition file: {value!r}",
                             "pass a JSON file path or inline JSON") from None
    except json.JSONDecodeError as exc:
        raise api.ERBGAError("INVALID_PARTITION", f"partition is not valid JSON: {exc}") from None


COMMANDS = {
    "detect": cmd_detect, "score": cmd_score, "summary": cmd_summary,
    "objectives": cmd_objectives, "mcp": cmd_mcp, "evolve": cmd_evolve,
}


if __name__ == "__main__":
    sys.exit(main())
