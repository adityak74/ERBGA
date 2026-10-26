"""MCP server exposing ERBGA as tools for LLM agents.

Run it with ``erbga mcp`` (stdio). Example client config::

    {"mcpServers": {"erbga": {"command": "erbga", "args": ["mcp"]}}}

Tools return the same dicts as ``erbga.api``. Errors come back as
``{"ok": false, "error": {"code", "message", "hint"}}`` rather than as
exceptions, so the model can read the hint and retry.
"""

from __future__ import annotations

from typing import Any

from . import api

try:  # mcp >= 2
    from mcp.server.mcpserver import MCPServer as _Server
except ImportError:  # mcp 1.x
    from mcp.server.fastmcp import FastMCP as _Server
from mcp.types import ToolAnnotations

INSTRUCTIONS = """\
ERBGA finds communities (densely connected groups of nodes) in a graph with a
genetic algorithm. Every community it returns is connected.

Typical flow: graph_summary -> detect_communities -> (optionally) score_partition.
A graph argument is a bundled dataset name (see list_datasets), a file path
(.gml, .json, or an edge list with one 'u v' per line), or an inline object
{"edges": [[u, v], ...]}. Pass a seed for reproducible results. Results on large
graphs are truncated to max_nodes_per_community nodes per community; sizes stay exact.
"""

READ_ONLY = ToolAnnotations(readOnlyHint=True, idempotentHint=True, openWorldHint=False)

GraphArg = str | dict[str, Any]

server = _Server("erbga", instructions=INSTRUCTIONS)


def _safe(fn, *args, **kwargs) -> dict:
    try:
        return fn(*args, **kwargs)
    except api.ERBGAError as err:
        return err.to_dict()


@server.tool(annotations=READ_ONLY)
def list_datasets() -> dict:
    """List the graphs bundled with ERBGA, with their sizes."""
    rows = []
    for name in api.BUNDLED_DATASETS:
        s = api.graph_summary(name)
        rows.append({"name": name, "nodes": s["nodes"], "edges": s["edges"]})
    return {"ok": True, "datasets": rows}


@server.tool(annotations=READ_ONLY)
def list_objectives() -> dict:
    """List the objectives detect_communities can optimise. Higher is better for all."""
    return api.list_objectives()


@server.tool(annotations=READ_ONLY)
def graph_summary(graph: GraphArg) -> dict:
    """Nodes, edges, density, degree stats and connected components of a graph.

    Use this first to check a graph loads and to judge how long detection will take
    (roughly 1-5 s for a few hundred nodes with default settings).
    """
    return _safe(api.graph_summary, graph)


@server.tool(annotations=READ_ONLY)
def detect_communities(
    graph: GraphArg,
    objective: str = "modularity",
    seed: int | None = None,
    pop_size: int = 100,
    generations: int = 200,
    local_search: bool = True,
    patience: int | None = 100,
    islands: int = 1,
    refinement: bool = False,
    max_nodes_per_community: int | None = 200,
) -> dict:
    """Find communities in a graph with the ERBGA genetic algorithm.

    Returns a summary sentence, the objective value ("fitness"), communities
    sorted by size (each with id, size and node ids), and run metadata.
    Raise pop_size/generations for harder graphs; set local_search=false to
    trade speed for a more exploratory search. islands=4 splits the population
    into 4 islands with periodic migration, which helps when repeated runs
    stall at slightly different scores (premature convergence).
    refinement=true adds Leiden-style refinement to local search: slightly
    better on large graphs (1000+ nodes), occasionally worse on small ones.
    """
    return _safe(api.detect_communities, graph, objective,
                 max_nodes_per_community=max_nodes_per_community, seed=seed,
                 pop_size=pop_size, generations=generations,
                 local_search=local_search, patience=patience, islands=islands,
                 refinement=refinement)


@server.tool(annotations=READ_ONLY)
def score_partition(
    graph: GraphArg,
    partition: dict[str, int | str] | list[list[int | str]],
    objectives: list[str] | None = None,
) -> dict:
    """Score an existing partition of a graph without running the GA.

    partition is {node_id: community} or a list of node-id lists; every node
    must appear exactly once. Also reports whether each community is connected.
    Useful for comparing ERBGA's result with another method's.
    """
    return _safe(api.score_partition, graph, partition, objectives)


def main(transport: str = "stdio") -> None:
    server.run(transport)


if __name__ == "__main__":
    main()
