import asyncio
import json

import pytest

pytest.importorskip("mcp")
from erbga.mcp_server import server  # noqa: E402


def call(name, args=None):
    result = asyncio.run(server.call_tool(name, args or {}))
    content = getattr(result, "content", result)  # mcp 2.x result object or 1.x list
    if isinstance(content, tuple):
        content = content[0]
    return json.loads(content[0].text)


def test_tools_listed():
    tools = asyncio.run(server.list_tools())
    assert {t.name for t in tools} == {"list_datasets", "list_objectives", "graph_summary",
                                       "detect_communities", "score_partition"}


def test_detect_via_mcp():
    out = call("detect_communities", {"graph": "karate", "seed": 0})
    assert out["ok"] and out["num_communities"] == 4


def test_inline_graph_and_score():
    out = call("score_partition", {"graph": {"edges": [[1, 2], [2, 3], [3, 1], [3, 4]]},
                                   "partition": {"1": 0, "2": 0, "3": 0, "4": 1}})
    assert out["connected"] and "modularity" in out["scores"]


def test_errors_are_returned_not_raised():
    out = call("graph_summary", {"graph": "missing.gml"})
    assert out["ok"] is False and out["error"]["code"] == "GRAPH_NOT_FOUND"


def test_list_datasets():
    names = [d["name"] for d in call("list_datasets")["datasets"]]
    assert names == ["karate", "dolphins", "polbooks", "football"]


def test_detect_with_islands_via_mcp():
    out = call("detect_communities", {"graph": "karate", "seed": 0, "islands": 4})
    assert out["ok"] and out["config"]["islands"] == 4
