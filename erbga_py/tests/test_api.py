import json

import pytest

from erbga import api
from erbga.cli import main


def test_load_graph_variants(tmp_path, karate):
    assert api.load_graph("karate").m == karate.m
    assert api.load_graph({"edges": [[1, 2], [2, 3]]}).n == 3
    assert api.load_graph('{"edges": [[1, 2]]}').m == 1
    node_link = {"nodes": [{"id": "a"}, {"id": "b"}], "links": [{"source": "a", "target": "b"}]}
    assert api.load_graph(node_link).m == 1
    edges = tmp_path / "g.txt"
    edges.write_text("# comment\n1 2\n2,3\n3\t1\n")
    assert api.load_graph(str(edges)).m == 3
    js = tmp_path / "g.json"
    js.write_text(json.dumps({"edges": [[1, 2]]}))
    assert api.load_graph(str(js)).name == "g"


@pytest.mark.parametrize("spec, code", [
    ("no-such-file.gml", "GRAPH_NOT_FOUND"),
    ({"nodes": [1, 2]}, "INVALID_GRAPH"),
    ({"edges": [[1, 1]]}, "EMPTY_GRAPH"),
    (42, "INVALID_GRAPH"),
])
def test_load_graph_errors(spec, code):
    with pytest.raises(api.ERBGAError) as info:
        api.load_graph(spec)
    assert info.value.code == code
    assert info.value.to_dict()["ok"] is False


def test_detect_output_is_json_and_compact():
    out = api.detect_communities("karate", seed=0, max_nodes_per_community=3)
    json.dumps(out)
    assert out["ok"] and out["num_communities"] == 4
    assert sum(c["size"] for c in out["communities"]) == 34
    assert all(len(c["nodes"]) <= 3 for c in out["communities"])
    assert out["communities"][0]["size"] >= out["communities"][-1]["size"]
    assert "history" not in out


@pytest.mark.parametrize("kwargs, code", [
    ({"objective": "nope"}, "UNKNOWN_OBJECTIVE"),
    ({"pop_size": 0}, "INVALID_CONFIG"),
    ({"popsize": 10}, "INVALID_CONFIG"),
])
def test_detect_errors(kwargs, code):
    with pytest.raises(api.ERBGAError) as info:
        api.detect_communities("karate", **kwargs)
    assert info.value.code == code


def test_score_partition_round_trip():
    found = api.detect_communities("karate", seed=0)
    partition = {str(n): c["id"] for c in found["communities"] for n in c["nodes"]}
    scored = api.score_partition("karate", partition)
    assert scored["scores"]["modularity"] == pytest.approx(found["fitness"])
    assert scored["connected"]


def test_score_partition_reports_disconnected():
    graph = {"edges": [[1, 2], [2, 3], [4, 5]]}
    out = api.score_partition(graph, [[1, 4], [2, 3, 5]], "modularity")
    assert not out["connected"]
    assert out["disconnected_communities"] == [0, 1]


@pytest.mark.parametrize("partition", [{"1": 0}, [[1, 2], [2, 3]], [[1, 2, 3, 99]], "x"])
def test_score_partition_errors(partition):
    with pytest.raises(api.ERBGAError) as info:
        api.score_partition({"edges": [[1, 2], [2, 3]]}, partition)
    assert info.value.code == "INVALID_PARTITION"


def test_summary_and_objectives():
    s = api.graph_summary("football")
    assert (s["nodes"], s["edges"], s["connected_components"]) == (115, 613, 1)
    names = [o["name"] for o in api.list_objectives()["objectives"]]
    assert names == ["modularity", "qs"]


def test_cli_json_mode(capsys):
    assert main(["detect", "karate", "--seed", "0", "--format", "json"]) == 0
    assert json.loads(capsys.readouterr().out)["num_communities"] == 4


def test_cli_json_error(capsys):
    assert main(["summary", "missing.gml", "--format", "json"]) == 1
    err = json.loads(capsys.readouterr().out)["error"]
    assert err["code"] == "GRAPH_NOT_FOUND" and err["hint"]


def test_cli_score_inline(capsys):
    assert main(["score", '{"edges": [[1,2],[2,3]]}', '[[1,2,3]]', "--format", "json"]) == 0
    assert json.loads(capsys.readouterr().out)["scores"]["modularity"] == pytest.approx(0.0)
    assert main(["score", "karate", "[" + ",".join(f"[{i}]" for i in range(1, 35)) + "]",
                 "--format", "json"]) == 0
    assert json.loads(capsys.readouterr().out)["num_communities"] == 34


def test_cli_without_subcommand_still_works(capsys):
    assert main(["karate", "--seed", "0", "--gen", "5"]) == 0
    assert "best modularity" in capsys.readouterr().out
