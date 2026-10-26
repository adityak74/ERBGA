from pathlib import Path

import pytest

from erbga import read_gml

DATA = Path(__file__).resolve().parent.parent / "erbga" / "data"


@pytest.fixture(scope="session")
def karate():
    return read_gml(DATA / "karate.gml")


@pytest.fixture(scope="session")
def football():
    return read_gml(DATA / "football.gml")
