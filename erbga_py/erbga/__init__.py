"""ERBGA: connected, objective-agnostic community detection with a genetic algorithm."""

from .ga import ERBGA, GAConfig, GAResult, detect_communities
from .graph import Graph, read_gml
from .objectives import OBJECTIVES, modularity, qs
from .partition import canonical_cut, connected_components, decode

__version__ = "1.0.0.dev0"

__all__ = [
    "ERBGA", "GAConfig", "GAResult", "detect_communities", "Graph", "read_gml",
    "OBJECTIVES", "modularity", "qs", "canonical_cut", "connected_components", "decode",
]
