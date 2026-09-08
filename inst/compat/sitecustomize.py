"""Compatibility shim for the bundled legacy Hierarchical HotNet source.

The installed code calls ``networkx.connected_component_subgraphs``, removed
from modern NetworkX. Python automatically imports this module when this
directory is present in ``PYTHONPATH``. It restores the former generator API
without editing the installed Hierarchical HotNet checkout.
"""

import networkx as nx
import h5py
import numpy as np


if not hasattr(nx, "connected_component_subgraphs"):
    def connected_component_subgraphs(graph, copy=True):
        for component in nx.connected_components(graph):
            subgraph = graph.subgraph(component)
            yield subgraph.copy() if copy else subgraph


    nx.connected_component_subgraphs = connected_component_subgraphs


if not hasattr(h5py.Dataset, "value"):
    h5py.Dataset.value = property(lambda dataset: dataset[()])


if "int" not in np.__dict__:
    np.int = int
if "bool" not in np.__dict__:
    np.bool = bool
