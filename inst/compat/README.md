# Local Hierarchical HotNet compatibility shim

The installed Hierarchical HotNet source is retained unchanged at
`/Users/li26191/software/hhotnet/hierarchical-hotnet/`. Its legacy NetworkX
call `connected_component_subgraphs()` was removed in NetworkX 3, and its
`h5py.Dataset.value` use was removed in h5py 3.
The legacy clustering code also uses NumPy's removed `np.int` and `np.bool`
aliases.

Set `PYTHONPATH` to this directory for local execution. Python imports
`sitecustomize.py` automatically and supplies the removed generator API using
the modern `networkx.connected_components()` implementation and the h5py
dataset-indexing API.

This shim is analysis-local and does not alter the installed program or the
system Python packages.
