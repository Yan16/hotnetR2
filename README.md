# hotnetR2

A planned standalone successor to hotnetR combining dataset preparation with
the v2_analysis2 and v3_analysis2 LDAK / hierarchical HotNet workflows.

**Status: design and planning revision; implementation and clean-machine
reproduction are not yet certified.** No existing analysis is modified.

Read in order:

1. [Workflow schematic](docs/workflow.md)
2. [Implementation plan](PLAN.md)
3. [Detailed task list](TASKS.md)
4. [API and coding contract](docs/api.md)
5. [Reproduction contract](docs/reproducibility.md)
6. [Architecture decisions and source findings](docs/decisions.md)

The references are `../v2_analysis2` and `../v3_analysis2` in the containing
project. Original hotnetR source is at `../../hotnetR`. These are development
references, never runtime dependencies of the finished package.

The first release must preserve names and arguments where their meaning is
unchanged, use explicit namespace-qualified tidyverse calls, provide Roxygen
documentation, and reproduce reference artifacts in isolated output trees.
It must install without hotnetR, analysis-local helper source trees, or
access to the original author's directories.

Public reference downloads, user-supplied GWAS inputs, software environments,
and output artifacts have separate versioned manifests. A new machine needs
the frozen inputs or verified source downloads, as well as the R package.

This repository will track reviewable implementation increments. Completed
work and pending acceptance gates are recorded in TASKS.md.

