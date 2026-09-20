"""Cut an observed hierarchy using the installed HHN implementation, no null model."""
import argparse
import sys

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--source", required=True)
parser.add_argument("--edges", required=True)
parser.add_argument("--index", required=True)
parser.add_argument("--delta", required=True, type=float)
parser.add_argument("--output", required=True)
args = parser.parse_args()
sys.path.insert(0, args.source)
from process_hierarchies import cut_hierarchy

clusters = sorted((sorted(c) for c in cut_hierarchy(args.edges, args.index, args.delta)),
                  key=lambda c: (-len(c), c))
with open(args.output, "w") as out:
    out.write("# Method: manual_delta_no_permutations\n")
    out.write("# Observed cut height: {}\n".format(args.delta))
    out.write("# Significance: not evaluated\n# Clusters:\n")
    for cluster in clusters:
        out.write("\t".join(cluster) + "\n")
