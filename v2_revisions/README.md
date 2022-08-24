# Assembling Overlapping Clusters Implementation

## Getting Started

We use Python3 (tested with `>= 3.7`) with dependencies listed in `requirements.txt`. The entry point is at `aoc.py`. To install dependencies:

```bash
python3 -m pip install -r requirements.txt
```

Then, given input an existing edge-list encoding a graph
and an existing clustering (`cluster_id node_id` pairs, or the output format of [IKC.py](https://github.com/chackoge/ERNIE_Plus/tree/master/Illinois/clustering/eleanor/code)),
AOC_m can be run as follows:

```bash
python3 aoc.py -g [graph_edgelist_file] -c [input_disjoint_clusters] --inclusion-criterion mcd
```

which will print the resulting clustering to stdout and use all non-singleton nodes in the input clustering as the candidates for augmenting the input disjoint clusters.

For running AOC_k with some predetermined $k$, the following command can be run:

```
python3 aoc.py -g [graph_edgelist_file] -c [input_disjoint_clusters] --inclusion-criterion k --min-k-core [k]
```

If the clusters are supposed to be written to a file. The output path can be specified by `-o output_clustering_path` which will write out the clustering and also some additional statistics.

## Other useful flags

**Use user-specified candidates**: the default candidate list is determined by all nodes that appear in the input clustering. If the following flag is set, the candidate list is determined by the user-specified list.

```bash
--candidate-file [list_of_candidate_nodes]
```

**Change output statistics directory**: the default behavior, if an output path of the clustering is set, is to write out the clustering and some additional statistics to a directory `{output_path}.stats` where `output_path` is the path of the output clustering. This flag can be used to change the output directory *for the statistic files*.

```bash
--stats-directory [output_directory]
```

For the full list of flags, run the script with ```--help```.

## Using with IKC

AOC augments an existing disjoint clustering. In our study, the disjoint clustering is produced by [IKC](https://direct.mit.edu/qss/article/3/1/289/109629/Center-periphery-structure-in-research-communities). We provide a full-stack example here of how to use AOC_m with [the IKC script](https://github.com/chackoge/ERNIE_Plus/tree/master/Illinois/clustering/eleanor/code).

```bash
# note that IKC.py has different dependencies; install them first
python3 IKC.py -e [graph_edgelist_file] -k [k] -o [ikc_clustering_path]
python3 aoc.py -g [graph_edgelist_file] -c [ikc_clustering_path] --inclusion-criterion mcd
```