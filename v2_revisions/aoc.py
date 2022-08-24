"""Assembling Overlapping Cluster main executable"""
import logging
import random
import os
import statistics
import sys
import time
from collections import Counter, defaultdict
from copy import deepcopy
import click
import networkx as nx
import numpy as np
import overlapping_clusters_stats as ocs


@click.command()
@click.option(
    "-c",
    "--clustering",
    required=True,
    type=click.Path(exists=True),
    help="Clustering as input for overlapping step; each line in the clustering file contains a whitespace-separated cluster_id node_id pair",
)
@click.option(
    "-g",
    "--network-file",
    required=True,
    type=click.Path(exists=True),
    help="Edgelist file (each line contains two nodes separated by whitespace), the background network of the input clustering",
)
@click.option(
    "-o",
    "--output-path",
    required=False,
    type=click.Path(),
    help="File path to save the output clustering; the directory of this path must already exist",
)
@click.option(
    "--min-k-core",
    required=False,
    type=int,
    help="Minimum k-value for candidate addition; only used if --inclusion-criterion is k",
)
@click.option(
    "--rank-type",
    required=False,
    type=click.Choice(["percent", "percentile"]),
    help="Ranking metric type for candidate consideration step; for advanced usage of non-standard candidate generation",
)
@click.option(
    "--rank-val",
    required=False,
    type=int,
    help="Rank value for rank type to set as threshold for candidate consideration; for advanced usage of non-standard candidate generation",
)
@click.option(
    "--inclusion-criterion",
    required=True,
    type=click.Choice(["k", "mcd"]),
    help="Criterion to include candidate nodes to a cluster; 'k' specifies AOC_k; 'mcd' specifies AOC_m",
)
@click.option(
    "--candidate-criterion",
    required=False,
    type=click.Choice(["total_degree", "indegree", "random"]),
    help="Criterion to generate candidates for the overlapping cluster step; if not specified, defaults to first --candidate-file and then the nodes appearing in --clustering",
)
@click.option(
    "--candidate-file",
    required=False,
    type=click.Path(exists=True),
    help="A newline separated list of custom candidate nodes to consider for augmenting the input clustering; has higher priority than --candidate-criterion",
)
@click.option(
    "--stats-directory",
    required=False,
    type=click.Path(),
    help="Directory to save statistics of the input and output clustering; if not specified, defaults to --output-path suffixed with .stats if --output-path is specified",
)
def main(
    clustering,
    network_file,
    output_path,
    min_k_core,
    rank_type,
    rank_val,
    inclusion_criterion,
    candidate_criterion,
    candidate_file,
    stats_directory,
):
    """Constructs an overlapping set of clusters from an input disjoint cluster."""
    logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.INFO)
    display_cluster_stats, save_output_stats = True, True

    output_prefix = (
        os.path.join(f"{output_path}.stats") if not stats_directory else stats_directory
    )

    if not output_path:
        # if the user does not specify an output path, we pretend the output should go to STDOUT.
        # In addition, we forgo outputting statistics since it could be that the user
        # is invoking AOC in a programmtic pipeline (reading its STDOUT) and only wants the clustering.
        save_output_stats = False

    if not os.path.isdir(output_prefix):
        if save_output_stats:
            os.makedirs(output_prefix)
            logging.info(f"Created output directory: {output_prefix}")

    if inclusion_criterion == "k":
        if min_k_core is None:
            raise ValueError(
                "Minimum k-core value must be provided for k-core criterion"
            )
        elif min_k_core < 1:
            raise ValueError("Minimum k-core value must be greater than 0")

    # Parse Network Data
    node_info = network_to_dict(network_file)
    G = nx.read_edgelist(
        network_file
    )  # TODO: we really don't need to read the graph twice
    logging.info("Finished Network Parsing")

    # Parse Clustering Data
    clusters, node_to_cluster_id = clustering_to_dict(
        clustering,
        min_k_core,
    )
    logging.info("Finished Cluster Parsing")

    # Generate Original Cluster by Cluster Stats
    if save_output_stats:
        _ = ocs.cluster_analysis(
            G,
            os.path.join(output_prefix, "original_cluster_stats.csv"),
            clusters,
            node_info,
            is_overlapping=False,
        )
        logging.info(
            "Finished Generating Cluster by Cluster Stats for Original Clusters"
        )

    # Generate Candidates
    if candidate_file is not None:
        candidates = generate_candidate_file(candidate_file)
    elif candidate_criterion is not None:
        candidates = generate_candidates(
            node_info, rank_val, candidate_criterion, rank_type
        )
    else:
        logging.info(
            "Did not specify candidates, now using all non-singleton nodes in the original cluster"
        )
        non_singleton_node_set = [
            clusters["Full Clusters"][key] for key in clusters["Full Clusters"]
        ]
        candidates = set().union(*non_singleton_node_set)
    logging.info(
        "Generated "
        + str(len(candidates))
        + " candidate(s) based on user chosen criteria"
    )
    candidates = sorted(
        candidates, key=lambda x: node_info["total_degree"][x], reverse=True
    )
    logging.info("Finished Candidate Generation")

    # Run Iterative OC Generation
    oc_start = time.time()
    overlapping_clusters = deepcopy(clusters)

    (
        overlapping_clusters,
        _overlapping_node_to_cluster_id,
    ) = overlapping_clusters_construction(
        clusters, G, node_info, node_to_cluster_id, inclusion_criterion, candidates
    )
    oc_end = time.time()
    logging.info(
        "Overlapping Clustering Construction Elapsed Time (s): "
        + str(oc_end - oc_start)
    )
    logging.info("Finished Overlapping Clustering Generation")

    # Generic cluster stats
    if display_cluster_stats:
        cluster_stats(
            overlapping_clusters,
            node_info,
            G,
            output_prefix,
            min_k_core,
            inclusion_criterion,
            save_output_stats,
        )
        logging.info("Finished Cluster Stats")

    if save_output_stats:
        # Generate Overlapping Cluster by Cluster Stats
        _ = ocs.cluster_analysis(
            G,
            os.path.join(output_prefix, "overlapping_cluster_stats.csv"),
            overlapping_clusters,
            node_info,
            is_overlapping=True,
        )
        logging.info(
            "Finished Generating Cluster by Cluster Stats for Overlapping Clusters"
        )

        # Generate Overlapping Cluster Intersection Stats
        _ = ocs.cluster_intersection_analysis(
            os.path.join(output_prefix, "intersection_stats.csv"),
            overlapping_clusters,
        )
        logging.info(
            "Finished Generating Cluster Cluster Intersection Stats for Overlapping Clusters"
        )
    # Save Final OC Output to either file or STDOUT
    if not output_path:
        overlapping_clusters_to_output(sys.stdout, overlapping_clusters)
    else:
        with open(output_path, "w+") as f:
            overlapping_clusters_to_output(f, overlapping_clusters)
    logging.info("Finished Outputting Final Overlapping Clustering")


def cluster_stats(
    clusters,
    node_info,
    G,
    output_prefix,
    min_k_core,
    inclusion_criterion,
    save_outputs,
):
    """
    Computes statistics for a given clustering

    Input:
    network_file str - path to network tsv file
    clusters {dict} - dictionary of clusters by cluster id
    node_info {dict} -  dictionary of info on nodes in network
    G networkx.Graph - networkx representation of input network
    candidates [list] - list of candidate node ids
    experiment_num - experiment number to save results to
    min_k_core int - minimum k value to parse
    experiment_name str - experiment name to identify specific run of OC generation
    inclusion_criterion - k or mcd used in candidate placement
    save_outputs bool - bool value to decide if stats should be saved
    include_markers bool - bool value to indicate if marker node analysis can be run

    Output:
    None - Cluster stats logging.infoed to console
    Format: Num Clusters, Num Singleton Nodes, Min Cluster Size, Median Cluster Size, Max Cluster Size, Node Coverage, Edge Coverage, Candidate Node Coverage, Candidate Edge Coverage
    """
    (
        num_clusters,
        num_singletons,
        min_size,
        max_size,
        median_size,
        node_coverage,
    ) = basic_cluster_info(clusters, node_info)
    edge_coverage = get_edge_coverage(G, clusters)

    logging.info(f"Num Clusters: {num_clusters}")
    logging.info(f"Num Singleton Nodes: {num_singletons}")
    logging.info(f"Min Cluster Size: {min_size}")
    logging.info(f"Median Cluster Size: {median_size}")
    logging.info(f"Max Cluster Size: {max_size}")
    logging.info(f"Node Coverage: {node_coverage}")
    logging.info(f"Edge Coverage: {edge_coverage}")

    if save_outputs:
        f = open(os.path.join(output_prefix, "cluster_basics.csv"), "a")
        f.write(
            ",".join(
                [
                    str(s)
                    for s in [
                        min_k_core,
                        output_prefix,
                        inclusion_criterion,
                        num_clusters,
                        num_singletons,
                        min_size,
                        median_size,
                        max_size,
                        node_coverage,
                        edge_coverage,
                    ]
                ]
            )
            + "\n"
        )
        f.close()


def basic_cluster_info(clusters, node_info):
    """
    Computes basic statistics for a given clustering

    Input:
    clusters {dict} - dictionary of clusters by cluster id
    node_info {dict} -  dictionary of info on nodes in network

    Output:
    num_clusters int
    num_singletons int
    min_cluster_size int
    max_cluster_size int
    median_cluster_size float
    node_coverage float
    """
    num_nodes = len(node_info["total_degree"].keys())
    num_clusters = len(clusters["Full Clusters"].keys())

    non_singleton_node_set = [
        clusters["Full Clusters"][key] for key in clusters["Full Clusters"]
    ]
    non_singleton_node_count = len(set().union(*non_singleton_node_set))
    cluster_size_array = [
        len(clusters["Full Clusters"][key]) for key in clusters["Full Clusters"]
    ]

    num_singletons = num_nodes - non_singleton_node_count
    min_cluster_size = min(cluster_size_array)
    max_cluster_size = max(cluster_size_array)
    median_cluster_size = statistics.median(cluster_size_array)
    node_coverage = non_singleton_node_count / num_nodes

    return (
        num_clusters,
        num_singletons,
        min_cluster_size,
        max_cluster_size,
        median_cluster_size,
        node_coverage,
    )


def get_edge_coverage(G, clusters):
    """
    Computes edge coverage of a given clustering

    Input:
    G networkx.Graph - networkx representation of input network
    clusters {dict} - dictionary of clusters by cluster id

    Output:
    edge_coverage float
    """
    total_edges = G.number_of_edges()
    counted_edges = set()
    for c in clusters["Full Clusters"].keys():
        nodes = list(clusters["Full Clusters"][c])
        G_prime = G.subgraph(nodes)
        counted_edges.update(set(G_prime.edges))
    edge_coverage = len(counted_edges) / total_edges
    return edge_coverage


def get_coverage_node_list(G, clusters, candidates):
    """
    Computes node and edge coverage of a given clustering for a specific list of nodes in the network

    Input:
    G networkx.Graph - networkx representation of input network
    clusters {dict} - dictionary of clusters by cluster id
    candidates [list] list of node ids of the subset of nodes to check coverage

    Output:
    node_coverage float
    edge_coverage float
    """
    candidates_in_cluster = set()
    total_edges = len(set(G.edges(candidates)))
    counted_edges = set()
    for c in clusters["Full Clusters"].keys():
        nodes = list(clusters["Full Clusters"][c])
        candidates_in_cluster.update(set(nodes).intersection(set(candidates)))
        G_prime = G.subgraph(nodes)
        counted_edges.update(set(G_prime.edges(candidates)))
    node_coverage = len(candidates_in_cluster) / len(candidates)
    edge_coverage = len(counted_edges) / total_edges
    return node_coverage, edge_coverage


def generate_candidates(node_info, rank_val, candidate_criterion, rank_type):
    """
    Returns list of candidate node ids given top n percent of nodes in network by total degree

    Input:
    node_info {dict} - dictionary of info on nodes in network
    rank_val int - integer representing rank value for a rank type (percent, percentile) to consider as minimum threshold
    candidate_criterion str - criterion to judge acceptability of node as a candidate
    rank_type str - percent or percentile used as rank metric for candidate consideration

    Output:
    candidates [list] - list of candidate node ids
    """
    if candidate_criterion == "random":
        candidates = list(node_info["total_degree"].keys())
        if rank_type == "percentile":
            rank_val = 100 - rank_val
            candidates = random.sample(
                candidates,
                int((rank_val * len(node_info["total_degree"].keys())) / 100),
            )
    else:
        degree_map = Counter(node_info[candidate_criterion])
        degrees = list(map(int, degree_map.values()))
        if rank_type == "percent":
            degree_cutoff = np.percentile(degrees, 100 - rank_val)
            candidates = [k for k, v in degree_map.items() if v >= degree_cutoff]
        else:
            bdegree_cutoff = np.percentile(degrees, max(rank_val - 0.5, 0))
            tdegree_cutoff = np.percentile(degrees, min(rank_val + 0.5, 100))
            candidates = [
                k
                for k, v in degree_map.items()
                if bdegree_cutoff <= v <= tdegree_cutoff
            ]

    candidates.sort(reverse=False, key=lambda n: node_info["indegree"][n])
    return candidates


def generate_candidate_file(candidate_file):
    """
    Returns list of candidate node ids given by custom candidate node file with one node id on each line

    Input:
    candidate_file str - file path to candidate file

    Output:
    candidates [list] - list of candidate node ids
    """
    candidates = []
    candidate_file_reader = open(candidate_file, "r")
    line = candidate_file_reader.readline()

    while line != "":
        node_id = str(line.strip())
        candidates.append(node_id)
        line = candidate_file_reader.readline()

    return candidates


def overlapping_clusters_construction(
    clusters, G, node_info, node_to_cluster_id, inclusion_criterion, candidates
):
    """
    Overlapping clustering method that adds candidate nodes to clusters based on some inclusion criterion

    Input:
    clusters {dict} - dictionary of clusters by cluster id
    G networkx.Graph - networkx representation of input network
    node_info {dict} - dictionary of info on nodes in network
    node_to_cluster_id {dict} - dictionary of nodes mapped to the disjoint cluster id they are part of
    inclusion_criterion str - criteria to include a node into a cluster {k, mcd}
    candidates [list] - list of node ids of the subset of nodes to check coverage

    Output:
    overlapping_clusters {dict} - dictionary of overlapping clusters by cluster id
    overlapping_node_to_cluster_id {dict} - dictionary of nodes mapped to the overlapping cluster ids they are part of
    """
    overlapping_clusters = deepcopy(clusters)
    overlapping_node_to_cluster_id = deepcopy(node_to_cluster_id)

    l = G.size()
    for cluster_id, cluster_nodes in overlapping_clusters["Core Node Clusters"].items():
        # Core Node Clusters mean the original IKC clusters
        subgraph = G.subgraph(cluster_nodes)
        ls = subgraph.size()
        ds = sum([node_info["total_degree"][node] for node in cluster_nodes])
        modularity = ls / l - (ds / (2 * l)) ** 2
        if modularity <= 0:
            logging.info(
                "Original Cluster " + str(cluster_id) + " has non-positive modularity!"
            )
        nodes_added = 0
        overlapping_clusters["mcd"][cluster_id] = get_mcd(G, cluster_nodes)
        added_candidates = set()  # Baqiao added
        for node in candidates:
            # Baqiao rewrote most this block.
            if (
                node not in cluster_nodes
                and len(node_info["neighbors"][node].intersection(cluster_nodes))
                >= overlapping_clusters[inclusion_criterion][
                    cluster_id
                ]  # cluster_nodes from IKC
            ):
                overlapping_init_core_count = len(
                    node_info["neighbors"][node].intersection(cluster_nodes)
                )  # number of neighbors to node within the original IKC cluster
                overlapping_added_core_count = len(
                    node_info["neighbors"][node].intersection(added_candidates)
                )  # number of neighbors to node within the growing cluster that have been added (i.e., not in the original IKC cluster)
                ls_prime = (
                    ls + overlapping_init_core_count + overlapping_added_core_count
                )  # Total number of edges in the growing cluster
                ds_prime = (
                    ds + node_info["total_degree"][node]
                )  # total degree of the nodes in the growing cluster
                new_modularity = recalculate_modularity(
                    ls_prime, ds_prime, l
                )  # the new modularity of the growing cluster
                if new_modularity > 0:
                    # this block indicates that the candidate node is accepted
                    # because the node has passed both tests
                    overlapping_clusters["Full Clusters"][cluster_id].add(
                        node
                    )  # adding the node to the overlapping clusters
                    overlapping_node_to_cluster_id[node].add(
                        cluster_id
                    )  # maintain the data structure of pointing a node to its set of cluster ids
                    added_candidates.add(node)  # Baqiao added
                    ls = ls_prime  # remember the updated ls for the next iteration
                    ds = ds_prime  # remember the updated ds for the next iteration
                    nodes_added += 1
                    overlapping_clusters["modularity"][cluster_id] = new_modularity

    return overlapping_clusters, overlapping_node_to_cluster_id


def recalculate_modularity(ls, ds, l):
    """
    Iteratively updates the modularity given a new node assigned to the cluster
    Modified from Eleanor Wedell's modularity code in eleanor/code/parsing_clusters_strict.py

    Input:
    ls int
    ds int
    l int
    modularity float - current modularity of cluster
    node_info {dict} - info on all nodes in network
    node int - node id of node added to cluster

    Output:
    modularity float - updated modularity of cluster
    """
    modularity = ls / l - (ds / (2 * l)) ** 2
    return modularity


def get_mcd(G, cluster):
    """
    Returns the mcd of a given cluster

    Input:
    G networkx.Graph - networkx representation of input network
    cluster [list] - list of node ids representing a single cluster

    Output:
    mcd int  - minimum core degree
    """
    G_prime = G.subgraph(cluster)
    a = [val for (node, val) in G_prime.degree()]
    if not a:
        return 0
    mcd = min(a)
    return mcd


def overlapping_clusters_to_output(io, overlapping_clusters):
    """
    Method to save overlapping cluster to a clustering file

    Input:
    io IO - IO object for output
    overlapping_clusters {dict} - dictionary of overlapping clusters by cluster id

    Output:
    None
    """

    for cluster_id, node_set in overlapping_clusters["Full Clusters"].items():
        for node in node_set:
            io.write(cluster_id + " " + node + "\n")


def network_to_dict(network_file):
    """
    Method to extract network data from network file

    Info Collected:
    neighbors - set of neighbors of a given node
    total_degree - total degree of a given node
    indegree - citations given to a specific node
    outdegree - references for a given node
    candidate_type - type of node {original, candidate}

    Input:
    network_file str - file path to network file to parse

    Output:
    node_info {dict} dictionary of info on each node in the network
    """
    node_info = defaultdict()
    node_info["neighbors"] = defaultdict(set)
    node_info["citing_neighbors"] = defaultdict(set)
    node_info["total_degree"] = defaultdict(int)
    node_info["indegree"] = defaultdict(int)
    node_info["outdegree"] = defaultdict(int)
    node_info["candidate_type"] = defaultdict(str)

    network_file_reader = open(network_file, "r")

    line = network_file_reader.readline()
    if "," in line:
        ikc_output_mode = True
        logging.info("Comma separated IKC clustering format detected as input")
    else:
        ikc_output_mode = False
    while line != "":
        if ikc_output_mode:
            v2, v1, _, _ = line.split(",")
            v2, v1 = v2.strip(), v1.strip()
        else:
            v1, v2 = line.split()
            v1, v2 = v1.strip(), v2.strip()
        node_info["total_degree"][v1] += 1
        node_info["outdegree"][v1] += 1
        node_info["indegree"][v1] += 0
        node_info["neighbors"][v1].add(v2)
        node_info["candidate_type"][v1] = "original"

        node_info["total_degree"][v2] += 1
        node_info["outdegree"][v2] += 0
        node_info["indegree"][v2] += 1
        node_info["neighbors"][v2].add(v1)
        node_info["candidate_type"][v2] = "original"
        node_info["citing_neighbors"][v2].add(v1)

        line = network_file_reader.readline()

    return node_info


def clustering_to_dict(clustering, min_k_core):
    """
    Method to extract cluster data from a clustering file

    Info Collected:
    Full Clusters - set of all nodes in a given cluster
    Core Node Clusters - set of all core nodes in a given cluster
    Periphery Node Clusters - set of all periphery nodes in a given cluster
    mcd - minimum core degree of a given cluster
    k - k value of a given cluster

    Input:
    clustering str - file path to clustering file to parse
    min_k_core int - minimum k value to parse

    Output:
    clusters {dict} - dictionary of clusters by cluster id
    node_to_cluster_id {dict} - dictionary of node ids mapped to their disjoint cluster id
    """
    clusters = defaultdict()
    clusters["Full Clusters"] = defaultdict(set)
    clusters["Core Node Clusters"] = defaultdict(set)
    clusters["mcd"] = defaultdict(int)
    clusters["k"] = defaultdict(int)
    clusters["modularity"] = defaultdict(float)
    node_to_cluster_id = defaultdict(set)
    with open(clustering, "r") as fh:
        for l in fh:
            cid_raw, nid_raw = l.strip().split()
            cid = cid_raw
            nid = nid_raw
            node_to_cluster_id[nid].add(cid)
            clusters["Full Clusters"][cid].add(nid)
            clusters["Core Node Clusters"][cid].add(nid)
            clusters["mcd"][cid] = min_k_core
            clusters["k"][cid] = min_k_core
            node_to_cluster_id[nid].add(cid)
    return clusters, node_to_cluster_id


if __name__ == "__main__":
    main()
