import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from overlapping_kmp_pipeline import network_to_dict, parse_marker_file, get_mcd
from cluster_id_matching import clustering_to_dict
import networkx as nx
from collections import defaultdict

pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", None)


def main():

    # clusters = clustering_to_dict('/shared/gc/experiment_55/equil_IKC_10.clustering_k')
    clusters_original = clustering_to_dict(
        "/shared/aj_manuscript_data/experiment_0/IKC_10_realignment.clustering"
    )
    clusters_k = clustering_to_dict(
        "/shared/aj_manuscript_data/experiment_1/IKC_10_km_totaldegree_1percent.clustering"
    )
    clusters_mcd = clustering_to_dict(
        "/shared/aj_manuscript_data/experiment_3/IKC_10_mcd_totaldegree_1percent.clustering"
    )

    # node_info = network_to_dict('/srv/local/shared/external/dbid/george/exosome_dimensions_wedell_retraction-depleted_jc250-corrected_no_header.tsv')
    network_file = "/srv/local/shared/external/dbid/george/exosome_dimensions_wedell_retraction-depleted_jc250-corrected_no_header.tsv"
    G = nx.read_edgelist(network_file, delimiter="\t")
    for i in range(1, 129):
        # index = str(i)
        print(i)
        index = i
        print(
            index,
            len(clusters_original[index]),
            get_mcd(G, clusters_original[index]),
            len(clusters_k[index]),
            get_mcd(G, clusters_k[index]),
            len(clusters_mcd[index]),
            get_mcd(G, clusters_mcd[index]),
        )
    return


def num_clusters_found(node, clusters):
    count = 0
    for k, v in clusters.items():
        if node in v:
            count += 1
    return count


def tier_1_clusters_found(node, clusters, node_info):
    count = 0
    for key, cluster in clusters.items():
        count += tier_1_count(node, cluster, node_info)

    return count


def tier_1_count(node, cluster, node_info):
    indegree_count = []
    node_value = 0
    for node_id in cluster:
        indegree_count.append(
            len(node_info["citing_neighbors"][node_id].intersection(cluster))
        )
        if node_id == node:
            node_value = len(
                node_info["citing_neighbors"][node_id].intersection(cluster)
            )

    cutoff = np.percentile(indegree_count, 90)
    return int(node_value >= cutoff)


if __name__ == "__main__":
    main()
