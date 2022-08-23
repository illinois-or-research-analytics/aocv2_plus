from collections import Counter, defaultdict
import networkx as nx
from overlapping_kmp_pipeline import network_to_dict
from cluster_id_matching import clustering_to_dict
import os
import numpy as np


def main():
    """
    original_clusters = defaultdict()
    k_clusters = defaultdict()
    mcd_clusters = defaultdict()

    for k in [10, 20, 30, 40, 50]:
      clusters = clustering_to_dict('/shared/aj_manuscript_data/experiment_0/IKC_' + str(k) + '_realignment.clustering')
      original_clusters[k] = set.union(*clusters.values())


    for file in os.listdir('/shared/aj_manuscript_data/experiment_1'):
      if file.endswith('.clustering'):
        k_clusters[file.split('.')[0]] = clustering_to_dict('/shared/aj_manuscript_data/experiment_1/'+ str(file))

    print('Finished EXP 1')

    for file in os.listdir('/shared/aj_manuscript_data/experiment_2'):
      if file.endswith('.clustering'):
        k_clusters[file.split('.')[0]] = clustering_to_dict('/shared/aj_manuscript_data/experiment_2/' + str(file))

    print('Finished EXP 2')

    for file in os.listdir('/shared/aj_manuscript_data/experiment_3'):
      if file.endswith('.clustering'):
        mcd_clusters[file.split('.')[0]] = clustering_to_dict('/shared/aj_manuscript_data/experiment_3/' + str(file))

    print('Finished EXP 3')
    """
    network_file = "/srv/local/shared/external/dbid/george/exosome_dimensions_wedell_retraction-depleted_jc250-corrected_no_header.tsv"
    node_info = network_to_dict(network_file)

    print("Finished Parsing Network")

    degree_map = Counter(node_info["indegree"])
    degrees = list(map(int, degree_map.values()))
    degree_cutoff = np.percentile(degrees, 99)
    print(degrees[:100])
    return


if __name__ == "__main__":
    main()
