from overlapping_kmp_pipeline import network_to_dict
from collections import Counter, defaultdict
from cluster_id_matching import clustering_to_dict
import numpy as np

def main():
  f = open('candidate_file.txt', 'w')
  node_info =  network_to_dict('/srv/local/shared/external/dbid/george/exosome_dimensions_wedell_retraction-depleted_jc250-corrected_no_header.tsv')
  degree_map = Counter(node_info['total_degree'])
  degrees = list(map(int, degree_map.values()))
  degree_cutoff = np.percentile(degrees, 99)
  candidates = set([k for k,v in degree_map.items() if v >= degree_cutoff])
  print(len(candidates))
  clusters = clustering_to_dict('/shared/aj_manuscript_data/experiment_0/IKC_10_realignment.clustering')
  cluster_nodes = set().union(*clusters.values())
  print(len(cluster_nodes))
  print(len(candidates))
  candidates = candidates - cluster_nodes
  print(len(candidates))
  for node in candidates:
    f.write(str(node) + '\n')

  f.close()
  return




if __name__ == '__main__':
  main()

