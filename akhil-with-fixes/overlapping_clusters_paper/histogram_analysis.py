from collections import Counter, defaultdict
import networkx as nx
from overlapping_kmp_pipeline import network_to_dict, parse_marker_file
from cluster_id_matching import clustering_to_dict
import os 
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

def main():
  candidates, _ = parse_marker_file('markers.csv')
  
  original_candidates, _ = parse_marker_file('original_markers.csv')
  candidates.append('1002157')
  candidates = set(candidates).union(original_candidates)
  print(len(candidates))
  
  original_clusters = defaultdict()
  k_clusters = defaultdict()
  mcd_clusters = defaultdict()

  for k in [10, 20, 30, 40, 50]:
    clusters = clustering_to_dict('/shared/aj_manuscript_data/experiment_0/IKC_' + str(k) + '_realignment.clustering')
    original_clusters[k] = set.union(*clusters.values())

  
  for file in os.listdir('/shared/aj_manuscript_data/experiment_1'):
    if file.endswith('.clustering'):
      clusters = clustering_to_dict('/shared/aj_manuscript_data/experiment_1/'+ str(file))
      k_clusters[file.split('.')[0]] = set.union(*clusters.values())


  print('Finished EXP 1')
  
  for file in os.listdir('/shared/aj_manuscript_data/experiment_2'):
    if file.endswith('.clustering'):
      clusters = clustering_to_dict('/shared/aj_manuscript_data/experiment_2/'+ str(file))
      k_clusters[file.split('.')[0]] = set.union(*clusters.values())
  
  print('Finished EXP 2')

  for file in os.listdir('/shared/aj_manuscript_data/experiment_3'):
    if file.endswith('.clustering'):
      clusters = clustering_to_dict('/shared/aj_manuscript_data/experiment_3/'+ str(file))
      mcd_clusters[file.split('.')[0]] = set.union(*clusters.values())

  print('Finished EXP 3')
  
  network_file = '/srv/local/shared/external/dbid/george/exosome_dimensions_wedell_retraction-depleted_jc250-corrected_no_header.tsv'
  node_info = network_to_dict(network_file)
  
  print(len(node_info['indegree'].keys()))
  print(sum(node_info['indegree'].values()))

  return
  print('Finished Parsing Network')

  degree_map = Counter(node_info['indegree'])
  degrees = list(map(int, degree_map.values()))
  degree_cutoff = np.percentile(degrees, 99)
  #candidates = set([k for k,v in degree_map.items() if v >= degree_cutoff])
  candidates, _ = parse_marker_file('markers.csv')
  
  original_candidates, _ = parse_marker_file('original_markers.csv')
  candidates.append('1002157')
  candidates = set(candidates).union(original_candidates)
  print(candidates)

  info = defaultdict()
  #Original Disjoint Clusters
  for k in [10, 20, 30, 40, 50]:
    info[k] = defaultdict(list)
    for node in candidates:
      num_nodes_original_cluster = len(original_clusters[k].intersection(node_info['citing_neighbors'][node]))
      total_nodes_network = len(node_info['citing_neighbors'][node])
      if total_nodes_network > 0:
        info[k]['k'].append(k)
        info[k]['node'].append(int(node))
        info[k]['indegree'].append(node_info['indegree'][node])
        info[k]['proportion'].append(float(num_nodes_original_cluster/total_nodes_network))
        info[k]['name'].append('original')
    print('Finished ', str(k))
  
  #Exp 1
  for file in os.listdir('/shared/aj_manuscript_data/experiment_1'):
    if file.endswith('.clustering'):
      new_k = int(file[4:6])
      for node in candidates:
        num_nodes_original_cluster = len(k_clusters[file.split('.')[0]].intersection(node_info['citing_neighbors'][node]))
        total_nodes_network = len(node_info['citing_neighbors'][node])
        if total_nodes_network > 0:
          info[new_k]['k'].append(new_k)
          info[new_k]['node'].append(int(node))
          info[new_k]['indegree'].append(node_info['indegree'][node])
          info[new_k]['proportion'].append(float(num_nodes_original_cluster/total_nodes_network))
          info[new_k]['name'].append('k_' + file.split('.')[0])
       
  #Exp 2
  for file in os.listdir('/shared/aj_manuscript_data/experiment_2'):
    if file.endswith('.clustering'):
      new_k = int(file[4:6])
      for node in candidates:
        num_nodes_original_cluster = len(k_clusters[file.split('.')[0]].intersection(node_info['citing_neighbors'][node]))
        total_nodes_network = len(node_info['citing_neighbors'][node])
        if total_nodes_network > 0:
          info[new_k]['k'].append(new_k)
          info[new_k]['node'].append(int(node))
          info[new_k]['indegree'].append(node_info['indegree'][node])
          info[new_k]['proportion'].append(float(num_nodes_original_cluster/total_nodes_network))
          info[new_k]['name'].append('k_' + file.split('.')[0])
  '''
  #Exp 3
  for file in os.listdir('/shared/aj_manuscript_data/experiment_3'):
    if file.endswith('.clustering'):
      new_k = int(file[4:6])
      for node in candidates:
        num_nodes_original_cluster = len(mcd_clusters[file.split('.')[0]].intersection(node_info['citing_neighbors'][node]))
        total_nodes_network = len(node_info['citing_neighbors'][node])
        info[new_k]['k'].append(new_k)
        info[new_k]['node'].append(int(node))
        info[new_k]['indegree'].append(node_info['indegree'][node])
        info[new_k]['proportion'].append(float(num_nodes_original_cluster/total_nodes_network))
        info[new_k]['name'].append('mcd_' + file.split('.')[0])
  '''

  info_df = pd.DataFrame.from_dict(info)
  for k in [10, 20, 30, 40, 50]:
    plt.figure()
    
    plt.figure(figsize=(20,20))
    ax = sns.boxplot(data=info_df[k], x='name', y='proportion', orient='v', width=0.5)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=10)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=40, ha="right")
    plt.savefig('./plots/k_boxplot_' + str(k) + '.png')

if __name__ == '__main__':
  main()

