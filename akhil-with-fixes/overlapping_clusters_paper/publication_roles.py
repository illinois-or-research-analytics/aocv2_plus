import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from overlapping_kmp_pipeline import network_to_dict, parse_marker_file, get_mcd
from cluster_id_matching import clustering_to_dict
import networkx as nx
from collections import defaultdict
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

def main():

  #clusters = clustering_to_dict('/shared/gc/experiment_55/equil_IKC_10.clustering_k')
  clusters_original = clustering_to_dict('/shared/aj_manuscript_data/experiment_0/IKC_10_realignment.clustering')
  clusters_k = clustering_to_dict('/shared/aj_manuscript_data/experiment_1/IKC_10_km_totaldegree_1percent.clustering')
  clusters_mcd = clustering_to_dict('/shared/aj_manuscript_data/experiment_3/IKC_10_mcd_totaldegree_1percent.clustering')

  #node_info = network_to_dict('/srv/local/shared/external/dbid/george/exosome_dimensions_wedell_retraction-depleted_jc250-corrected_no_header.tsv')
  network_file = '/srv/local/shared/external/dbid/george/exosome_dimensions_wedell_retraction-depleted_jc250-corrected_no_header.tsv'
  G = nx.read_edgelist(network_file, delimiter='\t')
  for i in range(1, 129):
    #index = str(i)
    print(i)
    index = i
    print(index, len(clusters_original[index]), get_mcd(G, clusters_original[index]), len(clusters_k[index]), get_mcd(G, clusters_k[index]), len(clusters_mcd[index]), get_mcd(G, clusters_mcd[index]))
  return
  df = defaultdict(list)
  cutoffs = defaultdict(int)
  for k,cluster in clusters.items():
    indegree_count = []
    node_value= 0
    for node_id in cluster:
      indegree_count.append(len(node_info['citing_neighbors'][node_id].intersection(cluster)))
    
    cutoff = np.percentile(indegree_count, 90)
    cutoffs[k] = cutoff
  
  print('Finished Cutoffs')

  cluster_nodes = set().union(*clusters.values())
  cluster_count = defaultdict(int)
  tier_count = defaultdict(int)
  node_map = defaultdict(set)
  for k, cluster in clusters.items():
    for node in cluster:
      node_map[node].add(k)
      cluster_count[node] += 1
      if len(node_info['citing_neighbors'][node].intersection(cluster)) >= cutoffs[k]:
        tier_count[node] += 1
  
  print('Finished Data Generation')

  for k in cluster_count.keys():
    for _id in node_map[k]:
      df['cluster_id'].append(_id)
      df['node_id'].append(k)
      df['cluster_count'].append(cluster_count[k])
      df['tier_1_count'].append(tier_count[k])
      df['indegree'].append(node_info['indegree'][k])

  df_pd = pd.DataFrame.from_dict(df)
  print(len(df_pd.index))
  df_pd.to_csv('/shared/aj_manuscript_data/tiervclustercount_singleton_k.csv')
  
  print(df_pd.groupby('cluster_id').describe())
  print()
  print(df_pd.groupby('cluster_count')['tier_1_count'].describe())
  plt.figure()
  plt.xlim((0,25))
  plt.ylim((0,25))
  
  df_pd['indegree'] = pd.cut(df_pd["indegree"], [0, 100, 1000, 10000, 50000, 100000], labels=["<100", "1000", "10000", "50000", ">100000"])
  g = sns.scatterplot(data=df_pd, x='cluster_count', y ='tier_1_count', hue='indegree', palette='CMRmap')
  g.legend(loc='center left', bbox_to_anchor=(1.25, 0.5), ncol=1)
  plt.savefig('./plots/tiervclustercount_singleton_k.png')
  
  return

def num_clusters_found(node, clusters):
  count = 0
  for k,v in clusters.items():
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
  node_value= 0
  for node_id in cluster:
    indegree_count.append(len(node_info['citing_neighbors'][node_id].intersection(cluster)))
    if node_id == node:
      node_value = len(node_info['citing_neighbors'][node_id].intersection(cluster))
  
  cutoff = np.percentile(indegree_count, 90)
  return int(node_value >= cutoff)



if __name__ == '__main__':
  main()
