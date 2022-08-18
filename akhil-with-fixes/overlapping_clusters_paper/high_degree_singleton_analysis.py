from collections import Counter, defaultdict
import networkx as nx
from overlapping_kmp_pipeline import network_to_dict
from cluster_id_matching import clustering_to_dict
import os 
import numpy as np


def main():
  '''
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
  '''
  network_file = '/srv/local/shared/external/dbid/george/exosome_dimensions_wedell_retraction-depleted_jc250-corrected_no_header.tsv'
  node_info = network_to_dict(network_file)
  
  print('Finished Parsing Network')

  degree_map = Counter(node_info['indegree'])
  degrees = list(map(int, degree_map.values()))
  degree_cutoff = np.percentile(degrees, 99)
  print(degrees[:100])
  return
  candidates = set([k for k,v in degree_map.items() if v >= degree_cutoff])
  
  for name, clusters in k_clusters.items():
    total_clusters = set()
    for c_id in clusters.keys():
      if len(clusters[c_id]) > 0:
        total_clusters.update(clusters[c_id])
    
    high_degree_singletons = candidates - total_clusters
    high_degree_singletons_sorted = sorted(high_degree_singletons, key=lambda item: node_info['indegree'][item], reverse=True)
    #low_degree_inclusions = candidates.union(total_clusters) - original_clusters[int(name[4:6])]
    #low_degree_inclusions_sorted = sorted(low_degree_inclusions, key=lambda item: node_info['indegree'][item], reverse=False)

    #print(len(low_degree_inclusions))

    writer = open('./k_singleton_data/k_' + str(name) + '_singleton_data.tsv', 'w')
    writer.write('node_id\tindegree\n')
    for node in high_degree_singletons_sorted:
      writer.write(str(node) + '\t' + str(node_info['indegree'][node]) + '\n')
    writer.close()
  
  for name, clusters in mcd_clusters.items():
    total_clusters = set()
    for c_id in clusters.keys():
      if len(clusters[c_id]) > 0:
        total_clusters.update(clusters[c_id])
    
    high_degree_singletons = candidates - total_clusters
    high_degree_singletons_sorted = sorted(high_degree_singletons, key=lambda item: node_info['indegree'][item], reverse=True)
    #low_degree_inclusions = candidates.union(total_clusters) - original_clusters[int(name[4:6])]
    #low_degree_inclusions_sorted = sorted(low_degree_inclusions, key=lambda item: node_info['indegree'][item], reverse=False)
    
    #print(len(low_degree_inclusions))

    writer = open('./mcd_singleton_data/mcd_' + str(name) + '_singleton_data.tsv', 'w')
    writer.write('node_id\tindegree\n')
    for node in high_degree_singletons_sorted:
      writer.write(str(node) + '\t' + str(node_info['indegree'][node]) + '\n')
    writer.close()

  doi_map = defaultdict(str)
  r = open('seed_map.csv', 'r')
  line = r.readline()
  while line != "":
    doi_info = line.split(',')
    doi_map[doi_info[0]] = str(doi_info[1])
    line = r.readline()


  node_ids = []
  degree_map = defaultdict(int)
  for path in ['./k_singleton_data', './mcd_singleton_data']:
    for file in os.listdir(path):
      reader = open(path + '/' + file, 'r')
      line = reader.readline()
      line = reader.readline()
      while line != "":
        info = line.split('\t')
        degree_map[info[0]] = info[1]
        node_ids.append(int(info[0]))
        line = reader.readline()
  
  c = Counter(node_ids)

  writer = open('singleton_counts.tsv', 'w')
  writer.write('rank\tnode_id\tcount\tdegree\tdoi\n')
  rank = 1
  for item in c.most_common():
    writer.write(str(rank) + '\t' + str(item[0]) + '\t' + str(item[1]) + '\t' + str(degree_map[str(item[0])]) + '\t' + doi_map[str(item[0])] + '\n')
    rank += 1
  writer.close()

if __name__ == '__main__':
  main()
