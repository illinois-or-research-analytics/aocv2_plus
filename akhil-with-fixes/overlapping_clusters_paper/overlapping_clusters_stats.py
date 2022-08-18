from collections import Counter, defaultdict
import overlapping_kmp_pipeline as okmp
import pandas as pd
import statistics


'''
Defunct method to generate basic cluster and size distribution statistics given a clustering

Input:
  clusters {dict} - dictionary of clusters by cluster id

Output:
  num_clusters int - number of clusters in clustering
  size_distribution pd.Series - size distribution and statistcs of all clusters in clustering
'''
def basic_cluster_stats(clusters):
  num_clusters = len(clusters['Full Clusters'].keys())
  size_distribution = pd.Series([len(nodes) for _, nodes in clusters['Full Clusters'].items()]).describe()
  return num_clusters, size_distribution

'''
Defunct method used to analyze neighbor coverage of overlapping clustering algorithm

Input:
  overlapping_clusters {dict} - dictionary of overlapping clusters by cluster id
  original_clusters {dict} - dictionary of original disjoint clusters by cluster id
  candidates [list] - list of candidate node ids
  node_info {dict} - dictionary of info on nodes in network
  overlapping_node_to_cluster_id {dict} - dictionary mapping node ids to all clusters they are in
  candidate_criterion str - criterion to judge acceptability of node as a candidate

Output:
  neighbor_info {dict} - dictionary of neighbor information
  num_non_added_nodes int - number of nodes not added to any cluster
'''
def neighbor_analysis(overlapping_clusters, original_clusters, candidates, node_info, overlapping_node_to_cluster_id, original_node_to_cluster_id, candidate_criterion):
  num_non_added_nodes = 0
  neighbor_info = defaultdict(list)

  for index, candidate in enumerate(candidates):
    all_cluster_ids = overlapping_node_to_cluster_id[candidate]
    original_cluster_id = None
    if candidate in original_node_to_cluster_id:
      original_cluster_id = list(original_node_to_cluster_id[candidate])[0]
    non_singleton_clusters = []

    for c_id in all_cluster_ids:
      if len(overlapping_clusters['Full Clusters'][c_id]) > 0:
        non_singleton_clusters.append(overlapping_clusters['Full Clusters'][c_id])
    if len(non_singleton_clusters) > 0 and original_cluster_id is not None:
      proportion_neighbors_original = proportion_neighbors_original_clusters(original_clusters['Full Clusters'][original_cluster_id], candidate, node_info['neighbors'][candidate])
      proportion_neighbors_overlapping = proportion_neighbors_all_clusters(non_singleton_clusters, candidate, node_info['neighbors'][candidate])
      neighbor_info['proportion_neighbors_overlapping_clusters'].append(proportion_neighbors_overlapping)
      neighbor_info['proportion_neighbors_original_clusters'].append(proportion_neighbors_original)
      neighbor_info['proportion_neighbor_coverage_increase'].append((proportion_neighbors_overlapping-proportion_neighbors_original)/proportion_neighbors_original)
      neighbor_info[candidate_criterion].append(node_info[candidate_criterion][candidate])
      neighbor_info['degree_rank'].append(index)
    else:
      num_non_added_nodes += 1

  return neighbor_info, num_non_added_nodes

'''
Helper function to calculate proportion of neighbors of a node in its original cluster

Input:
  original_cluster  [list] - list of all nodes in original cluster with given node
  node int - node id of given node
  neighbors [list] - list of all neighbors of given node

Output:
  proportion int - proportion of all neighbors of given node in original cluster
'''
def proportion_neighbors_original_clusters(original_cluster, node, neighbors):
  return float(len(neighbors.intersection(original_cluster))/len(neighbors))

'''
Helper function to calculate proportion of neighbors of a node in its overlapping clusters

Input:
  all_clusters  [list] - list of all nodes in all overlapping clusters that contain the given node
  node int - node id of given node
  neighbors [list] - list of all neighbors of given node

Output:
  proportion int - proportion of all neighbors of given node in all its overlapping clusters
'''
def proportion_neighbors_all_clusters(all_clusters, node, neighbors):
  total_clusters = set().union(*all_clusters)
  return float(len(neighbors.intersection(total_clusters))/len(neighbors))

'''
Stores and outputs node placement frequency in multiple clusters

Input:
  clusters {dict} - dictionary of clusters by cluster id
  candidates [list] -  list of candidate node ids
  mapping {dict} - mapping of node id to doi
  node_placement_file str - file path to store node placement data

Output:
  placement_frequency [list] - list of candidates by placement into overlapping clusters frequency
'''
def node_placement_frequency(clusters, candidates, mapping, node_placement_file):
  writer = open(node_placement_file, 'w')
  
  placement_frequency = defaultdict(int)
  for node in candidates:
    placement_frequency[node] = 0
  
  for c_id, node_set in clusters['Full Clusters'].items():
    for node in node_set:
      if node in candidates:
        placement_frequency[node] += 1
  
  placement_frequency_sorted = dict(sorted(placement_frequency.items(), key=lambda item: item[1], reverse=True))
  
  for k,v in placement_frequency_sorted.items():
    writer.write(str(mapping[k]) +  ': ' + str(v) + '\n')

  writer.close()
  return [placement_frequency[n_id] for n_id in placement_frequency.keys() if n_id in candidates]
  
'''
Helper function to calculate the number of edges added into clusters when a node is added
to a specific cluster

Input:
  node int - node id of given node
  cluster [list] - list of node ids in the given cluster that the given node is being added to
  node_info {dict} - dictionary of info on nodes in network

Output:
  edges_added int - new edges added inside clusters when the given node is added to the cluster
'''
def edges_added(node, cluster, node_info):
  return len(node_info['neighbors'][node].intersection(cluster))

'''
Tier analysis similar to the one outlined in Chandrashekaran et al.

Input:
  clusters {dict} - dictionary of clulsters by cluster id
  node_info {dict} - dictionary of info on nodes in network
  candidates [list] - list of candidate node ids

Output:
  cluster_tier_info {dict} - tier info for nodes grouped by clusters included in
  aggregate_node_tier_info {dict} - aggregate info for node tiers by tier
'''
def tier_analysis(clusters, node_info, candidates):
  cluster_tier_info = defaultdict(list)
  node_tier_info = defaultdict(list)
  for cluster_id in clusters['Full Clusters']:
    citations = [(n, node_info['indegree'][n]) for n in clusters['Full Clusters'][cluster_id]]
    citations.sort(reverse=True, key = lambda n: n[1])
    threshold_count = int(0.1*len(clusters['Full Clusters'][cluster_id]))
    threshold_value = citations[threshold_count][1]
    greater_than_threshold_count = 0
    for node in clusters['Full Clusters'][cluster_id]:
      if node_info['indegree'][node] >= threshold_value:
        greater_than_threshold_count += 1
    if greater_than_threshold_count > threshold_count:
      threshold_value += 1

    node_tiers = [0, 0, 0]
    for node in clusters['Full Clusters'][cluster_id].intersection(set(candidates)):
        if node_info['indegree'][node] >= threshold_value:
          node_tier_info[node].append(1)
          node_tiers[0] += 1
        elif node_info['indegree'][node] == 0:
          #print(node_info['total_degree'][node], node_info['outdegree'][node])
          node_tier_info[node].append(3)
          node_tiers[2] += 1
        else:
          node_tier_info[node].append(2)
          node_tiers[1] += 1


    cluster_tier_info['Tier 1'].append(node_tiers[0])
    cluster_tier_info['Tier 2'].append(node_tiers[1])
    cluster_tier_info['Tier 3'].append(node_tiers[2])
    cluster_tier_info['Cluster Size'].append(len(clusters['Full Clusters'][cluster_id]))

  aggregate_node_tier_info = defaultdict(list)

  for node, tiers in node_tier_info.items():
    tier_frequency = Counter(tiers)
    aggregate_node_tier_info['node_id'].append(node)
    aggregate_node_tier_info['degree'].append(node_info['total_degree'][node])
    aggregate_node_tier_info['tier_1_count'].append(tier_frequency[1])
    aggregate_node_tier_info['tier_2_count'].append(tier_frequency[2])
    aggregate_node_tier_info['tier_3_count'].append(tier_frequency[3])
    aggregate_node_tier_info['average_tier'].append(sum(tiers)/len(tiers))

    aggregate_node_tier_info['node_type'].append(node_info['candidate_type'][node])
  
  #print('Tier 1 Count:', sum(aggregate_node_tier_info['tier_1_count']), 
        #'Tier 2 Count:', sum(aggregate_node_tier_info['tier_2_count']), 
        #'Tier 3 Count:', sum(aggregate_node_tier_info['tier_3_count']))

  return cluster_tier_info, aggregate_node_tier_info

'''
Method to generate statistics on the intersection points across overlapping clusters

Input:
  output_path str - output file path to save intersection data to
  clusters {dict} - dictionary of clusters by cluster id

Output:
  intersection_stats {dict} - statistics of overlapping clusters intersections
'''
def cluster_intersection_analysis(output_path, clusters):
  intersection_stats = defaultdict(list)

  for index1, c_id in enumerate(clusters['Full Clusters']):
    for index2, c_id2 in enumerate(clusters['Full Clusters']):
      if index1 < index2:
        c1 = clusters['Full Clusters'][c_id]
        c2 = clusters['Full Clusters'][c_id2]
        intersect = c1.intersection(c2)
        intersection_stats['c1'].append(c_id)
        intersection_stats['c2'].append(c_id2)
        intersection_stats['c1_size'].append(len(c1))
        intersection_stats['c2_size'].append(len(c2))

        intersection_stats['intersection_size'].append(len(intersect))
        intersection_stats['jaccard_similarity'].append(jaccard_similarity(c1, c2))
        intersection_stats['f1_score'].append(f1_score(c1, c2))

  df_intersection = pd.DataFrame.from_dict(intersection_stats)
  df_intersection.to_csv(output_path)
        
  #print('Min Intersection Size', min(intersection_stats['intersection_size']))
  #print('Median Intersection Size', statistics.median(intersection_stats['intersection_size']))
  #print('Max Intersection Size', max(intersection_stats['intersection_size']))
  return intersection_stats

'''
Helper function to compute the jaccard similarity of two clusters

Input:
  c1 [list] - list of node ids in cluster 1
  c2 [list] - list of node ids in cluster 2

Output:
  jaccard_similarity float - jaccard similarity calculated for c1 and c2
'''
def jaccard_similarity(c1, c2):
  return float(len(c1.intersection(c2)) / len(c1.union(c2)))

'''
Helper function to compute the f1 score of two clusters

Input:
  c1 [list] - list of node ids in cluster 1
  c2 [list] - list of node ids in cluster 2

Output:
  jaccard_similarity float - jaccard similarity calculated for c1 and c2
'''
def f1_score(c1, c2):
  return float(2*len(c1.intersection(c2)) / (len(c1) + len(c2)))

'''
Method to analyze clusters both before and after the overlapping step

Input:
  G networkx.Graph - networkx representation of input graph
  output_path str - output file path to save cluster analysis to
  clusters {dict} - dictionary of clusters by cluster id
  node_info {dict} - dictionary of info on nodes in network
  is_overlapping bool - if clusters analyzed are overlapping or not

Output:
  cluster_stats {dict} - dictionary of cluster stats including cluster id, cluster size, modularity, and average total degree of 95th percentile of nodes in cluster
'''
def cluster_analysis(G, output_path, clusters, node_info, is_overlapping=False):
  cluster_stats = defaultdict(list)

  for index, c_id in enumerate(clusters['Full Clusters']):
        c1 = clusters['Full Clusters'][c_id]
        cluster_stats['cluster_id'].append(c_id)
        cluster_stats['cluster_size'].append(len(c1))
        cluster_stats['modularity'].append(clusters['modularity'][c_id])
        cluster_stats['average_95_percentile_total_degree'].append(average_top_n_total_degree(c1, node_info, 95)) 
        cluster_stats['mcd'].append(okmp.get_mcd(G, c1))
        if is_overlapping:
          cluster_stats['is_overlapping'].append(True)
        else:
          cluster_stats['is_overlapping'].append(False)

  df_clusters = pd.DataFrame.from_dict(cluster_stats)
  df_clusters.to_csv(output_path)
        
  return cluster_stats


'''
Method to calculate the average total degree of the top n percent nodes in a given cluster

Input:
  cluster [list] - list of node ids in a given cluster
  node_info {dict} - dictionary of info on nodes in network
  n int - top n percent nodes in cluster to consider

Output:
  average_top_n_total_degree float - average total degree of top n percent nodes in cluster
'''
def average_top_n_total_degree(cluster, node_info, n):
  c_sorted = sorted(cluster, key=lambda x:node_info['total_degree'][x], reverse=True)
  c_sorted_n = c_sorted[:int(((100-n)/100)*len(c_sorted))+1]
  c_sorted_n_degree = [node_info['total_degree'][c] for c in c_sorted_n]
  return sum(c_sorted_n_degree)/len(c_sorted_n_degree)


