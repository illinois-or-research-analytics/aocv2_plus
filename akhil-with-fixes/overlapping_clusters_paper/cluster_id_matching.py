import pandas as pd
import numpy as np
import os
import overlapping_kmp_pipeline as okmp
from collections import defaultdict

def main():
  files = [40]
  original = []
  overlapping = []
  
  for k in files:
    overlapping2 = './experiment_2/IKC_' + str(k) + '_km_indegree_1percent.clustering'
    overlapping = './experiment_1/IKC_' + str(k) + '_km_totaldegree_1percent.clustering'
    original = './experiment_0/IKC_' + str(k) + '_realignment.clustering'
    
    clusters = clustering_to_dict(original)
    overlapping_clusters = clustering_to_dict(overlapping)
    overlapping_clusters2 = clustering_to_dict(overlapping2)

    for key_ in clusters.keys():
      if len(clusters[key_]) > 1:
        print(key_, len(clusters[key_]), len(overlapping_clusters[key_]), len(overlapping_clusters2[key_]))
  return 0


def clustering_to_dict(clustering):
  clusters = defaultdict(set)

  clustering_reader = open(clustering, 'r')
  line = clustering_reader.readline()
  
  while line != "":
    node_cluster_info = line.split(' ')
    clusters[int(node_cluster_info[0])].add(node_cluster_info[1].strip())
    line = clustering_reader.readline()

  return clusters

def clustering_to_dict2(clustering):
  clusters = defaultdict(set)

  clustering_reader = open(clustering, 'r')
  line = clustering_reader.readline()
  
  while line != "":
    node_cluster_info = line.split(',')
    clusters[int(node_cluster_info[1])].add(node_cluster_info[0].strip())
    line = clustering_reader.readline()

  return clusters



if __name__ == '__main__':
  main()
