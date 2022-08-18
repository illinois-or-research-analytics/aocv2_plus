#import click
import numpy as np
import networkit as nk
import networkx as nx
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from operator import itemgetter

'''
Generates a dictionary representation of a cluster given a filename

Input:
  filename [string]

Output:
  clusters [dictionary]
'''
def generate_clusters(filename):
  clusters = dict()
  f = open(filename, 'r')
  for line in f.readlines():
      cluster_id, node_id = line.split(' ')
      cluster_id, node_id = int(cluster_id.strip()), int(node_id.strip())
      
      if cluster_id not in clusters:
        clusters[cluster_id] = [node_id]
      else:
        clusters[cluster_id].append(node_id)
  # Remove singleton clusters from dictionary
  clusters = {k:v for k, v in clusters.items() if len(v) > 1}
  return clusters

'''
Generates a networkx graph from a filename given a tab spaced edgelist

Input:
  filename [string]

Output:
  graph [nx.DiGraph]
'''
def generate_graph(filename):
  return nx.read_edgelist(filename, delimiter='\t', create_using=nx.DiGraph)


'''
Generates all degree tuple pairs in a list that have a degree value in the percentile
group specified or greater

Input:
  degree_tup [list of tuples]
  percentile [float] range:0.0-100.0

Output:
  degree_tup [list of tuples]

'''
def generate_percentile_plus(degree_tup, percentile):
  degree_tup = list(sorted(degree_tup, key = lambda x: x[1]))
  return degree_tup[int((percentile/100)*len(degree_tup)):]


'''
Compares the similarity of two lists using set intersection

Input:
  l1 [list]
  l2 [list]

Output:
  percentage [float]
'''
def similarity_percentage(l1, l2):
  total = 0
  '''
  citation_increase = []
  for k1,v1 in l1:
    for k2,v2 in l2:
      if k1 == k2:
        total += 1
        citation_increase.append(v2-v1)
  '''
  l1_set = set(X[0] for X in l1)
  l2_set = set(X[0] for X in l2)
  len_intersection = len(l1_set.intersection(l2_set))
  len_l1 = len(l1_set)
  return float(len_intersection/len_l1), len_intersection, len_l1
  #return float(total/len(l1)), max(citation_increase), min(citation_increase), float(sum(citation_increase)/len(citation_increase))
'''
Generates a list of induced subgraphs given a list of clusters and a base network

Input:
  network [nx.DiGraph]
  clusters [dict]

Output:
  induced_subgraphs [list]
'''
def induce_subgraph(G, clusters):
  induced_subgraphs = []
  print('Inducing Subgraphs')
  for _, v in clusters.items():
    s = nx.DiGraph()
    s.add_edges_from((n, nbr, d) for n, nbrs in G.adj.items() if n in v for nbr, d in nbrs.items() if nbr in v)
    if (nx.classes.function.is_empty(s)) is False:
      print('Adding Real Subgraph of Size' + str(s.number_of_nodes()))
      induced_subgraphs.append(s)
  return induced_subgraphs


#def get_connected(G, G_s):
  

def main():
  for i in range(0, 31, 5):
    
  #connected_dict = {'Cluster_Length':[], 'is_strongly_connected':[], 'is_weakly_connected':[]}
    network_file = './peelback_tsv/peelback' + str(1990+i) + '.tsv'
    filename = './peelback_clusters/peelback' + str(1990+i) + '.clusters'
    G = generate_graph(network_file)
    #clusters = generate_clusters(filename)
    #induced_subgraphs = induce_subgraph(G, clusters)
    print(1990+i, G.number_of_nodes()/G.number_of_edges(), G.number_of_nodes(), G.number_of_edges())
  '''
  for s in induced_subgraphs:
    if (nx.classes.function.is_empty(s)) is not True :
      connected_dict['Cluster_Length'].append(s.number_of_nodes())
      connected_dict['is_strongly_connected'].append(nx.is_strongly_connected(s))
      connected_dict['is_weakly_connected'].append(nx.is_weakly_connected(s))
  df = pd.DataFrame.from_dict(connected_dict)
  plt.figure() 
  sns.scatterplot(x=df['is_strongly_connected'], y=df['is_weakly_connected'], hue=df['Cluster_Length'])
  plt.savefig('../plots/connectivity_correlation_cluster_size' + str(1990+i) + '.png')
  '''

  '''
  original_cluster_file = open(kmp_file, 'r')
  for i in range(0,31, 5):
    f = open(prefix + str(1990+i) + suffix, 'r')
    f_write = open('./peelback_clusters_converted/peelback' + str(1990+i) + '_kmp.clustering', 'w')
    linked_nodes = generate_linked_nodes(f)
    print(len(linked_nodes))
    for line in original_cluster_file.readlines():
      cluster_id, node_id = line.split(' ')
      cluster_id, node_id = cluster_id.strip(), node_id.strip()
    
      if node_id in linked_nodes:
        f_write.write(line)

    print('Done w/ '+str(1990+i))
  '''
  #full_network = '/srv/local/shared/external/Dimensions/deduplicated_datasets_minhyuk2/exosome_1900_2010_sabpq/citing_cited_network.integer.tsv'
  '''
  G, G2 = None, None
  for i in range(0, 30, 5):
    network = "./peelback_tsv/peelback" + str(1990+i) + ".tsv"
    network2 = "./peelback_tsv/peelback" + str(1990+i+5) + ".tsv"
    if i > 0:
       G = G2
       G2 = generate_graph(network2)
    else:
      G = generate_graph(network)
      G2 = generate_graph(network2)
    #outdegree = [X for X in G.out_degree() if X[1] < 1000]
    indegree = G.in_degree()
    #outdegree = [X for X in G.out_degree() if X[1] < 1000]
    indegree2 = G2.in_degree()

    for percentile in np.arange(95, 100, 0.5):
      percentile_list = generate_percentile_plus(indegree, percentile)
      percentile_list2 = generate_percentile_plus(indegree2, percentile)
      print(percentile)
      print(similarity_percentage(percentile_list, percentile_list2))
    
  
  for i in range(30,31,5):
    network = "./peelback_tsv/peelback" + str(1990+i) + ".tsv"
    #cluster_filename = './peelback_clusters_converted/peelback' + str(1990+0) + '_kmp.clustering'
    #main_cluster = 'peelback_clusters_converted/original_kmp_cluster_50_2.clustering'
    G = generate_graph(network)
    #nx.readwrite.gpickle.write_gpickle(G, 'full_network.pickle')
    outdegree = [X[1] for X in G.out_degree() if X[1] < 1000]
    indegree = [X[1] for X in G.in_degree()]
    #print(len(outdegree), len(indegree), max(outdegree), max(indegree))
    #print(max(G.out_degree(),key=itemgetter(1))[0])
    plt.figure()
    sns.histplot(data=outdegree, bins=int(len(outdegree)/10000))
    plt.yscale('symlog')
    plt.xlabel('Out-Degree Value', fontsize=14)
    plt.ylabel('Count', fontsize=16)
    plt.savefig('../plots/outdegree_distribution_' + str(1990+i) + '.png')
    print('Finished out-degree graph')
    plt.figure()
    sns.histplot(data=indegree, bins=int(len(indegree)/10000))
    plt.yscale('symlog')
    plt.xlabel('In-Degree Value', fontsize=14)
    plt.ylabel('Count', fontsize=16)
    plt.savefig('../plots/indegree_distribution_' + str(1990+i) + '.png')
    print('Finished in-degree graph')
    '''
if __name__ == "__main__":
  main()

