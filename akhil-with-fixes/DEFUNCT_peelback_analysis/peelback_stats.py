import matplotlib.pyplot as plt
import networkit as nk
import numpy as np
import pandas as pd
import seaborn as sns
import subprocess

def wccount(filename):
  out = subprocess.Popen(['wc', '-l', filename],
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT
  ).communicate()[0]
  return int(out.partition(b' ')[0])


def modularity():
    return 0

def get_node_edge_count(filepath_node, filepath_edge):
  f = open(filepath_edge, 'r')
  out_degree = {}
  in_degree = {}
  unique_nodes = set()

  for line in f.readlines():
    nodes = line.split("\t")
    unique_nodes.update([nodes[0].strip(), nodes[1].strip()])

    if nodes[0] not in out_degree:
      out_degree[nodes[0].strip()] = 1
    else:
      out_degree[nodes[0].strip()] += 1
  
    if nodes[1] not in in_degree:
        in_degree[nodes[1].strip()] = 1
    else:
        in_degree[nodes[1].strip()] += 1
  return wccount(filepath_node), wccount(filepath_edge), out_degree, in_degree, unique_nodes



def save_edge_node(node_count, edge_count, year, pd_dict):
  pd_dict['count'].append(node_count)
  pd_dict['type'].append('node')
  pd_dict['year'].append(year)
  pd_dict['count'].append(edge_count)
  pd_dict['type'].append('edge')
  pd_dict['year'].append(year)

def save_degree(out_degree, in_degree, unique_nodes):
  degree_dct = {'node_id':[], 'out_degree':[], 'in_degree':[]}
  for degree in unique_nodes:
    degree_dct['node_id'].append(degree)
    if degree not in out_degree:
      degree_dct['out_degree'].append(0)
    else:
      degree_dct['out_degree'].append(out_degree[degree])
    
    if degree not in in_degree:
      degree_dct['in_degree'].append(0)
    else:
      degree_dct['in_degree'].append(in_degree[degree])
  
  return degree_dct


def generate_full_network_clusters(input_filepath, output_filepath):
  f = open(input_filepath, 'r')
  f_out = open(output_filepath, 'r+')
  unique_nodes = set()
  for line in f.readlines():
    nodes = line.split("\t")
    unique_nodes.update(nodes)
  
  f.close()
  for n in unique_nodes:
    f_out.write("0\t" + str(n) + '\n')

  for line in f_out:
    if not line.isspace():
      f_out.write(line)
  f_out.close()


def generate_peelback_clusters(network_file, cluster_file, out_file):
  network = open(network_file, 'r')
  cluster = open(cluster_file, 'r')
  out = open(out_file, 'w')
  unique_nodes = set()

  for line in network.readlines():
    n1, n2 = line.split('\t')
    n1, n2 = n1.strip(), n2.strip()
    unique_nodes.add(n1)
    unique_nodes.add(n2)
  network.close()
  
  for line in cluster.readlines():
    n1, n2 = line.split(' ')
    n1, n2 = n1.strip(), n2.strip()
    if n2 in unique_nodes:
      out.write(line)

  cluster.close()




def main():
  for i in range(0,31,5):
    network_file = './peelback_tsv/peelback' + str(1990+i) + '.tsv'
    cluster_file = '../kmp_files/final_kmp_converted.clustering'
    out_file = './peelback_clusters/peelback' + str(1990+i) + '.clusters'
    generate_peelback_clusters(network_file, cluster_file, out_file)
    print('Finished ' + str(1990 + i) + ' peelback cluster')
  '''  
  print(pd_dict)
  df = pd.DataFrame.from_dict(pd_dict)
  sns.barplot(x="year", y="count", hue="type", data=df)
  plt.savefig('./plots/node_edge_plot.png')
  '''
if __name__ == "__main__":
  main()
