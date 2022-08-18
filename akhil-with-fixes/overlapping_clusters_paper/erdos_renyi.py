import networkx as nx
import os
from overlapping_kmp_pipeline import network_to_dict, clustering_to_dict
#from cluster_id_matching import clustering_to_dict
import seaborn as sns
import matplotlib.pyplot as plt


'''
Program to generate and analyze erdos renyi graphs with the same number of vertices and 
edges as the original CEN
'''

def main():
  for k in [10, 20, 30, 40, 50]:
    for file in os.listdir('./erdos_renyi_graphs'):
      print(file, k)
      run_ikc('./erdos_renyi_graphs/' + file, './erdos_renyi_ikc_clustering/', k)
  
  for file in os.listdir('./erdos_renyi_graphs'):
    print(file)
    node_info = network_to_dict('./erdos_renyi_graphs/' + file)
    plt.figure()
    #node_info = network_to_dict('/srv/local/shared/external/dbid/george/exosome_dimensions_wedell_retraction-depleted_jc250-corrected_no_header.tsv')
    plt.hist(node_info['indegree'].values(), density=True, bins=100)
    plt.ylabel('Count Density')
    plt.xlabel('Indegree Value');
    f = file.split('.')[0]
    plt.savefig('./erdos_renyi_indegree_distribution/' + str(f) + '_distribution.png')
  
  for file in os.listdir('./erdos_renyi_ikc_clustering'):
    clusters, _ = clustering_to_dict('./erdos_renyi_ikc_clustering/' + file, 10)
    print(len([k for k in clusters['Full Clusters'].keys() if len(clusters['Full Clusters'][k]) > 1]))



  
def generate_errdos_renyi_graphs():
  for i in range(100):
    #G = nx.gnm_random_graph(n=139,m=920,seed=i,directed=True)
    G = nx.gnm_random_graph(n=13989436,m=92051051,seed=i,directed=True)
    file_path = './erdos_renyi_graphs/erdos_renyi_' + str(i) + '.tsv'
    fp = open(file_path, 'w')
    for u,v in G.edges():
      fp.write(str(u) + '\t' + str(v) + '\n')

    fp.close()
    G.clear()



def run_ikc(network_file, output_path, k):
  filename = network_file.split('/')[2]
  filename = filename.split('.')[0]
  command = 'pipenv run python3 modified_IKC.py -e ' + str(network_file) + ' -o ' + str(output_path) + 'k_'+ str(k) + '_' + str(filename) + '.clustering -k ' + str(k)
  print(command)
  os.system(command)


if __name__ == '__main__':
  main()
